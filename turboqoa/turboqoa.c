#include "turboqoa.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>
#include <unistd.h>

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

// QOA constants
#define QOA_MAGIC "qoaf"
#define QOA_FILE_HEADER_SIZE 8
#define QOA_FRAME_HEADER_SIZE 8
#define QOA_LMS_STATE_SIZE 16
#define QOA_SLICE_SIZE 8
#define QOA_SAMPLES_PER_SLICE 20
#define QOA_SLICE_PER_CHANNEL_PER_FRAME 256


static void turboqoa_work_buffer_init(struct TurboQOAWorkBuffer *work_buffer)
{
    work_buffer->data = NULL;
    work_buffer->capacity = 0;
    work_buffer->size = 0;
}

static void turboqoa_work_buffer_destroy(struct TurboQOAWorkBuffer *work_buffer)
{
    free(work_buffer->data);
}

static void turboqoa_work_buffer_resize(struct TurboQOAWorkBuffer *work_buffer, size_t new_size)
{
    if(work_buffer->capacity < new_size) {
        work_buffer->data = realloc(work_buffer->data, new_size);
        work_buffer->capacity = new_size;
    }
}

static void turboqoa_work_buffer_wipe_and_copy(struct TurboQOAWorkBuffer *work_buffer, const uint8_t *data, size_t size)
{
    turboqoa_work_buffer_resize(work_buffer, size);
    memcpy(work_buffer->data, data, size);
    work_buffer->size = size;
}

static void turboqoa_work_buffer_append(struct TurboQOAWorkBuffer *work_buffer, const uint8_t *data, size_t size)
{
    size_t new_size = work_buffer->size + size;
    turboqoa_work_buffer_resize(work_buffer, new_size);
    memcpy(work_buffer->data + work_buffer->size, data, size);
    work_buffer->size = new_size;
}

static uint8_t* turboqoa_work_buffer_push_or_passthrough(struct TurboQOAWorkBuffer *work_buffer, size_t req_size, const uint8_t *data, size_t size, size_t* taken_from_data)
{
    if(work_buffer->size == 0) {
        if(size >= req_size) {
            *taken_from_data = req_size;
            return (uint8_t*)data;
        }
        else {
            *taken_from_data = size;
            turboqoa_work_buffer_wipe_and_copy(work_buffer, data, size);
            return NULL;
        }
    }
    else {
        size_t remaining = req_size - work_buffer->size;
        if(size >= remaining) {
            // Should not happen
            fprintf(stderr, "Error: buffer already contains more data than requested\n");
            abort();
        }
        else {
            *taken_from_data = remaining;
            turboqoa_work_buffer_append(work_buffer, data, size);
            return NULL;
        }
    }
}

static void turboqoa_work_buffer_clear(struct TurboQOAWorkBuffer *work_buffer)
{
    work_buffer->size = 0;
}

struct TurboQOADecoder *turboqoa_decoder_create()
{
    struct TurboQOADecoder *decoder = malloc(sizeof(struct TurboQOADecoder));
    decoder->state = TURBOQOA_DECODER_STATE_EXPECT_FILE_HEADER;
    decoder->total_decoded_samples_per_channel = 0;
    decoder->expected_channel_count = (size_t)-1;
    decoder->expected_sample_rate = (size_t)-1;
    turboqoa_work_buffer_init(&decoder->work_buffer);
    return decoder;
}

void turboqoa_decoder_destroy(struct TurboQOADecoder *decoder)
{
    turboqoa_work_buffer_destroy(&decoder->work_buffer);
    free(decoder->lms_states);
    free(decoder);
}

uint32_t loadu32be(const uint8_t *data)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return (uint32_t)data[0] << 24 | (uint32_t)data[1] << 16 | (uint32_t)data[2] << 8 | (uint32_t)data[3];
#else
    return *(uint32_t*)data;
#endif
}

uint32_t loadu24be(const uint8_t *data)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return (uint32_t)data[0] << 16 | (uint32_t)data[1] << 8 | (uint32_t)data[2];
#else
    return *(uint32_t*)data & 0xFFFFFF;
#endif
}

uint16_t loadu16be(const uint8_t *data)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return (uint16_t)data[0] << 8 | (uint16_t)data[1];
#else
    return *(uint16_t*)data;
#endif
}

uint64_t loadu64be(const uint8_t *data)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return (uint64_t)data[0] << 56 | (uint64_t)data[1] << 48 | (uint64_t)data[2] << 40 | (uint64_t)data[3] << 32 | (uint64_t)data[4] << 24 | (uint64_t)data[5] << 16 | (uint64_t)data[6] << 8 | (uint64_t)data[7];
#else
    return *(uint64_t*)data;
#endif
}

static enum TurboQOADecoderError turboqoa_decoder_decode_file_header(struct TurboQOADecoder *self, const uint8_t *data)
{
    if(memcmp(QOA_MAGIC, data, 4) != 0) {
        self->state = TURBOQOA_DECODER_STATE_EXPECT_FILE_HEADER;
        return TURBOQOA_DECODER_ERROR_INVALID_HEADER;
    }

    uint32_t sameples_per_channel = loadu32be(data + 4);
    self->total_samples_per_channel = sameples_per_channel;
    self->state = TURBOQOA_DECODER_STATE_EXPECT_FRAME_HEADER;
    self->lms_states = NULL;
    return TURBOQOA_DECODER_ERROR_NONE;
}

static const int qoa_dequant_tab[16][8] = {
	{   1,    -1,    3,    -3,    5,    -5,     7,     -7},
	{   5,    -5,   18,   -18,   32,   -32,    49,    -49},
	{  16,   -16,   53,   -53,   95,   -95,   147,   -147},
	{  34,   -34,  113,  -113,  203,  -203,   315,   -315},
	{  63,   -63,  210,  -210,  378,  -378,   588,   -588},
	{ 104,  -104,  345,  -345,  621,  -621,   966,   -966},
	{ 158,  -158,  528,  -528,  950,  -950,  1477,  -1477},
	{ 228,  -228,  760,  -760, 1368, -1368,  2128,  -2128},
	{ 316,  -316, 1053, -1053, 1895, -1895,  2947,  -2947},
	{ 422,  -422, 1405, -1405, 2529, -2529,  3934,  -3934},
	{ 548,  -548, 1828, -1828, 3290, -3290,  5117,  -5117},
	{ 696,  -696, 2320, -2320, 4176, -4176,  6496,  -6496},
	{ 868,  -868, 2893, -2893, 5207, -5207,  8099,  -8099},
	{1064, -1064, 3548, -3548, 6386, -6386,  9933,  -9933},
	{1286, -1286, 4288, -4288, 7718, -7718, 12005, -12005},
	{1536, -1536, 5120, -5120, 9216, -9216, 14336, -14336},
};

enum TurboQOADecoderError turboqoa_decoder_decode_step(struct TurboQOADecoder *self, const uint8_t *data, size_t size, size_t* input_consumed, int16_t* output, size_t max_output_size, size_t* samples_written, enum TruboQOADecoderWants* wants)
{
    *samples_written = 0;
    *input_consumed = 0;
    if(self->state == TURBOQOA_DECODER_STATE_EXPECT_FILE_HEADER) {
        uint8_t* ptr = turboqoa_work_buffer_push_or_passthrough(&self->work_buffer, QOA_FILE_HEADER_SIZE, data, size, input_consumed);
        if(ptr == NULL) {
            *wants = TURBOQOA_DECODER_WANTS_MORE_DATA;
            return TURBOQOA_DECODER_ERROR_NONE;
        }
        else {
            enum TurboQOADecoderError e = turboqoa_decoder_decode_file_header(self, ptr);
            turboqoa_work_buffer_clear(&self->work_buffer);
            if(e != TURBOQOA_DECODER_ERROR_NONE) {
                return e;
            }
            self->state = TURBOQOA_DECODER_STATE_EXPECT_FRAME_HEADER;
            *wants = TURBOQOA_DECODER_WANTS_CONTINUE_DECODING;
            return TURBOQOA_DECODER_ERROR_NONE;
        }
    }
    else if(self->state == TURBOQOA_DECODER_STATE_EXPECT_FRAME_HEADER) {    
        uint8_t* ptr = turboqoa_work_buffer_push_or_passthrough(&self->work_buffer, QOA_FRAME_HEADER_SIZE, data, size, input_consumed);
        if(ptr == NULL) {
            *wants = TURBOQOA_DECODER_WANTS_MORE_DATA;
            return TURBOQOA_DECODER_ERROR_NONE;
        }
        else {
            uint8_t num_channels = ptr[0];
            uint32_t samplerate = loadu24be(ptr + 1);
            uint16_t samples_per_channel = loadu16be(ptr + 4);
            uint16_t frame_size = loadu16be(ptr + 6);
            turboqoa_work_buffer_clear(&self->work_buffer);

            self->current_frame_header.num_channels = num_channels;
            self->current_frame_header.samplerate = samplerate;
            self->current_frame_header.samples_per_channel = samples_per_channel;
            self->current_frame_header.frame_size = frame_size;
            if(num_channels == 0 || samplerate == 0 || samples_per_channel == 0 || frame_size > QOA_FRAME_HEADER_SIZE + QOA_LMS_STATE_SIZE * num_channels + QOA_SLICE_SIZE * QOA_SLICE_PER_CHANNEL_PER_FRAME * num_channels) {
                return TURBOQOA_DECODER_ERROR_INVALID_FRAME_HEADER;
            }

            // Protect against inconsistent number of channels
            if(self->expected_channel_count != (size_t)-1) {
                if(self->expected_channel_count != self->current_frame_header.num_channels) {
                    self->state = TURBOQOA_DECODER_STATE_EXPECT_FRAME_HEADER;
                    return TURBOQOA_DECODER_ERROR_INCONSISTENT_NUMBER_OF_CHANNELS;
                }
            }
            else {
                self->expected_channel_count = self->current_frame_header.num_channels;
            }

            // same for sample rate
            if(self->expected_sample_rate != (size_t)-1) {
                if(self->expected_sample_rate != self->current_frame_header.samplerate) {
                    self->state = TURBOQOA_DECODER_STATE_EXPECT_FRAME_HEADER;
                    return TURBOQOA_DECODER_ERROR_INCONSISTENT_SAMPLE_RATE;
                }
            }
            else {
                self->expected_sample_rate = self->current_frame_header.samplerate;
            }

            self->state = TURBOQOA_DECODER_STATE_EXPECT_LMS_STATE;
            self->bytes_in_frame_read = QOA_FRAME_HEADER_SIZE;
            *wants = TURBOQOA_DECODER_WANTS_CONTINUE_DECODING;
            return TURBOQOA_DECODER_ERROR_NONE;
        }
    }
    else if(self->state == TURBOQOA_DECODER_STATE_EXPECT_LMS_STATE) {
        uint8_t* ptr = turboqoa_work_buffer_push_or_passthrough(&self->work_buffer, QOA_LMS_STATE_SIZE * self->current_frame_header.num_channels, data, size, input_consumed);
        if(ptr == NULL) {
            *wants = TURBOQOA_DECODER_WANTS_MORE_DATA;
            return TURBOQOA_DECODER_ERROR_NONE;
        }
        else {
            // TODO: Optimize this free + malloc away, it's not necessary in most cases
            free(self->lms_states);
            self->lms_states = malloc(sizeof(struct TurboQOALMSState) * self->current_frame_header.num_channels);

            for(int i = 0; i < self->current_frame_header.num_channels; i++) {
                struct TurboQOALMSState *state = &self->lms_states[i];
                uint8_t* lms_ptr = ptr + i * 16;
                for(int j = 0; j < 4; j++) {
                    state->history[j] = loadu16be(lms_ptr + j * 2);
                    state->weights[j] = loadu16be(lms_ptr + 8 + j * 2);
                }
            }

            turboqoa_work_buffer_clear(&self->work_buffer);
            self->state = TURBOQOA_DECODER_STATE_EXPECT_SLICES;
            *wants = TURBOQOA_DECODER_WANTS_CONTINUE_DECODING;
            self->slices_decoded = 0;
            self->bytes_in_frame_read += QOA_LMS_STATE_SIZE * self->current_frame_header.num_channels;
            return TURBOQOA_DECODER_ERROR_NONE;
        }
    }
    else if(self->state == TURBOQOA_DECODER_STATE_EXPECT_SLICES) {
        // Each slice contains 20 samples per channel
        size_t output_max_slices = max_output_size / self->current_frame_header.num_channels / QOA_SAMPLES_PER_SLICE;
        if(output_max_slices == 0) {
            *wants = TURBOQOA_DECODER_WANTS_MORE_OUTPUT_BUFFER;
            return TURBOQOA_DECODER_ERROR_NONE;
        }

        size_t buf_contain_slices = (self->work_buffer.size + size) / QOA_SLICE_SIZE;
        if(buf_contain_slices == 0) {
            turboqoa_work_buffer_append(&self->work_buffer, data, size);
            *input_consumed = size;
            *wants = TURBOQOA_DECODER_WANTS_MORE_DATA;
            return TURBOQOA_DECODER_ERROR_NONE;
        }

        size_t slices_remain_in_frame = (self->current_frame_header.frame_size - self->bytes_in_frame_read) / QOA_SLICE_SIZE;
        size_t slices_can_be_decoded = MIN(MIN(output_max_slices, slices_remain_in_frame), buf_contain_slices);

        if(slices_can_be_decoded < self->current_frame_header.num_channels) {
            turboqoa_work_buffer_append(&self->work_buffer, data, size);
            *input_consumed = size;
            *wants = TURBOQOA_DECODER_WANTS_MORE_DATA;
            return TURBOQOA_DECODER_ERROR_NONE;
        }

        size_t slices_to_decode = slices_can_be_decoded - slices_can_be_decoded % self->current_frame_header.num_channels;
        size_t bytes_to_decode = slices_to_decode * 8;
        uint8_t* ptr = turboqoa_work_buffer_push_or_passthrough(&self->work_buffer, bytes_to_decode, data, size, input_consumed);
        assert(ptr != NULL); // this case should be handled by the previous if

        size_t slice_per_channel = slices_to_decode / self->current_frame_header.num_channels;
        for(int i = 0; i < slice_per_channel; i++) {
            for(int j = 0; j < self->current_frame_header.num_channels; j++) {
                struct TurboQOALMSState *state = &self->lms_states[j];
                // top 4 bits of the first byte are the quantized scaling factor
                uint8_t sf_quant = ptr[0] >> 4;
                float sf = roundf(powf(sf_quant + 1, 2.75f));

                // afterwards every 3 bits are the index of the dequantization table
                uint64_t slice = (uint64_t)loadu64be(ptr);
                ptr += 8;

                for(int k = 0; k < QOA_SAMPLES_PER_SLICE /*20*/; k++) {
                    uint8_t quant = (slice >> (57 - k * 3)) & 0b111;
                    int r = qoa_dequant_tab[sf_quant][quant];
                    int32_t p = 0;
                    for(int l = 0; l < 4; l++) {
                        p += state->weights[l] * state->history[l];
                    }
                    p >>= 13;
                    short s = MIN(MAX(r + p, SHRT_MIN), SHRT_MAX);

                    // LMS update
                    int32_t delta = r >> 4;
                    for(int l = 0; l < 4; l++) {
                        state->weights[l] += state->history[l] < 0 ? -delta : delta;
                    }
                    for(int l = 0; l < 3; l++) {
                        state->history[l] = state->history[l + 1];
                    }
                    state->history[3] = s;

                    output[(i * 20 + k) * self->current_frame_header.num_channels + j] = s;
                }
                self->slices_decoded++;
            }
        }
        turboqoa_work_buffer_clear(&self->work_buffer);
        self->bytes_in_frame_read += bytes_to_decode;
        size_t decoded_samples = slices_to_decode * QOA_SAMPLES_PER_SLICE;
        if(self->total_decoded_samples_per_channel + decoded_samples / self->current_frame_header.num_channels > self->total_samples_per_channel) {
            decoded_samples = (self->total_samples_per_channel - self->total_decoded_samples_per_channel) * self->current_frame_header.num_channels;
        }
        *samples_written = decoded_samples;
        self->total_decoded_samples_per_channel += decoded_samples / self->current_frame_header.num_channels;

        if(self->slices_decoded == QOA_SLICE_PER_CHANNEL_PER_FRAME * self->current_frame_header.num_channels) {
            self->state = TURBOQOA_DECODER_STATE_EXPECT_FRAME_HEADER;
            *wants = TURBOQOA_DECODER_WANTS_CONTINUE_DECODING;
        }
        else {
            *wants = TURBOQOA_DECODER_WANTS_MORE_OUTPUT_BUFFER;
        }
        return TURBOQOA_DECODER_ERROR_NONE;
    }
    return TURBOQOA_DECODER_ERROR_INTERNAL_ERROR;
}

enum TurboQOADecoderError turboqoa_decoder_decode(struct TurboQOADecoder *self, const uint8_t *data, size_t size, size_t* input_consumed, int16_t* output, size_t max_output_size, size_t* samples_written, enum TruboQOADecoderWants* wants)
{
    *input_consumed = 0;
    *samples_written = 0;
    size_t accum_written = 0;
    size_t accum_consumed = 0;
    while(!turboqoa_decoder_decode_done(self)) {
        *samples_written = 0;
        enum TurboQOADecoderError e = turboqoa_decoder_decode_step(self, 
            data + accum_consumed,
            size - accum_consumed,
            input_consumed,
            output + accum_written,
            max_output_size - accum_written,
            samples_written,
            wants);
        accum_consumed += *input_consumed;
        accum_written += *samples_written;
        if(e != TURBOQOA_DECODER_ERROR_NONE) {
            *input_consumed = accum_consumed;
            *samples_written = accum_written;
            return e;
        }
        if(*wants != TURBOQOA_DECODER_WANTS_CONTINUE_DECODING) {
            *input_consumed = accum_consumed;
            *samples_written = accum_written;
            return TURBOQOA_DECODER_ERROR_NONE;
        }
    }

    *input_consumed = accum_consumed;
    *samples_written = accum_written;
    return TURBOQOA_DECODER_ERROR_NONE;
}

int16_t* turboqoa_decode_buffer(const uint8_t* data, size_t size, uint8_t* num_channels, uint32_t* sample_rate)
{
    struct TurboQOADecoder *decoder = turboqoa_decoder_create();
    size_t samples_written;
    enum TruboQOADecoderWants wants;


    size_t total_consumed = 0;
    size_t total_output_written = 0;
    size_t output_buffer_size = size * 16 / 3; // QOA is ~3.2bps vs 16bps for PCM
    int16_t* output = malloc(output_buffer_size * sizeof(int16_t));
    while(!turboqoa_decoder_decode_done(decoder)) {
        size_t consumed = 0;
        size_t output_written = 0;
        enum TurboQOADecoderError error = turboqoa_decoder_decode(
            decoder,
            data + total_consumed,
            size - total_consumed,
            &consumed,
            output + total_output_written,
            output_buffer_size - total_output_written,
            &samples_written,
            &wants);
        if(error != TURBOQOA_DECODER_ERROR_NONE) {
            free(output);
            turboqoa_decoder_destroy(decoder);
            return NULL;
        }
        if(wants == TURBOQOA_DECODER_WANTS_MORE_DATA) {
            free(output);
            turboqoa_decoder_destroy(decoder);
            return NULL;
        }

        // Resize output buffer if needed
        if(decoder->total_samples_per_channel != 0 && decoder->expected_channel_count != (size_t)-1
            && output_buffer_size != decoder->total_samples_per_channel * decoder->expected_channel_count) {
            size_t n = decoder->total_samples_per_channel * decoder->expected_channel_count;
            output_buffer_size = n;
            output = realloc(output, output_buffer_size * sizeof(int16_t));
        }
        if(wants == TURBOQOA_DECODER_WANTS_MORE_OUTPUT_BUFFER) {
            output_buffer_size += output_buffer_size / 2;
            output = realloc(output, output_buffer_size * sizeof(int16_t));
        }

        total_consumed += consumed;
        total_output_written += samples_written;
        assert(consumed <= size);
    }

    turboqoa_decoder_destroy(decoder);
    return NULL;
}

int turboqoa_decoder_decode_done(struct TurboQOADecoder *decoder)
{
    return decoder->total_decoded_samples_per_channel >= decoder->total_samples_per_channel
        || decoder->total_samples_per_channel == 0; // streaming mode
}

static void storeu16be(uint8_t *data, uint16_t value)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    data[0] = value >> 8;
    data[1] = value & 0xFF;
#else
    *(uint16_t*)data = value;
#endif
}

static void storeu24be(uint8_t *data, uint32_t value)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    data[0] = value >> 16;
    data[1] = (value >> 8) & 0xFF;
    data[2] = value & 0xFF;
#else
    *(uint32_t*)data = value;
#endif
}

static void storeu32be(uint8_t *data, uint32_t value)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    data[0] = value >> 24;
    data[1] = (value >> 16) & 0xFF;
    data[2] = (value >> 8) & 0xFF;
    data[3] = value & 0xFF;
#else
    *(uint32_t*)data = value;
#endif
}

static void storeu64be(uint8_t *data, uint64_t value)
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    data[0] = value >> 56;
    data[1] = (value >> 48) & 0xFF;
    data[2] = (value >> 40) & 0xFF;
    data[3] = (value >> 32) & 0xFF;
    data[4] = (value >> 24) & 0xFF;
    data[5] = (value >> 16) & 0xFF;
    data[6] = (value >> 8) & 0xFF;
    data[7] = value & 0xFF;
#else
    *(uint64_t*)data = value;
#endif
}

static void turboqoa_encoder_write(struct TurboQOAEncoder *self, const void *data, size_t size)
{
    self->write(self->io_user_data, data, size);
}

struct TurboQOAEncoder *turboqoa_encoder_create(int32_t sample_rate, uint8_t num_channels, uint32_t total_samples_per_channel, void (*write_cb)(void* user_data, const uint8_t* data, size_t size), void* user_data)
{
    struct TurboQOAEncoder *encoder = malloc(sizeof(struct TurboQOAEncoder));
    encoder->sample_rate = sample_rate;
    encoder->num_channels = num_channels;
    encoder->total_samples_per_channel = total_samples_per_channel;
    encoder->write = write_cb;
    encoder->io_user_data = user_data;
    encoder->total_samples_written_per_channel = 0;
    turboqoa_work_buffer_init(&encoder->work_buffer);

    encoder->lms_states = malloc(sizeof(struct TurboQOALMSState) * num_channels);
    for(int i = 0; i < num_channels; i++) {
        struct TurboQOALMSState *lms = &encoder->lms_states[i];
        lms->weights[0] = 0;
        lms->weights[1] = 0;
        lms->weights[2] = -(1<<13);
        lms->weights[3] =  (1<<14);
        for(int j = 0; j < 4; j++) {
            lms->history[j] = 0;
        }
    }

    // Write file header
    uint8_t file_header[8];
    memcpy(file_header, QOA_MAGIC, 4);
    storeu32be(file_header + 4, total_samples_per_channel);
    turboqoa_encoder_write(encoder, file_header, 8);

    return encoder;
}

void turboqoa_encoder_destroy(struct TurboQOAEncoder *encoder)
{
    turboqoa_work_buffer_destroy(&encoder->work_buffer);
    free(encoder->lms_states);
    free(encoder);
}

static const int qoa_reciprocal_tab[16] = {
	65536, 9363, 3121, 1457, 781, 475, 311, 216, 156, 117, 90, 71, 57, 47, 39, 32
};

static inline int qoa_div(int v, int scalefactor) {
	int reciprocal = qoa_reciprocal_tab[scalefactor];
	int n = (v * reciprocal + (1 << 15)) >> 16;
	n = n + ((v > 0) - (v < 0)) - ((n > 0) - (n < 0)); /* round away from 0 */
	return n;
}

static const int quant_tab[17] = {
	7, 7, 7, 5, 5, 3, 3, 1, /* -8..-1 */
	0,                      /*  0     */
	0, 2, 2, 4, 4, 6, 6, 6  /*  1.. 8 */
};


enum TurboQOAEncoderError turboqoa_encoder_encode_step(struct TurboQOAEncoder *self, const int16_t* data, size_t size, size_t* input_consumed, enum TurboQOAEncoderWants* wants)
{
    const size_t samples_per_frame = QOA_SLICE_PER_CHANNEL_PER_FRAME * QOA_SAMPLES_PER_SLICE * self->num_channels;
    const size_t remaining_samples_per_channel = self->total_samples_per_channel == 0 ? (size_t)-1 : self->total_samples_per_channel - self->total_samples_written_per_channel;
    const size_t encoding_samples = MIN(samples_per_frame, remaining_samples_per_channel * self->num_channels);
    uint8_t* ptr = turboqoa_work_buffer_push_or_passthrough(&self->work_buffer, encoding_samples * 2, (const uint8_t*)data, size * 2, input_consumed);
    *input_consumed /= 2;
    if(ptr == NULL) {
        *wants = TURBOQOA_ENCODER_WANTS_MORE_DATA;
        return TURBOQOA_ENCODER_ERROR_NONE;
    }

    // now ptr points to the first sample of the frame. Write frame header
    const size_t encode_sample_per_channel = encoding_samples / self->num_channels;
    const size_t encode_slices = encode_sample_per_channel / QOA_SAMPLES_PER_SLICE + (encode_sample_per_channel % QOA_SAMPLES_PER_SLICE != 0);
    uint8_t frame_header[8];
    frame_header[0] = self->num_channels;
    storeu24be(frame_header + 1, self->sample_rate);
    storeu16be(frame_header + 4, encoding_samples / self->num_channels);
    storeu16be(frame_header + 6, 8 + 16 * self->num_channels + 8 * encode_slices * self->num_channels); // frame size (including header)
    turboqoa_encoder_write(self, frame_header, 8);

    // store LMS states
    for(int i = 0; i < self->num_channels; i++) {
        struct TurboQOALMSState *lms = &self->lms_states[i];
        uint8_t buf[16];
        for(int j = 0; j < 4; j++) {
            storeu16be(buf + j * 2, lms->history[j]);
            storeu16be(buf + 8 + j * 2, lms->weights[j]);
        }
        turboqoa_encoder_write(self, buf, 16);
    }
    int prevsf[self->num_channels];
    memset(prevsf, 0, sizeof(int) * self->num_channels);

    // Compute LMS states
    int16_t* input_samples = (int16_t*)ptr;
    for(size_t sample_idx_in_channel = 0; sample_idx_in_channel < encoding_samples / self->num_channels; sample_idx_in_channel+=QOA_SAMPLES_PER_SLICE) {
        for(size_t channel = 0; channel < self->num_channels; channel++) {
            // search through all scale factors to find the best one
            int32_t best_sf = 0;
            uint64_t best_slice = 0;
            uint64_t best_error = (uint64_t)-1;
            size_t best_rank = (size_t)-1;

            struct TurboQOALMSState best_lms;
            const size_t slice_samples = MIN(QOA_SAMPLES_PER_SLICE, encoding_samples - sample_idx_in_channel * self->num_channels - channel);
            for(int sfi = 0; sfi < 16; sfi++) {
                size_t sample_begin = sample_idx_in_channel * self->num_channels + channel;
                // scale factor is heavily correlated with the previous scale factor so we start searching near the previous one
                // -1 because it might also be useful to try first
                int32_t scalefactor = (sfi + prevsf[channel]) % 16 - 1;
                if(scalefactor < 0) {
                    scalefactor += 16;
                }
                uint64_t slice = scalefactor;
                struct TurboQOALMSState lms = self->lms_states[channel];
                int64_t current_rank = 0;
                uint64_t current_error = 0;
                for(size_t sample_in_slice = 0; sample_in_slice < slice_samples; sample_in_slice++) {
                    size_t idx = sample_begin + sample_in_slice * self->num_channels;
                    int32_t sample = input_samples[idx];
                    int32_t prediction = 0;
                    for(int i = 0; i < 4; i++) {
                        prediction += lms.weights[i] * lms.history[i];
                    }
                    prediction >>= 13;
                    int32_t residual = sample - prediction;

                    int scaled = qoa_div(residual, scalefactor);
                    int clamped = MIN(MAX(scaled, -8), 8);
                    int quantized = quant_tab[clamped + 8];
                    int dequantized = qoa_dequant_tab[scalefactor][quantized];
                    int reconstructed = MIN(MAX(dequantized + prediction, SHRT_MIN), SHRT_MAX);

                    int weights_penalty = ((
						lms.weights[0] * lms.weights[0] + 
						lms.weights[1] * lms.weights[1] + 
						lms.weights[2] * lms.weights[2] + 
						lms.weights[3] * lms.weights[3]
					) >> 18) - 0x8ff;

                    weights_penalty = MAX(weights_penalty, 0);

                    int64_t error = sample - reconstructed;
                    uint64_t l2error = error * error;

                    current_rank += l2error + weights_penalty * weights_penalty;
                    current_error += l2error;

                    if (current_rank > best_rank) {
						break;
					}

                    int delta = dequantized >> 4;
                    for(int i = 0; i < 4; i++) {
                        lms.weights[i] += lms.history[i] < 0 ? -delta : delta;
                    }
                    for(int i = 0; i < 3; i++) {
                        lms.history[i] = lms.history[i + 1];
                    }
                    lms.history[3] = reconstructed;
                    slice = (slice << 3) | quantized;

                }
                if(current_rank < best_rank) {
                    best_rank = current_rank;
                    best_sf = scalefactor;
                    best_slice = slice;
                    best_lms = lms;
                    best_error = current_error;
                }
            }
            self->lms_states[channel] = best_lms;
            prevsf[channel] = best_sf;
            uint8_t slice_buf[8];
            if(slice_samples != QOA_SAMPLES_PER_SLICE) {
                best_slice <<= (QOA_SAMPLES_PER_SLICE - slice_samples) * 3;
            }
            storeu64be(slice_buf, best_slice);
            turboqoa_encoder_write(self, slice_buf, 8);
        }
    }
    turboqoa_work_buffer_clear(&self->work_buffer);

    self->total_samples_written_per_channel += encode_sample_per_channel;
    *wants = TURBOQOA_ENCODER_WANTS_CONTINUE_ENCODING;
    return TURBOQOA_ENCODER_ERROR_NONE;
}

enum TurboQOAEncoderError turboqoa_encoder_encode(struct TurboQOAEncoder *self, const int16_t* data, size_t size, size_t* input_consumed, enum TurboQOAEncoderWants* wants)
{
    *input_consumed = 0;
    *wants = TURBOQOA_ENCODER_WANTS_CONTINUE_ENCODING;
    size_t accum_consumed = 0;
    while(!turboqoa_encoder_encode_done(self)) {
        size_t consumed = 0;
        enum TurboQOAEncoderError e = turboqoa_encoder_encode_step(self, data + accum_consumed, size - accum_consumed, &consumed, wants);
        accum_consumed += consumed;
        if(e != TURBOQOA_ENCODER_ERROR_NONE) {
            *input_consumed = accum_consumed;
            return e;
        }
        if(*wants != TURBOQOA_ENCODER_WANTS_CONTINUE_ENCODING) {
            *input_consumed = accum_consumed;
            return TURBOQOA_ENCODER_ERROR_NONE;
        }
    }
    *input_consumed = accum_consumed;
    return TURBOQOA_ENCODER_ERROR_NONE;
}

int turboqoa_encoder_encode_done(struct TurboQOAEncoder *encoder)
{
    return encoder->total_samples_written_per_channel >= encoder->total_samples_per_channel || encoder->total_samples_per_channel == 0;
}

struct EncoderBufferWriter {
    uint8_t* buffer;
    size_t size;
    size_t offset;
};

void write_to_buffer(void* user_data, const uint8_t* data, size_t size)
{
    struct EncoderBufferWriter *writer = (struct EncoderBufferWriter*)user_data;
    if(writer->offset + size > writer->size) {
        writer->size = writer->offset + size;
        writer->buffer = realloc(writer->buffer, writer->size);
    }

    memcpy(writer->buffer + writer->offset, data, size);
}

uint8_t* turboqoa_encode_buffer(const int16_t* data, size_t size, uint8_t num_channels, uint32_t sample_rate, size_t* out_size)
{
    assert(size % num_channels == 0);
    struct EncoderBufferWriter writer;
    writer.buffer = malloc(1024);
    writer.size = 1024;
    writer.offset = 0;
    struct TurboQOAEncoder *encoder = turboqoa_encoder_create(sample_rate, num_channels, size / num_channels, write_to_buffer, &writer);

    size_t total_consumed = 0;
    while(!turboqoa_encoder_encode_done(encoder)) {
        size_t consumed = 0;
        enum TurboQOAEncoderError error = turboqoa_encoder_encode(encoder, data + total_consumed, size - total_consumed, &consumed, NULL);
        if(error != TURBOQOA_ENCODER_ERROR_NONE) {
            free(writer.buffer);
            turboqoa_encoder_destroy(encoder);
            return NULL;
        }
        total_consumed += consumed;
    }
    return writer.buffer;
}