#pragma once

#include <stddef.h>
#include <stdint.h>
#include <unistd.h>

#ifdef __cplusplus
extern "C" {
#endif

enum TurboQOADecoderState
{
    TURBOQOA_DECODER_STATE_EXPECT_FILE_HEADER,
    TURBOQOA_DECODER_STATE_EXPECT_FRAME_HEADER,
    TURBOQOA_DECODER_STATE_EXPECT_LMS_STATE,
    TURBOQOA_DECODER_STATE_EXPECT_SLICES
};

enum TruboQOADecoderWants
{
    TURBOQOA_DECODER_WANTS_CONTINUE_DECODING,
    TURBOQOA_DECODER_WANTS_MORE_DATA,
    TURBOQOA_DECODER_WANTS_MORE_OUTPUT_BUFFER,
};

enum TurboQOADecoderError
{
    TURBOQOA_DECODER_ERROR_NONE,
    TURBOQOA_DECODER_ERROR_INVALID_HEADER,
    TURBOQOA_DECODER_ERROR_INVALID_FRAME_HEADER,
    TURBOQOA_DECODER_ERROR_INTERNAL_ERROR,
    TURBOQOA_DECODER_ERROR_INCONSISTENT_NUMBER_OF_CHANNELS, // QOA supports different number of channels per frame. I don't want to expose user to horror of dealing with this
    TURBOQOA_DECODER_ERROR_INCONSISTENT_SAMPLE_RATE,        // DITTO as above
};

struct TurboQOAWorkBuffer
{
    uint8_t* data;
    size_t capacity;
    size_t size;
};

struct TurboQOAFrameHeader
{
    uint8_t num_channels;
    uint32_t samplerate;  // stored as 24-bit
    uint16_t samples_per_channel;
    uint16_t frame_size;
};

struct TurboQOALMSState
{
    int16_t history[4];
    int16_t weights[4];
};

struct TurboQOADecoder
{
    enum TurboQOADecoderState state;

    struct TurboQOAWorkBuffer work_buffer;
    uint32_t total_samples_per_channel;

    struct TurboQOAFrameHeader current_frame_header;
    struct TurboQOALMSState* lms_states;
    size_t slices_decoded;
    size_t bytes_in_frame_read;
    size_t total_decoded_samples_per_channel;

    size_t expected_channel_count;
    size_t expected_sample_rate;
};

struct TurboQOADecoder *turboqoa_decoder_create();
void turboqoa_decoder_destroy(struct TurboQOADecoder *decoder);
enum TurboQOADecoderError turboqoa_decoder_decode_step(struct TurboQOADecoder *decoder,
    const uint8_t *data,
    size_t size,
    size_t* consumed_input,
    int16_t* output,
    size_t max_output_size,
    size_t* samples_written,
    enum TruboQOADecoderWants* wants);
enum TurboQOADecoderError turboqoa_decoder_decode(struct TurboQOADecoder *decoder,
    const uint8_t *data,
    size_t size,
    size_t* consumed_input,
    int16_t* output,
    size_t max_output_size,
    size_t* samples_written,
    enum TruboQOADecoderWants* wants);
int turboqoa_decoder_decode_done(struct TurboQOADecoder *decoder);

struct TurboQOAEncoder
{
    void (*write)(void* user_data, const uint8_t* data, size_t size);
    void* io_user_data;

    struct TurboQOALMSState* lms_states;

    int sample_rate;
    int num_channels;
    uint32_t total_samples_per_channel;
    uint32_t total_samples_written_per_channel;

    struct TurboQOAWorkBuffer work_buffer;
    size_t n_slices_in_frame;
};

enum TurboQOAEncoderError
{
    TURBOQOA_ENCODER_ERROR_NONE,
    TURBOQOA_ENCODER_ERROR_TOO_MANY_SAMPLES,
};


enum TurboQOAEncoderWants
{
    TURBOQOA_ENCODER_WANTS_MORE_DATA,
    TURBOQOA_ENCODER_WANTS_CONTINUE_ENCODING,
};

struct TurboQOAEncoder *turboqoa_encoder_create(int32_t sample_rate, uint8_t num_channels, uint32_t total_samples_per_channel, void (*write_cb)(void* user_data, const uint8_t* data, size_t size), void* user_data);
void turboqoa_encoder_destroy(struct TurboQOAEncoder *encoder);
enum TurboQOAEncoderError turboqoa_encoder_encode(struct TurboQOAEncoder *encoder, const int16_t* data, size_t size, size_t* consumed_input, enum TurboQOAEncoderWants* wants);
int turboqoa_encoder_encode_done(struct TurboQOAEncoder *encoder);

#ifdef __cplusplus
}
#endif