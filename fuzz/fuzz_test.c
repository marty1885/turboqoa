#include "turboqoa/turboqoa.h"
#include <stddef.h>
#include <stdint.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 1) {
        return 0;
    }

    struct TurboQOADecoder *decoder = turboqoa_decoder_create();
    if (!decoder) {
        return 0;
    }

    int16_t output[4096];
    size_t samples_written;
    enum TruboQOADecoderWants wants;
    size_t consumed;

    turboqoa_decoder_decode(decoder, data, size, &consumed, output, 4096, &samples_written, &wants);

    turboqoa_decoder_destroy(decoder);

    return 0;
}
