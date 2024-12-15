#include "turboqoa/turboqoa.h"
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

std::vector<uint8_t> read_file(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if(!in.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    in.seekg(0, std::ios::end);
    size_t size = in.tellg();
    in.seekg(0, std::ios::beg);

    std::vector<uint8_t> data(size);
    in.read((char*)data.data(), size);
    return data;
}

int main(int argc, char** argv)
{

    std::string input_path;
    std::string output_path = "encoding_test.qoa";
    if(argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input.pcm> [<output.qoa>]" << std::endl;
        return 1;
    }
    input_path = argv[1];
    if(argc > 2) {
        output_path = argv[2];
    }
    std::vector<uint8_t> data = read_file(input_path);
    size_t channels = 2;
    size_t sample_rate = 48000;
    size_t samples_per_channel = data.size() / 2 / channels;

    std::ofstream out(output_path, std::ios::binary);

    struct TurboQOAEncoder* encoder = turboqoa_encoder_create(sample_rate, channels, samples_per_channel, [](void* user_data, const uint8_t* data, size_t size) {
        std::ofstream* out = (std::ofstream*)user_data;
        out->write((char*)data, size);
    }, &out);

    size_t consumed = 0;
    enum TurboQOAEncoderWants wants;

    int16_t* data_ptr = (int16_t*)data.data();
    size_t num_samples = data.size() / 2;
    size_t total_consumed = 0;
    while(!turboqoa_encoder_encode_done(encoder)) {
        enum TurboQOAEncoderError error = turboqoa_encoder_encode(encoder, data_ptr, num_samples, &consumed, &wants);
        if(error != TURBOQOA_ENCODER_ERROR_NONE) {
            std::cout << "Error while encoding. Code: " << error << std::endl;
            break;
        }
        if(wants == TURBOQOA_ENCODER_WANTS_MORE_DATA) {
            std::cout << "Encoder wants more data" << std::endl;
            break;
        }

        data_ptr += consumed;
        num_samples -= consumed;
        total_consumed += consumed;

        printf("Processed %lu/%lu samples\n", total_consumed, samples_per_channel * channels);
    }

    turboqoa_encoder_destroy(encoder);
}