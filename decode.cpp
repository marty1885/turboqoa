#include "turboqoa/turboqoa.h"
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
    if(argc < 3) {
        std::cout << "Usage: " << argv[0] << " <input.qoa> [<output.pcm>]" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];
    std::string output_path = argc > 2 ? argv[2] : "output.pcm";
    std::vector<uint8_t> data = read_file(input_path);

    struct TurboQOADecoder *decoder = turboqoa_decoder_create();
    int16_t output[1024];
    size_t samples_written;
    enum TruboQOADecoderWants wants;

    uint8_t* data_ptr = data.data();
    size_t consumed = 0;

    std::ofstream out(output_path, std::ios::binary);
    size_t total_samples = 0;
    size_t total_consumed = 0;
    while(!turboqoa_decoder_decode_done(decoder)) {
        enum TurboQOADecoderError error = turboqoa_decoder_decode(decoder, data.data() + total_consumed, data.size() - total_consumed, &consumed, output, 1024, &samples_written, &wants);
        if(error != TURBOQOA_DECODER_ERROR_NONE) {
            std::cout << "Error while decoding. Code: " << error << std::endl;
            break;
        }
        if(wants == TURBOQOA_DECODER_WANTS_MORE_DATA) {
            std::cout << "Decoder wants more data" << std::endl;
            break;
        }
        total_consumed += consumed;
        assert(consumed <= data.size());

        if(samples_written > 0) {
            out.write((char*)output, samples_written * sizeof(int16_t));
            total_samples += samples_written;
        }

        data_ptr += consumed;
        std::cout << "Progress: " << total_consumed << " / " << data.size() << std::endl;
    }
    std::cout << "Decoded " << total_samples << " samples" << std::endl;

    turboqoa_decoder_destroy(decoder);
}