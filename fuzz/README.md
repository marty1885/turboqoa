# TurboQOA Fuzzing Test

This directory contains a fuzzing test for the TurboQOA decoder. The test uses libFuzzer and AddressSanitizer to detect memory corruption vulnerabilities.

## Building and Running the Fuzzer

The fuzzer must be built with Clang.

1.  **Install Clang and the fuzzer runtime libraries:**
    ```bash
    sudo apt-get update
    sudo apt-get install -y clang libclang-rt-18-dev
    ```

2.  **Build the fuzzer:**
    From the root of the project:
    ```bash
    mkdir build
    cd build
    CC=clang CXX=clang++ cmake .. -DBUILD_FUZZ_TEST=ON -DBUILD_EXAMPLES=OFF
    make
    ```
    *Note: The C++ examples in this project may not compile correctly with Clang in all environments. If you encounter build errors, ensure you have disabled the examples with `-DBUILD_EXAMPLES=OFF`.*

3.  **Run the fuzzer:**
    ```bash
    ./fuzz/fuzz_test -max_total_time=60
    ```
