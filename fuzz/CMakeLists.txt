add_executable(fuzz_test fuzz_test.c)
target_compile_options(fuzz_test PRIVATE -fsanitize=fuzzer,address -g)
target_link_libraries(fuzz_test PRIVATE turboqoa -fsanitize=fuzzer,address)
