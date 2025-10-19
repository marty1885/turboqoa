cmake_minimum_required(VERSION 3.15)
project(turboqoa)

add_library(turboqoa STATIC turboqoa/turboqoa.c)
target_include_directories(turboqoa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

option(BUILD_EXAMPLES "Build the example executables" ON)
option(BUILD_FUZZ_TEST "Build fuzzing test" OFF)

if(BUILD_EXAMPLES)
    add_executable(decode decode.cpp)
    target_link_libraries(decode turboqoa)

    add_executable(encode encode.cpp)
    target_link_libraries(encode turboqoa)
endif()

if(BUILD_FUZZ_TEST)
    if(NOT CMAKE_C_COMPILER_ID MATCHES "Clang")
        message(FATAL_ERROR "Building the fuzz test requires the Clang compiler.")
    endif()
    add_subdirectory(fuzz)
endif()
