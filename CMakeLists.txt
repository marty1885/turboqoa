cmake_minimum_required(VERSION 3.15)
project(turboqoa)

add_library(turboqoa STATIC turboqoa/turboqoa.c)
target_include_directories(turboqoa PUBLIC turboqoa)

add_executable(decode decode.cpp)
target_link_libraries(decode turboqoa)

add_executable(encode encode.cpp)
target_link_libraries(encode turboqoa)