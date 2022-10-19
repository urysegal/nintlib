cmake_minimum_required(VERSION 3.16)

project(fourc)



set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(COMPILE_FLAGS -Wall -Werror -fsanitize=address  -fno-omit-frame-pointer -pthread)
set(LINK_FLAGS -fsanitize=address -lm)

find_package(PkgConfig)

include_directories(${fourc_SOURCE_DIR} ${fourc_SOURCE_DIR}/include)


add_executable(fourc src/nintlib.c include/nintlib.h src/fourc.cpp)
target_include_directories(fourc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(fourc PUBLIC ${COMPILE_FLAGS})
target_link_options(fourc PUBLIC ${LINK_FLAGS})

