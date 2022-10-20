cmake_minimum_required(VERSION 3.16)

project(fourc)


set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(COMPILE_FLAGS -Wall -Werror -fsanitize=address  -fno-omit-frame-pointer -pthread)
set(LINK_FLAGS -fsanitize=address -lm -llibslater)

find_package(PkgConfig)

include(ExternalProject)

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

ExternalProject_Add(libslater
        GIT_REPOSITORY https://github.com/urysegal/libslater.git
        GIT_TAG main
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
        )

include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

include_directories(${fourc_SOURCE_DIR} ${fourc_SOURCE_DIR}/include )


add_executable(fourc src/nintlib.cpp include/nintlib.h src/fourc.cpp)

add_dependencies(fourc libslater)

target_include_directories(fourc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(fourc PUBLIC ${COMPILE_FLAGS})
target_link_options(fourc PUBLIC ${LINK_FLAGS})

target_link_libraries(fourc -llibslater -lgcov)