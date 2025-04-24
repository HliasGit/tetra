cmake_minimum_required(VERSION 3.10)

cmake_policy(SET CMP0048 NEW)
project(PYTHAGORAS VERSION 0.1)

set(LANGUAGES C CXX)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(VERBOSE "Enable verbose mode")

# ====== Set the build type ======
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ====== Set the compiler flags ======
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode")
    # AddressSanitizer settings
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer -static-libasan -O0 -g -Wall -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")

    # Set environment variable to catch memory leaks
    set(ENV{ASAN_OPTIONS} "detect_leaks=1")

    message(STATUS "AddressSanitizer enabled for memory leak detection")
elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
    message(STATUS "Building in profile mode")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")
else()
    message(STATUS "Building in release mode")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ffast-math")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math")
endif()



# ====== Now build the library ======
file(GLOB_RECURSE TETRA_SOURCES src/*.c)

add_library(TETRA ${TETRA_SOURCES})
target_include_directories(TETRA PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

if(VERBOSE)
    target_compile_definitions(TETRA PRIVATE VERBOSE)
endif()

add_subdirectory(test)