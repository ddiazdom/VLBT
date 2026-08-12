cmake_minimum_required(VERSION 3.16)
project(VLBT-scripts)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MESSAGE_LOG_LEVEL DEBUG)

# modules to detect SIMD compilations flags
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
include(CheckSIMDInstructions)
#

# check the available CPU instructions
set(SIMD_FLAGS)
check_neon_flags(SIMD_FLAGS)
check_sse42_flags(SIMD_FLAGS)
#

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    #for the file_system functionality on Linux
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL "9.1")
        message(STATUS "GCC version is <=9.1, adding the compiler flag stdc++fs")
        add_compile_options(stdc++fs)
    endif()
endif()

#SOME useful scripts for developing/testing VLBT

#generate random patterns from a text
add_executable(gen_patterns gen_patterns.cpp)
target_compile_options(gen_patterns PUBLIC -O3)

#transform the patterns to fasta file
add_executable(patfile2fasta patfile2fasta.cpp)
target_compile_options(patfile2fasta PUBLIC -Wall -pedantic -O3 -g)

#filter patterns by frequency
add_executable(filter_pat_set filter_patset.cpp)
target_compile_options(filter_pat_set PUBLIC -Wall -pedantic -O3 -g ${SIMD_FLAGS})

#run-length encode a BWT
add_executable(plain2rle plain2vbyte_rle.cpp)
target_compile_options(plain2rle PUBLIC -Wall -pedantic -O3 -g)