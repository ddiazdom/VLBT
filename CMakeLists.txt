cmake_minimum_required(VERSION 3.16)

project(VLBT)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MESSAGE_LOG_LEVEL DEBUG)

# modules to detect SIMD compilations flags
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CheckSIMDInstructions)
#

# check the available CPU instructions
set(SIMD_FLAGS)
check_neon_flags(SIMD_FLAGS)
check_sse42_flags(SIMD_FLAGS)
#

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    add_compile_options(-march=native)
    #for the file_system functionality on Linux
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL "9.1")
        message(STATUS "GCC version is <=9.1, adding the compiler flag stdc++fs")
        add_compile_options(stdc++fs)
    endif()
endif()

#the command line interface of the library
add_executable(vlbt-cli vlbt-cli.cpp)
target_compile_options(vlbt-cli
        PUBLIC
        -Wall -pedantic -O3 -g
        ${SIMD_FLAGS}
)