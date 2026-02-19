cmake_minimum_required(VERSION 3.16)

project(VLBT)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MESSAGE_LOG_LEVEL DEBUG)

# modules to detect SIMD compilations flags
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CheckSIMDInstructions)
#

# check the available CPU instructions
set(CMP_FLAGS)#compilation flags
check_neon_flags(CMP_FLAGS)
check_sse42_flags(CMP_FLAGS)
check_bmi2_flags(CMP_FLAGS)
#

if(CMP_FLAGS)
    message(STATUS "Extra compilation flags ${CMP_FLAGS}")
endif()

#the command line interface of the library
add_executable(vlbt-cli vlbt-cli.cpp)
target_compile_options(vlbt-cli
        PUBLIC
        -Wall -pedantic -O3 -g
        ${CMP_FLAGS}
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(vlbt-cli PUBLIC -march=native)

    # for the file_system functionality on Linux
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL "9.1")
        message(STATUS "GCC version is <=9.1, linking stdc++fs library")
        target_link_libraries(vlbt-cli PUBLIC stdc++fs)
    endif()
endif()
