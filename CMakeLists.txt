cmake_minimum_required(VERSION 3.16.3)

project(VLBT)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_BCK "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -mavx2 -msse4.2")

include(CheckCXXSourceCompiles)
set(SSE42_CXX_CODE "
        #include <x86intrin.h>
        int main(){
        #if __SSE4_2__
        return 0;
        #else
        #error \"SSE4.2 is not supported\"
        #endif
        }
")

set(AVX2_CXX_CODE "
        #include <immintrin.h>
        int main(){
        #if __AVX2__
        return 0;
        #else
        #error \"AVX2 is not supported\"
        #endif
        }
")

set(NEON_CXX_CODE "
        #include <arm_neon.h>
        int main(){
        #if __ARM_NEON__
        return 0;
        #else
        #error \"NEON is not supported\"
        #endif
        }
")

set(SIMD_FLAGS)
check_cxx_source_compiles("${AVX2_CXX_CODE}" AVX2_SUPPORT)
#if(AVX2_SUPPORT)
#    list(APPEND SIMD_FLAGS -mavx2)
#endif()

check_cxx_source_compiles("${SSE42_CXX_CODE}" SSE42_SUPPORT)
if(SSE42_SUPPORT)
    list(APPEND SIMD_FLAGS -msse4.2)
endif()
message("Compiling with ${SIMD_FLAGS}")

check_cxx_source_compiles("${NEON_CXX_CODE}" NEON_SUPPORT)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BCK}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/test_suite/cmake/Modules")
find_package(LibSDSL REQUIRED)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    add_compile_options(-march=native)
    #for the file_system functionality on Linux
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL "9.1")
        message(STATUS "GCC version is <=9.1, adding the compiler flag stdc++fs")
        add_compile_options(stdc++fs)
    endif()
endif()

#test the VLBT suite
add_executable(vlbt_test_functions test_suite/vlbt_test_functions.cpp)
target_compile_options(vlbt_test_functions PUBLIC -Wall -pedantic -O3 -ffast-math -g ${SIMD_FLAGS} )
target_link_libraries(vlbt_test_functions LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(vlbt_test_functions PUBLIC ${LIBSDSL_INCLUDE_DIRS} scripts/r-index/internal)

#the CLI of the library
add_executable(vlbt-cli vlbt-cli.cpp)
target_compile_options(vlbt-cli
        PUBLIC
        -Wall -pedantic -O3 -ffast-math -g
        ${SIMD_FLAGS}
)