cmake_minimum_required(VERSION 3.16.3)
project(bwt_dts_benchmarks)

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
if(AVX2_SUPPORT)
    list(APPEND SIMD_FLAGS -mavx2)
endif()

check_cxx_source_compiles("${SSE42_CXX_CODE}" SSE42_SUPPORT)
if(SSE42_SUPPORT)
    list(APPEND SIMD_FLAGS -msse4.2)
endif()

check_cxx_source_compiles("${NEON_CXX_CODE}" NEON_SUPPORT)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BCK}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Modules")

find_package(LibSDSL REQUIRED)

include_directories(SYSTEM ${LIBSDSL_INCLUDE_DIRS} fb_wt/wt-fbb-0.1.0)

#add_executable(test_bwt_dts test_bwt_dts.cpp ../cdt_common.cpp)
#target_compile_options(test_bwt_dts
#        PUBLIC
#        -Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function #-Winline
#        -O3 -funroll-loops -fomit-frame-pointer -ffast-math -mavx2
#)
#target_link_libraries(test_bwt_dts LINK_PUBLIC ${LIBSDSL_LIBRARIES})
#target_include_directories(test_bwts PUBLIC ${LIBSDSL_INCLUDE_DIRS})
#if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
#    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
#    target_compile_options(test_bwt_dts PUBLIC -march=native)
#endif()

#Build the BWT data structures
add_executable(build_bwt_dts build_bwt_dts.cpp)
target_compile_options(build_bwt_dts
        PUBLIC
        #-Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function
        #-O3 -funroll-loops -fomit-frame-pointer -ffast-math
        #-O3 -funroll-loops -mavx2
        #-O3
        ${SIMD_FLAGS}
)
target_link_libraries(build_bwt_dts LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(build_bwt_dts PUBLIC ${LIBSDSL_INCLUDE_DIRS})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(build_bwt_dts PUBLIC -march=native)
    #for the file_system functionality on Linux
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL "9.1")
        message(STATUS "GCC version is <=9.1, adding the compiler flag stdc++fs")
        target_compile_options(build_bwt_dts PUBLIC stdc++fs)
    endif()
endif()