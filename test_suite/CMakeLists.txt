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
check_cxx_source_compiles("${SSE42_CXX_CODE}" SSE42_SUPPORT)
if(SSE42_SUPPORT)
    list(APPEND SIMD_FLAGS -msse4.2)
endif()
message("Compiling with ${SIMD_FLAGS}")

check_cxx_source_compiles("${NEON_CXX_CODE}" NEON_SUPPORT)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BCK}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
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
add_executable(vlbt_test_functions vlbt_test_functions.cpp)
target_compile_options(vlbt_test_functions PUBLIC -Wall -pedantic -O3 -ffast-math -g ${SIMD_FLAGS})
target_link_libraries(vlbt_test_functions LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(vlbt_test_functions PUBLIC ${LIBSDSL_INCLUDE_DIRS} r-index/internal)

#benchmark cache misses in the VLBT
add_executable(vlbt_bench_cmisses vlbt_bench_cmisses.cpp)
target_compile_options(vlbt_bench_cmisses PUBLIC -Wall -pedantic -O3 -ffast-math -g ${SIMD_FLAGS})

#benchmark cache misses in the RLBWTs of the SDSL
add_executable(rlbwt_bench_cmisses rlbwt_bench_cmisses.cpp)
target_compile_options(rlbwt_bench_cmisses PUBLIC -Wall -pedantic -O3 -ffast-math -g)
target_link_libraries(rlbwt_bench_cmisses LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(rlbwt_bench_cmisses PUBLIC ${LIBSDSL_INCLUDE_DIRS})

#benchmark cache misses in the r-index
add_executable(ri_bench_cmisses ri_bench_cmisses.cpp)
target_compile_options(ri_bench_cmisses PUBLIC -Wall -pedantic -O3 -ffast-math -g)
target_link_libraries(ri_bench_cmisses LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(ri_bench_cmisses PUBLIC ${LIBSDSL_INCLUDE_DIRS} r-index/internal)

#benchmark cache misses in the r-index
add_executable(sri_bench_cmisses sri_bench_cmisses.cpp)
target_compile_options(sri_bench_cmisses PUBLIC -Wall -pedantic -O3 -ffast-math -g)
target_link_libraries(sri_bench_cmisses LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(sri_bench_cmisses PUBLIC ${LIBSDSL_INCLUDE_DIRS} sr-index/include)
target_include_directories(sri_bench_cmisses SYSTEM PUBLIC sr-index/include/json/include)

