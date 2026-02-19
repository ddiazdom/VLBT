include(CheckCXXSourceCompiles)
include(CheckCXXCompilerFlag)

function(check_sse42_flags SIMD_FLAGS)
    set(CMAKE_CXX_FLAGS_BCK "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wall -msse4.2")
    set(SSE42_CXX_CODE "
        #include <x86intrin.h>
        int main(){
        #if __SSE4_2__
        return 0;
        #else
        #error \"SSE4.2 is not supported\"
        #endif
        }
        "
    )
    check_cxx_source_compiles("${SSE42_CXX_CODE}" SSE42_SUPPORT)
    if(SSE42_SUPPORT)
        list(APPEND SIMD_FLAGS -msse4.2)
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BCK}")
endfunction()

function(check_avx2_flags SIMD_FLAGS)
    set(CMAKE_CXX_FLAGS_BCK "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wall -mavx2")
    set(AVX2_CXX_CODE "
        #include <immintrin.h>
        int main(){
        #if __AVX2__
        return 0;
        #else
        #error \"AVX2 is not supported\"
        #endif
        }"
    )
    check_cxx_source_compiles("${AVX2_CXX_CODE}" AVX2_SUPPORT)
    if(AVX2_SUPPORT)
        list(APPEND SIMD_FLAGS -mavx2)
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BCK}")
endfunction()

function(check_bmi2_flags SIMD_FLAGS)
    check_cxx_compiler_flag("-mbmi2" BMI2_SUPPORT)
    if(BMI2_SUPPORT)
        list(APPEND SIMD_FLAGS -mbmi2)
        message("Adding BMI2 as compilation flag")
    endif()
endfunction()

function(check_neon_flags SIMD_FLAGS)
    set(CMAKE_CXX_FLAGS_BCK "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "")
    set(NEON_CXX_CODE "
        #include <arm_neon.h>
        int main(){
        #if __ARM_NEON__
        return 0;
        #else
        #error \"NEON is not supported\"
        #endif
        }"
    )
    check_cxx_source_compiles("${NEON_CXX_CODE}" NEON_SUPPORT)
    ##TODO do something here
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BCK}")
endfunction()
