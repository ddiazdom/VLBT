cmake_minimum_required(VERSION 3.16)

project(VLBT-test)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
find_package(LibSDSL REQUIRED)

include(CheckSIMDInstructions)
set(SIMD_FLAGS)
check_neon_flags(SIMD_FLAGS)
check_sse42_flags(SIMD_FLAGS)

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
target_compile_options(vlbt_test_functions PUBLIC -Wall -pedantic -O3 -g ${SIMD_FLAGS} -Wno-deprecated-declarations)
target_link_libraries(vlbt_test_functions LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(vlbt_test_functions SYSTEM PUBLIC ${LIBSDSL_INCLUDE_DIRS} r-index/internal)

#benchmark cache misses in the VLBT
add_executable(vlbt_bench_cmisses vlbt_bench_cmisses.cpp)
target_compile_options(vlbt_bench_cmisses PUBLIC -Wall -pedantic -O3 -g ${SIMD_FLAGS})

#benchmark cache misses in the RLBWTs of the SDSL
add_executable(rlbwt_bench_cmisses rlbwt_bench_cmisses.cpp)
target_compile_options(rlbwt_bench_cmisses PUBLIC -Wall -pedantic -O3 -g -Wno-deprecated-declarations)
target_link_libraries(rlbwt_bench_cmisses LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(rlbwt_bench_cmisses SYSTEM PUBLIC ${LIBSDSL_INCLUDE_DIRS})

#benchmark cache misses in the r-index
add_executable(ri_bench_cmisses ri_bench_cmisses.cpp)
target_compile_options(ri_bench_cmisses PUBLIC -Wall -pedantic -O3 -g -Wno-deprecated-declarations)
target_link_libraries(ri_bench_cmisses LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(ri_bench_cmisses SYSTEM PUBLIC ${LIBSDSL_INCLUDE_DIRS} r-index/internal)

#benchmark cache misses in the r-index
add_executable(sri_bench_cmisses sri_bench_cmisses.cpp)
target_compile_options(sri_bench_cmisses PUBLIC -Wall -pedantic -O3 -g -Wno-deprecated-declarations -Wno-unused-but-set-variable -Wno-unused-lambda-capture -Wno-unused-variable)
target_link_libraries(sri_bench_cmisses LINK_PUBLIC ${LIBSDSL_LIBRARIES})
target_include_directories(sri_bench_cmisses SYSTEM PUBLIC ${LIBSDSL_INCLUDE_DIRS} sr-index/include sr-index/include/json/include)