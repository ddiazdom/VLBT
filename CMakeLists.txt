cmake_minimum_required(VERSION 3.27)
project(simple_rl_bwt_2024)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Modules/")
find_package(LibSDSL REQUIRED)

add_executable(simple_rl_bwt_2024
        main.cpp
        cdt_common.hpp
        cdt_common.cpp)

target_compile_options(simple_rl_bwt_2024
        PUBLIC
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math -std=c++20
)

target_include_directories(simple_rl_bwt_2024
        PUBLIC
        ${LIBSDSL_INCLUDE_DIRS}
)


target_link_libraries(simple_rl_bwt_2024 LINK_PUBLIC ${LIBSDSL_LIBRARIES})
