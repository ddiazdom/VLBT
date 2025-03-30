cmake_minimum_required(VERSION 3.16.3)
project(bwt_dts_benchmarks)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Modules")

#find_package(LibSDSL REQUIRED)

include_directories(SYSTEM ${LIBSDSL_INCLUDE_DIRS}  ${CMAKE_CURRENT_SOURCE_DIR}/fb_wt/wt-fbb-0.1.0)

add_executable(test_bwt_dts test_bwt_dts.cpp
        #../cdt_common.hpp
        #../cdt_common.cpp
)
target_compile_options(test_bwt_dts
        PUBLIC
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function #-Winline
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math -mavx2
)
#target_include_directories(test_bwts PUBLIC ${LIBSDSL_INCLUDE_DIRS})
target_link_libraries(test_bwt_dts LINK_PUBLIC ${LIBSDSL_LIBRARIES})

#build BWT data structures
add_executable(build_bwt_dts build_bwt_dts.cpp
        #../cdt_common.hpp
        ../cdt_common.cpp
)
target_compile_options(build_bwt_dts
        PUBLIC
        #-Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function
        #-O3 -funroll-loops -fomit-frame-pointer -ffast-math
        #-O3 -funroll-loops -mavx2
)

#target_link_libraries(build_bwt_dts LINK_PUBLIC ${LIBSDSL_LIBRARIES})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
    target_compile_options(test_bwt_dts PUBLIC -march=native)
    target_compile_options(build_bwt_dts PUBLIC -march=native)
endif()


