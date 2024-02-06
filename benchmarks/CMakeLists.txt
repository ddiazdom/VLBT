cmake_minimum_required(VERSION 3.21)
project(testing_bwts)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Modules")

message(${CMAKE_MODULE_PATH})

find_package(LibSDSL REQUIRED)

add_executable(test_bwts test_bwts.cpp ../cdt_common.hpp ../cdt_common.cpp)
target_compile_options(test_bwts
        PUBLIC
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function -Winline
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math
)
target_include_directories(test_bwts PUBLIC ${LIBSDSL_INCLUDE_DIRS})
target_link_libraries(test_bwts LINK_PUBLIC ${LIBSDSL_LIBRARIES})


#build wavelet-tree data structures
add_executable(build_wt_dts build_bwt_dts.cpp ../cdt_common.hpp ../cdt_common.cpp)
target_compile_options(build_wt_dts
        PUBLIC
        -Wall -Wextra -Wno-ignored-qualifiers -pedantic -Wno-deprecated-copy -Wno-unused-function -Winline
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math
)
target_include_directories(build_wt_dts PUBLIC ${LIBSDSL_INCLUDE_DIRS}  ${CMAKE_CURRENT_SOURCE_DIR}/fb_wt/wt-fbb-0.1.0)
target_link_libraries(build_wt_dts LINK_PUBLIC ${LIBSDSL_LIBRARIES})


