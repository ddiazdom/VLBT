/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */
#include "../include/bwt_streams.h"
#include <filesystem>

int main(int argc, char** argv){

    if(argc!=3){
        std::cout<<"usage: ./plain2rle <input_plain_bwt> <output_vbyte_rle_bwt>"<<std::endl;
        exit(1);
    }

    auto plain_bwt = std::string(argv[1]);
    auto grlbwt_fmt_bwt = std::filesystem::path(std::string(argv[2])).replace_extension(".vbyte_rle");
    plain2vbyte_rle(plain_bwt, grlbwt_fmt_bwt);
}