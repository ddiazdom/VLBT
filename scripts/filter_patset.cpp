/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */
#include <iostream>
#include <fstream>
#include "utils.h"
#include "../include/vlbt_bwt.h"

template<class dt_type>
void filter_int3(const std::string& input_index, const std::string& pat_file,
                 size_t max_freq, std::string output_pat_file){
    dt_type dt;
    load_from_file(input_index, dt);

    uint64_t n_pats, pat_len, n_del_pats=0;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    for(auto &p : pat_list) {
        std::pair<uint64_t, uint64_t> ans = dt.count(p);
        size_t count=ans.second-ans.first+1;
        if (count>max_freq) {
            p.clear();
            n_del_pats++;
        }
    }

    std::string header = "# number="+std::to_string(n_pats-n_del_pats)+" length="+std::to_string(pat_len)+" file=missing forbidden=missing\n";
    std::ofstream ofs(output_pat_file, std::ios::out);
    ofs.write(header.data(), static_cast<std::streamsize>(header.size()));

    for (auto &p : pat_list) {
        if (!p.empty()) {
            ofs.write(p.data(), pat_len);
        }
    }

    std::cout<<"We kept "<<n_pats-n_del_pats<<" patterns of length "<<pat_len<<" from file "<<pat_file<<std::endl;
}

template<size_t b_size>
void filter_int2(temp_param_t& tp, const std::string& input_index, const std::string& pat_file,
                 size_t max_freq, std::string output_pat_file) {
    switch (tp.tag) {
        case RLBWT:
            filter_int3<vlbt_rlbwt<b_size>>(input_index, pat_file, max_freq, output_pat_file);
            break;
        case RLBWT_WITH_TOEHOLDS:
            filter_int3<vlbt_rlbwt_th<b_size>>(input_index, pat_file, max_freq, output_pat_file);
            break;
        default:
            exit(1);
    }
}

void filter_int(const std::string& input_index, const std::string& pat_file,
                size_t max_freq, std::string output_pat_file) {

    temp_param_t tp = read_template_param(input_index);
    switch (tp.b_size) {
        case 1024:
            filter_int2<1024>(tp, input_index, pat_file, max_freq, output_pat_file);
            break;
        case 4096:
            filter_int2<4096>(tp, input_index, pat_file, max_freq, output_pat_file);
            break;
        case 16384:
            filter_int2<16384>(tp, input_index, pat_file, max_freq,  output_pat_file);
            break;
        case 65536:
            filter_int2<65536>(tp, input_index, pat_file, max_freq, output_pat_file);
            break;
        case 262144:
            filter_int2<262144>(tp, input_index, pat_file, max_freq, output_pat_file);
            break;
        case 1048576:
            filter_int2<1048576>(tp, input_index, pat_file, max_freq, output_pat_file);
            break;
        default:
            exit(1);
    }
}

int main(int argc, char** argv){

    if(argc!=5){
        std::cout<<"usage: ./filter_patset <input_pat_file> <input_rl_bwt> <max_freq> <output_filtered_pat_file>"<<std::endl;
        exit(1);
    }

    auto input_pat_file = std::string(argv[1]);
    auto rlbwt_index_file = std::string(argv[2]);

    char *pend;
    long int max_freq = strtol(argv[3], &pend, 10);
    auto output_pat_file = std::string(argv[4]);

    filter_int(rlbwt_index_file, input_pat_file, max_freq, output_pat_file);
}