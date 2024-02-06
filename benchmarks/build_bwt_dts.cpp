//
// Created by Diaz, Diego on 17.10.2022.
//

#include <iostream>
#include <ostream>

#include <sdsl/wt_huff.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/wt_blcd.hpp>
#include <sdsl/wt_int.hpp>
#include <sdsl/wt_rlmn.hpp>
#include "fb_wt/wt-fbb-0.1.0/wt_fbb.hpp"
#include "../simple_rl_bwt.h"

#define build_dt(dt, suffix) \
{                            \
std::cout<<"Building "<<suffix<<std::endl;\
dt instance;\
sdsl::construct(instance, plain_input_file, 1);\
sdsl::store_to_file(instance, output_file+"."+suffix); \
}\

/*
#define TESTED_DTS \
build_dt(sdsl::wt_huff<>, "wt_huff_bv");\
build_dt(sdsl::wt_huff<sdsl::rrr_vector<>>, "wt_huff_rrr");\
build_dt(sdsl::wt_huff<sdsl::hyb_vector<>>, "wt_huff_hyb");\
build_dt(sdsl::wt_huff<sdsl::bit_vector_il<>>, "wt_huff_il");\
build_dt(sdsl::wt_blcd<>, "wt_blcd_bv");\
build_dt(sdsl::wt_blcd<sdsl::rrr_vector<>>, "wt_blcd_rrr");\
build_dt(sdsl::wt_blcd<sdsl::hyb_vector<>>, "wt_blcd_hyb");\
build_dt(sdsl::wt_blcd<sdsl::bit_vector_il<>>, "wt_blcd_il");\
build_dt(sdsl::wt_rlmn<>, "rlmn");      \
build_dt(wt_fbb<>, "wt_fbb");           \
*/

/*#define TESTED_DTS \
build_dt(wt_fbb<sdsl::bit_vector>, "wt_fbb_bv"); \
build_dt(wt_fbb<sdsl::rrr_vector<>>, "wt_fbb_rrr"); \
build_dt(wt_fbb<sdsl::hyb_vector<>>, "wt_fbb_hyb"); \
build_dt(wt_fbb<sdsl::bit_vector_il<>>, "wt_fbb_il"); \
build_dt(sdsl::wt_huff<>, "wt_huff_bv");\
build_dt(sdsl::rlmn<>, "wt_rlmn");               \*/
                   \
#define TESTED_DTS \
build_dt(wt_fbb<sdsl::bit_vector>, "wt_fbb_bv"); \
build_dt(wt_fbb<sdsl::rrr_vector<>>, "wt_fbb_rrr"); \
build_dt(wt_fbb<sdsl::hyb_vector<>>, "wt_fbb_hyb"); \
build_dt(wt_fbb<sdsl::bit_vector_il<>>, "wt_fbb_il"); \
build_dt(sdsl::wt_huff<>, "wt_huff_bv");\
build_dt(sdsl::wt_rlmn<>, "wt_rlmn");\

void rl2plain(std::string& rl_file, std::string& output_plain_file){

    std::ofstream ofs(output_plain_file, std::ios::out | std::ios::binary);
    uint8_t buffer[1024]={0};
    bwt_buff_reader bwt_reader(rl_file);
    size_t sym, freq, k=0;
    size_t sym_freqs[256]={0};
    for(size_t i=0;i<bwt_reader.size();i++){
        bwt_reader.read_run(i, sym, freq);
        for(size_t j=0;j<freq;j++){
            buffer[k++] = sym;
            if(k==1024){
                ofs.write((char *)buffer, 1024);
                k=0;
            }
        }
        sym_freqs[sym]+=freq;
    }
    if(k!=0){
        ofs.write((char *)buffer, (std::streamsize)k);
    }
    bwt_reader.close();
    ofs.close();
}

int main(int argc, char** argv){

    if(argc!=3){
        std::cout<<"usage: ./build-wt-huff-bwt plain_bwt.txt output_file"<<std::endl;
        exit(1);
    }

    std::string input_file = std::string(argv[1]);
    std::string output_file = std::string(argv[2]);

    std::string plain_input_file = "tmp_plain.txt";
    rl2plain(input_file, plain_input_file);

    std::cout<<"Creating wavelet trees for "<<input_file<<std::endl;
    TESTED_DTS

    std::cout<<"Creating my bwt "<<input_file<<std::endl;
    simple_rl_bwt<5> my_bwt(input_file);
    store_to_file(output_file+".my_simple_bwt", my_bwt);
}