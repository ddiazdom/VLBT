//
// Created by Diaz, Diego on 17.10.2022.
//

#include <iostream>
#include <ostream>
#include "perf_utils.h"
#include <sdsl/wt_huff.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/wt_blcd.hpp>
#include <sdsl/wt_int.hpp>
#include <sdsl/wt_rlmn.hpp>
#include "fb_wt/wt-fbb-0.1.0/wt_fbb.hpp"

#include "../rlbwt_small_alpha.h"
#include "../rlbwt_vlb.h"

#define build_dt(dt, suffix) \
{                            \
dt instance;                 \
auto t1 = std::chrono::high_resolution_clock::now();\
sdsl::construct(instance, plain_input_file, 1);\
auto t2 = std::chrono::high_resolution_clock::now();\
sdsl::store_to_file(instance, output_prefix+"."+suffix); \
std::cout<<suffix<<" "<<report_time(t1, t2, 0)<<",  space_usage: "<<double(sdsl::size_in_bytes(instance))/1000000<<" MB ("<<float(sdsl::size_in_bytes(instance)*8)/float(instance.size())<<" bps)"<<std::endl;\
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
build_dt(sdsl::wt_rlmn<>, "wt_rlmn");\
build_dt(sdsl::wt_huff<>, "wt_huff_bv");\
build_dt(wt_fbb<sdsl::bit_vector>, "wt_fbb_bv");\
build_dt(wt_fbb<sdsl::rrr_vector<>>, "wt_fbb_rrr");\
build_dt(wt_fbb<sdsl::hyb_vector<>>, "wt_fbb_hyb");\
//build_dt(wt_fbb<sdsl::bit_vector_il<>>, "wt_fbb_il");\

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

/*template<class bwt_type>
void build_my_bwt(std::string& input_file, std::string& output_file){
    auto t1 = std::chrono::high_resolution_clock::now();
    bwt_type my_bwt(input_file);
    auto t2 = std::chrono::high_resolution_clock::now();
    size_t written_bytes = store_to_file(output_file+".my_simple_bwt", my_bwt);

    int status;
    const char * dt_name = typeid(bwt_type).name();
    char* dem_name = abi::__cxa_demangle(dt_name, nullptr, nullptr, &status);

    std::cout<<dem_name<<"  build_time:"<<report_time(t1, t2, 0)<<",  space_usage:"<<float(written_bytes*8)/float(my_bwt.size())<<" bps"<<std::endl;

    //todo testing
    / *std::vector<uint8_t> cs(6,0);
    uint64_t k;
    std::vector<uint64_t> rank_c_i(6,0);
    std::vector<uint64_t> rank_c_j(6, 0);
    my_bwt.interval_symbols(8756, 23065, k, cs, rank_c_i, rank_c_j);

    rlbwt_small_alpha<sigma> tmp_bwt;
    load_from_file(output_file+".my_simple_bwt", tmp_bwt);
    std::vector<uint8_t> cs2(6,0);
    uint64_t k2;
    std::vector<uint64_t> rank_c_i2(6,0);
    std::vector<uint64_t> rank_c_j2(6, 0);
    std::cout<<sigma<<std::endl;
    tmp_bwt.interval_symbols(8756, 23065, k2, cs2, rank_c_i2, rank_c_j2); * /
    //
    my_bwt.stats();

    switch (alphabet) {
        case 3:
            build_my_bwt<rlbwt_small_alpha<3>>(input_file, output_file);
            break;
        case 4:
            build_my_bwt<rlbwt_small_alpha<4>>(input_file, output_file);
            break;
        case 5:
            build_my_bwt<rlbwt_small_alpha<5>>(input_file, output_file);
            break;
        case 6:
            build_my_bwt<rlbwt_small_alpha<6>>(input_file, output_file);
            break;
        case 7:
            build_my_bwt<rlbwt_small_alpha<7>>(input_file, output_file);
            break;
        case 8:
            build_my_bwt<rlbwt_small_alpha<8>>(input_file, output_file);
            break;
        case 9:
            build_my_bwt<rlbwt_small_alpha<9>>(input_file, output_file);
            break;
        case 10:
            build_my_bwt<rlbwt_small_alpha<10>>(input_file, output_file);
            break;
        case 11:
            build_my_bwt<rlbwt_small_alpha<11>>(input_file, output_file);
            break;
        case 12:
            build_my_bwt<rlbwt_small_alpha<12>>(input_file, output_file);
            break;
        case 13:
            build_my_bwt<rlbwt_small_alpha<13>>(input_file, output_file);
            break;
        case 14:
            build_my_bwt<rlbwt_small_alpha<14>>(input_file, output_file);
            break;
        case 15:
            build_my_bwt<rlbwt_small_alpha<15>>(input_file, output_file);
            break;
        case 16:
            build_my_bwt<rlbwt_small_alpha<16>>(input_file, output_file);
            break;
        default:
            std::cout<<"Alphabet size not supported"<<std::endl;
            exit(1);
    }
}*/

int main(int argc, char** argv){

    if(argc!=4){
        std::cout<<"usage: ./build_bwt_dts input_bwt.rl_bwt build_other_dts=0|1 output_prefix"<<std::endl;
        exit(1);
    }

    std::string input_file = std::string(argv[1]);
    char *pend;
    long int other_dts = strtol(argv[2], &pend, 10);
    assert(other_dts>=0 && other_dts<=1);

    if(other_dts){
        std::cout<<"We will build other DTs for the BWT..."<<std::endl;
        std::string output_prefix = std::string(argv[3]);
        std::string plain_input_file = "tmp_plain.txt";
        rl2plain(input_file, plain_input_file);
        std::cout<<"Creating wavelet trees for "<<input_file<<std::endl;
        TESTED_DTS
    }else{
        std::cout<<"We will build only my BWT..."<<std::endl;
    }

    using bwt_type = rlbwt_vlb<65536, 64, 4>;
    bwt_type bwt;
    build_rlbwt_vlb<bwt_type>(bwt, input_file, INPUT_FORMAT::GRL_BWT);
}