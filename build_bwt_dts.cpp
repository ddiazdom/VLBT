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

#include "../../rlbwt_small_alpha.h"
#include "../../construct_rlbwt_vlb.h"
#include "../../rlbwt_vlb.h"

#include <unordered_set>
#include <vector>
#include <random>

std::vector<uint64_t> sample_unique(uint64_t n, uint64_t x) {
    if (x > n) throw std::invalid_argument("x cannot be larger than n");

    std::unordered_set<uint64_t> seen;
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist(0, n - 1);

    while (seen.size() < x) {
        uint64_t val = dist(rng);
        seen.insert(val);  // insert is a no-op if val already exists
    }

    return {seen.begin(), seen.end()};
}

namespace std {
    template<typename X, typename Y>
    struct hash<std::pair<X, Y>> {
        std::size_t operator()(const std::pair<X, Y> &pair) const {
            return std::hash<X>()(pair.first) ^ std::hash<Y>()(pair.second);
        }
    };
}

std::vector<std::pair<uint64_t, uint8_t>> compute_random_rank_queries(uint64_t text_size, uint8_t alphabet_size, uint64_t n_samps) {

    if(n_samps > text_size*alphabet_size) throw std::invalid_argument("invalid sample size");

    std::unordered_set<std::pair<uint64_t, uint8_t>> seen;
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist_idx(0, text_size - 1);
    std::uniform_int_distribution<uint64_t> dist_sym(0, alphabet_size- 1);

    while (seen.size() < n_samps) {
        uint64_t idx = dist_idx(rng);
        uint64_t symbol = dist_sym(rng);
        seen.insert({idx, symbol});
    }
    return {seen.begin(), seen.end()};
}

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
//build_dt(wt_fbb<sdsl::rrr_vector<>>, "wt_fbb_rrr");\
//build_dt(wt_fbb<sdsl::hyb_vector<>>, "wt_fbb_hyb");\
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

#define MEASURE(query, time_answer, query_answer) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query_answer = query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();\
}

template<class bwt_type>
void test_access(bwt_type& my_dt, std::string& input_file){
    //size_t samp_size = (wt_dt.size()*10)/100;

    sdsl::wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_file+".wt_rlmn");
    /*for(size_t i=0;i<wt_rlmn.size();i++){
        assert(wt_rlmn[i]==my_dt[i]);
    }*/

    size_t samp_size = 1000000;
    std::vector<uint64_t> samples = sample_unique(wt_rlmn.size(), samp_size);

    double acc_time=0;
    std::vector<uint8_t> my_dt_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(my_dt[samples[j]], acc_time, my_dt_ans[j]);
    }
    std::cout<<"access rlbwt_vlb:";
    std::cout<<acc_time/double(samp_size)<<" nanoseconds"<<std::endl;

    /*acc_time=0;
    std::vector<uint8_t> wt_huff_ans(samp_size);
    sdsl::wt_huff<> wt_huff;
    sdsl::load_from_file(wt_huff, input_file+".wt_huff_bv");
    for(size_t j=0;j<samples.size();j++){
        MEASURE(wt_huff[samples[j]], acc_time, wt_huff_ans[j]);
    }
    std::cout<<"access wt_huff:";
    std::cout<<acc_time/double(samp_size)<<" nanoseconds"<<std::endl;*/

    acc_time=0;
    std::vector<uint8_t> wt_rlmn_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(wt_rlmn[samples[j]], acc_time, wt_rlmn_ans[j]);
    }
    std::cout<<"access wt_rlmn:";
    std::cout<<acc_time/double(samp_size)<<" nanoseconds"<<std::endl;

    for(size_t j=0;j<samples.size();j++){
        if(wt_rlmn_ans[j]!=my_dt_ans[j]){
            std::cout<<"wt_huff idx: "<<samples[j]<<" |\t sym: "<<int(wt_rlmn_ans[j])<<std::endl;
            std::cout<<"my_dt   idx: "<<samples[j]<<" |\t sym: "<<int(my_dt_ans[j])<<std::endl;
        }
        assert(wt_rlmn_ans[j]==my_dt_ans[j]);
    }

    /*for(size_t j=0;j<my_dt.size();j++){
        //size_t p = rand() % my_dt.size();
        size_t p = j;

        auto res1= wt_huff.inverse_select(p);
        auto res2 = my_dt.inverse_select(p);

        if(res1.first!=res2.first || res1.second!=my_dt.eff2byte(res2.second)){
            std::cout<<"wt_huff idx: "<<p<<" |\t sym: "<<int(res1.second)<<" rank: "<<res1.first<<std::endl;
            std::cout<<"my_dt   idx: "<<p<<" |\t sym: "<<int(my_dt.eff2byte(res2.second))<<" rank: "<<res2.first<<"\n"<<std::endl;
        }
        assert(res1.first==res2.first && res1.second==my_dt.eff2byte(res2.second));
        //MEASURE(my_dt.inverse_select(p), acc_time, answers[0]);
    }*/
}

template<class bwt_type>
void test_rank(bwt_type& my_dt, std::string& input_file){

    sdsl::wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_file+".wt_rlmn");

    size_t samp_size = 1000000;
    std::vector<std::pair<uint64_t, uint8_t>> tests = compute_random_rank_queries(wt_rlmn.size(), wt_rlmn.sigma, samp_size);

    double acc_time=0;
    std::vector<int64_t> my_dt_ans(samp_size);
    for(size_t j=0;j<tests.size();j++){
        //std::cout<<tests[j].first<<" "<<int(tests[j].second)<<std::endl;
        MEASURE(my_dt.rank(tests[j].first, tests[j].second), acc_time, my_dt_ans[j]);
    }
    std::cout<<"rank rlbwt_vlb:";
    std::cout<<acc_time/double(samp_size)<<" nanoseconds"<<std::endl;

    acc_time=0;
    std::vector<int64_t> wt_rlmn_ans(samp_size);
    for(size_t j=0;j<tests.size();j++){
        MEASURE(wt_rlmn.rank(tests[j].first, my_dt.eff2byte(tests[j].second)), acc_time, wt_rlmn_ans[j]);
    }
    std::cout<<"rank wt_rlmn:";
    std::cout<<acc_time/double(samp_size)<<" nanoseconds"<<std::endl;

    for(size_t j=0;j<tests.size();j++){
        if(my_dt_ans[j]<0) continue;
        if(wt_rlmn_ans[j]!=uint64_t(my_dt_ans[j])){
            std::cout<<"query:  idx:"<<tests[j].first<<", sym:"<<int(tests[j].second)<<std::endl;
            std::cout<<"wt_huff rank answer: "<<wt_rlmn_ans[j]<<std::endl;
            std::cout<<"my_dt   rank answer: "<<my_dt_ans[j]<<"\n"<<std::endl;
        }
        assert(wt_rlmn_ans[j]==uint64_t(my_dt_ans[j]));
    }
}

template<class bwt_type>
void test_inverse_select(bwt_type& my_dt, std::string& input_file){
    //size_t samp_size = (wt_dt.size()*10)/100;

    sdsl::wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_file+".wt_rlmn");

    size_t samp_size = 1000000;
    std::vector<uint64_t> samples = sample_unique(wt_rlmn.size(), samp_size);

    double acc_time=0;
    std::vector<std::pair<uint8_t, uint64_t>> my_dt_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(my_dt.inverse_select(samples[j]), acc_time, my_dt_ans[j]);
    }
    std::cout<<"inverse_select rlbwt_vlb:";
    std::cout<<acc_time/double(samp_size)<<" nanoseconds"<<std::endl;

    acc_time=0;
    std::vector<std::pair<uint8_t, uint64_t>> wt_rlmn_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(wt_rlmn.inverse_select(samples[j]), acc_time, wt_rlmn_ans[j]);
    }
    std::cout<<"inverse_select wt_rlmn:";
    std::cout<<acc_time/double(samp_size)<<" nanoseconds"<<std::endl;

    //acc_time=0;
    //sdsl::wt_huff<> wt_huff;
    //sdsl::load_from_file(wt_huff, input_file+".wt_huff_bv");
    //std::vector<std::pair<uint8_t, uint64_t>> wt_huff_ans(samp_size);
    //for(size_t j=0;j<samples.size();j++){
    //    MEASURE(wt_huff.inverse_select(samples[j]), acc_time, wt_huff_ans[j]);
    //}
    //std::cout<<"inverse_select wt_huff:";
    //std::cout<<acc_time/double(samp_size)<<" nanoseconds"<<std::endl;

    for(size_t j=0;j<samples.size();j++){
        if(wt_rlmn_ans[j].first!=my_dt_ans[j].first ||
           wt_rlmn_ans[j].second!=my_dt.eff2byte(my_dt_ans[j].second)){
            std::cout<<"wt_huff idx: "<<samples[j]<<" |\t sym: "<<int(wt_rlmn_ans[j].second)<<" rank: "<<int(wt_rlmn_ans[j].first)<<std::endl;
            std::cout<<"my_dt   idx: "<<samples[j]<<" |\t sym: "<<int(my_dt.eff2byte(my_dt_ans[j].second))<<" rank: "<<int(my_dt_ans[j].first)<<"\n"<<std::endl;
        }
        assert(wt_rlmn_ans[j].first==my_dt_ans[j].first &&
               wt_rlmn_ans[j].second==my_dt.eff2byte(my_dt_ans[j].second));
    }

    /*for(size_t j=0;j<my_dt.size();j++){
        //size_t p = rand() % my_dt.size();
        size_t p = j;

        auto res1= wt_huff.inverse_select(p);
        auto res2 = my_dt.inverse_select(p);

        if(res1.first!=res2.first || res1.second!=my_dt.eff2byte(res2.second)){
            std::cout<<"wt_huff idx: "<<p<<" |\t sym: "<<int(res1.second)<<" rank: "<<res1.first<<std::endl;
            std::cout<<"my_dt   idx: "<<p<<" |\t sym: "<<int(my_dt.eff2byte(res2.second))<<" rank: "<<res2.first<<"\n"<<std::endl;
        }
        assert(res1.first==res2.first && res1.second==my_dt.eff2byte(res2.second));
        //MEASURE(my_dt.inverse_select(p), acc_time, answers[0]);
    }*/
}

int main(int argc, char** argv){

    if(argc!=4){
        std::cout<<"usage: ./build_bwt_dts input_bwt.rl_bwt build_other_dts=0|1 output_prefix"<<std::endl;
        exit(1);
    }

    std::string input_file = std::string(argv[1]);
    char *pend;
    long int other_dts = strtol(argv[2], &pend, 10);
    assert(other_dts>=0 && other_dts<=1);
    std::string output_prefix = std::string(argv[3]);

    if(other_dts){
        std::cout<<"We will build other DTs for the BWT..."<<std::endl;
        std::string plain_input_file = "tmp_plain.txt";
        rl2plain(input_file, plain_input_file);
        std::cout<<"Creating wavelet trees for "<<input_file<<std::endl;
        TESTED_DTS
    }else{
        std::cout<<"We will build only my BWT..."<<std::endl;
    }

    using bwt_type = rlbwt_vlb<4096, 64, 4>;
    bwt_type bwt_dt;
    build_rlbwt_vlb<bwt_type>(bwt_dt, input_file, INPUT_FORMAT::GRL_BWT);
    std::string output_file = output_prefix+".rlbwt_vlb";
    size_t written_bytes = store_to_file(output_file, bwt_dt);
    std::cout<<"We store "<<written_bytes<<" in "<<output_file<<std::endl;

    //std::cout<<bwt_dt.rank(250020433,126)<<std::endl;
    //sdsl::wt_rlmn<> wt_rlmn;
    //sdsl::load_from_file(wt_rlmn, output_prefix+".wt_rlmn");
    //std::cout<<bwt_dt.rank(467616716, 57)<<std::endl;
    //std::cout<<wt_rlmn.rank(467616716,bwt_dt.eff2byte(57))<<std::endl;

    test_inverse_select(bwt_dt, output_prefix);
    test_access(bwt_dt, output_prefix);//not implemented in SSE4.2 or AVX2
    test_rank(bwt_dt, output_prefix);

    /*uint64_t r1;
    for(size_t i=0;i<wt_rlmn.sigma;i++){
        r1 = wt_rlmn.rank(4200, bwt_dt.eff2byte(i));
        std::cout<<r1<<std::endl;
    }*/
    //auto r1 = wt_rlmn.rank(4200, bwt_dt.eff2byte(2));
    //auto r2 = bwt_dt.rank(4200, 2);
    //std::cout<<r1<<" "<<r2<<std::endl;
    //std::cout<<wt_rlmn.rank(137527296, bwt_dt.eff2byte(1))<<std::endl;
    //bwt_dt.rank(122589194, 1);
    //bwt_dt.rank(10, 1);
}