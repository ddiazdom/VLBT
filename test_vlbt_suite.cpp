//
// Created by Diaz, Diego on 17.10.2022.
//

#include <iostream>
#include <ostream>
//#include "perf_utils.h"
//#include <sdsl/wt_huff.hpp>
//#include <sdsl/construct.hpp>
//#include <sdsl/wt_rlmn.hpp>
//#include <sdsl/wt_int.hpp>
//#include <sdsl/wt_rlmn.hpp>
//#include <sdsl/suffix_arrays.hpp>
//#include <sdsl/suffix_array_algorithm.hpp>
//#include "fb_wt/wt-fbb-0.1.0/wt_fbb.hpp"

//the framework
#include "include/vlbt_build_bwt.h"
#include "include/vlbt_build_phi.h"
#include "include/vlbt_build_sr_index.h"
//#include "include/construct_vlbt_bwt_th.h"

#include "include/vlbt_bwt.h"
#include "include/vlbt_bwt_th.h"
#include "include/vlbt_phi.h"
//#include "include/vlbt_sr_index.h"
//

#include "scripts/fm_index.h"
#include "scripts/custom_wt_rlmn.hpp"
#include <unordered_set>
#include <vector>
#include <random>

using ulint = uint64_t;
void header_error(){
    std::cout << "Error: malformed header in patterns file" << std::endl;
    std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
    exit(0);
}

ulint get_number_of_patterns(std::string header){

    ulint start_pos = header.find("number=");
    if (start_pos == std::string::npos or start_pos+7>=header.size())
        header_error();

    start_pos += 7;

    ulint end_pos = header.substr(start_pos).find(" ");
    if (end_pos == std::string::npos)
        header_error();

    ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());
    return n;
}

ulint get_patterns_length(std::string header){

    ulint start_pos = header.find("length=");
    if (start_pos == std::string::npos or start_pos+7>=header.size())
        header_error();

    start_pos += 7;

    ulint end_pos = header.substr(start_pos).find(" ");
    if (end_pos == std::string::npos)
        header_error();

    ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

    return n;
}

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
sdsl::store_to_file(instance, input_prefix+"."+suffix); \
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
//build_dt(sdsl::wt_huff<>, "wt_huff_bv");\
//build_dt(wt_fbb<sdsl::bit_vector>, "wt_fbb_bv");\
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
#define MEASURE(query, time_answer, query_answer, time_unit) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query_answer = query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<time_unit>( t2 - t1 ).count();\
}

template<class bwt_type>
void test_count(bwt_type& my_dt, std::string& input_file, std::string my_dt_name){

    std::cout<<"Testing count (nanosecs/pat and nanosecs/occ)"<<std::endl;

    sdsl::custom_wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_file+".wt_rlmn");

    std::vector<uint64_t> C(wt_rlmn.sigma+1, 0);
    for(size_t sym=0;sym<wt_rlmn.sigma;sym++){
        C[sym] = wt_rlmn.rank(wt_rlmn.size(), my_dt.eff2byte(sym));
    }

    size_t acc=0, tmp;
    for(size_t i=0;i<wt_rlmn.sigma;i++){
        tmp = C[i];
        C[i] = acc;
        acc+=tmp;
    }
    C[wt_rlmn.sigma] = acc;

    fm_index<sdsl::custom_wt_rlmn<>> csa_rlmn(wt_rlmn, C, my_dt.get_packed_alpha(), my_dt.get_unpacked_alpha());
    //TODO checking for errors
    //std::string pattern = "wart ";
    //csa_rlmn.backward_search(pattern);
    //std::cout<<my_dt.rank(151244695, 'w')<<std::endl;
    //std::cout<<my_dt.rank(151208662, 'w')<<std::endl;
    //std::cout<<wt_rlmn.rank(151208662, 'w')<<std::endl;
    //std::cout<<wt_rlmn.rank(151244695, 'w')<<std::endl;
    //csa_mydt.backward_search(pattern);
    //return;
    //std::cout<<wt_rlmn.rank(228579272, 'b')<<std::endl;
    //exit(1);
    //

    //std::string pat_file = input_file+".pats";
    std::string pat_file = input_file+".pats";
    std::ifstream ifs(pat_file);
    std::string header;
    std::getline(ifs, header);
    ulint n_pats = get_number_of_patterns(header);
    ulint pat_len = get_patterns_length(header);
    std::cout<<"\t"<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    std::vector<std::string> pat_list(n_pats);
    for(ulint i=0;i<n_pats;++i){
        pat_list[i].reserve(pat_len);
        for(ulint j=0;j<pat_len;++j){
            char c;
            ifs.get(c);
            pat_list[i].push_back(c);
        }
    }

    size_t acc_count=0;
    size_t j=0;
    double rlmn_acc_time=0;
    std::vector<std::pair<uint64_t, uint64_t>> rlmn_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(csa_rlmn.backward_search(p), rlmn_acc_time, rlmn_ans[j], std::chrono::nanoseconds)
        acc_count+=rlmn_ans[j].second-rlmn_ans[j].first+1;
        j++;
    }
    //std::cout<<"\tTotal number of occurrences "<<acc_count<<" avg:"<<double(acc_count)/double(n_pats)<<std::endl;

    double my_acc_time=0;
    j=0;
    std::vector<std::pair<uint64_t, uint64_t>> my_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(my_dt.count(p), my_acc_time, my_ans[j], std::chrono::nanoseconds)
        j++;
    }
    std::cout<<"\t"<<my_dt_name<<": ("<<my_acc_time/double(n_pats)<<", "<<my_acc_time/double(acc_count)<<"), ";
    std::cout<<"wt_rlmn: ("<<rlmn_acc_time/double(n_pats)<<", "<<rlmn_acc_time/double(acc_count)<<")"<<std::endl;

    size_t n_errors=0, acc_occ=0;
    for(size_t i=0;i<pat_list.size();i++){
        if(my_ans[i].first!=rlmn_ans[i].first || my_ans[i].second!=rlmn_ans[i].second){
            std::cout<<"Pattern["<<i<<"]: \""<<pat_list[i]<<"\" coords:"<<my_ans[i].first<<"!="<<rlmn_ans[i].first <<" or "<<my_ans[i].second<<"!="<<rlmn_ans[i].second<<std::endl;
            n_errors++;
        }
        acc_occ+=rlmn_ans[i].second-rlmn_ans[i].first+1;
    }
    if(n_errors>0){
        std::cout<<"There are "<<n_errors<<"/"<<pat_list.size()<<" errors "<<std::endl;
    }
    assert(n_errors==0);
}

template<class bwt_type>
void test_locate(bwt_type& my_dt, std::string& input_prefix, std::string my_dt_name){

    std::cout<<"Testing locate (nanosecs/pat and nanosecs/occ)"<<std::endl;
    std::string samp_sa_file = input_prefix+".sa_samples";

    sdsl::custom_wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_prefix+".wt_rlmn");

    std::vector<uint64_t> C(wt_rlmn.sigma+1, 0);
    for(size_t sym=0;sym<wt_rlmn.sigma;sym++){
        C[sym] = wt_rlmn.rank(wt_rlmn.size(), my_dt.eff2byte(sym));
    }

    size_t acc=0, tmp;
    for(size_t i=0;i<wt_rlmn.sigma;i++){
        tmp = C[i];
        C[i] = acc;
        acc+=tmp;
    }
    C[wt_rlmn.sigma] = acc;

    fm_index<sdsl::custom_wt_rlmn<>> csa_rlmn(wt_rlmn, C, my_dt.get_packed_alpha(), my_dt.get_unpacked_alpha());
    std::string pat_file = input_prefix+".pats";
    std::ifstream ifs(pat_file);
    std::string header;
    std::getline(ifs, header);
    ulint n_pats = get_number_of_patterns(header);
    ulint pat_len = get_patterns_length(header);
    std::cout<<"\t"<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    std::vector<std::string> pat_list(n_pats);
    for(ulint i=0;i<n_pats;++i){
        pat_list[i].reserve(pat_len);
        for(ulint j=0;j<pat_len;++j){
            char c;
            ifs.get(c);
            pat_list[i].push_back(c);
        }
    }

    size_t acc_count=0;
    size_t j=0;
    double rlmn_acc_time=0;
    std::vector<std::pair<uint64_t, uint64_t>> rlmn_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(csa_rlmn.backward_search(p), rlmn_acc_time, rlmn_ans[j], std::chrono::nanoseconds)
        acc_count+=rlmn_ans[j].second-rlmn_ans[j].first+1;
        j++;
    }
    //std::cout<<"\tTotal number of occurrences "<<acc_count<<" avg:"<<double(acc_count)/double(n_pats)<<std::endl;

    double my_acc_time=0;
    j=0;
    std::vector<std::pair<uint64_t, uint64_t>> my_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(my_dt.count(p), my_acc_time, my_ans[j], std::chrono::nanoseconds)
        j++;
    }
    std::cout<<"\t"<<my_dt_name<<": ("<<my_acc_time/double(n_pats)<<", "<<my_acc_time/double(acc_count)<<"), ";
    std::cout<<"wt_rlmn: ("<<rlmn_acc_time/double(n_pats)<<", "<<rlmn_acc_time/double(acc_count)<<")"<<std::endl;

    size_t n_errors=0, acc_occ=0;
    for(size_t i=0;i<pat_list.size();i++){
        if(my_ans[i].first!=rlmn_ans[i].first || my_ans[i].second!=rlmn_ans[i].second){
            std::cout<<"Pattern["<<i<<"]: \""<<pat_list[i]<<"\" coords:"<<my_ans[i].first<<"!="<<rlmn_ans[i].first <<" or "<<my_ans[i].second<<"!="<<rlmn_ans[i].second<<std::endl;
            n_errors++;
        }
        acc_occ+=rlmn_ans[i].second-rlmn_ans[i].first+1;
    }
    if(n_errors>0){
        std::cout<<"There are "<<n_errors<<"/"<<pat_list.size()<<" errors "<<std::endl;
    }
    assert(n_errors==0);
}

template<class bwt_type>
void test_access(bwt_type& my_dt, std::string& input_file, std::string my_dt_name){

    std::cout<<"Testing access (avg_time in nanoseconds)"<<std::endl;
    sdsl::custom_wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_file+".wt_rlmn");

    size_t samp_size = 1000000;
    std::vector<uint64_t> samples = sample_unique(wt_rlmn.size(), samp_size);

    double acc_time=0;
    std::vector<uint8_t> my_dt_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(my_dt[samples[j]], acc_time, my_dt_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"\t"<<my_dt_name<<": "<<acc_time/double(samp_size)<<", ";

    acc_time=0;
    std::vector<uint8_t> wt_rlmn_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(wt_rlmn[samples[j]], acc_time, wt_rlmn_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"wt_rlmn: "<<acc_time/double(samp_size)<<std::endl;

    for(size_t j=0;j<samples.size();j++){
        if(wt_rlmn_ans[j]!=my_dt_ans[j]){
            std::cout<<"wt_huff idx: "<<samples[j]<<" |\t sym: "<<int(wt_rlmn_ans[j])<<std::endl;
            std::cout<<"my_dt   idx: "<<samples[j]<<" |\t sym: "<<int(my_dt_ans[j])<<std::endl;
        }
        assert(wt_rlmn_ans[j]==my_dt_ans[j]);
    }
}

template<class bwt_type>
void test_rank(bwt_type& my_dt, std::string& input_file, std::string my_dt_name){

    std::cout<<"Testing rank (avg_time in nanoseconds)"<<std::endl;
    sdsl::custom_wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_file+".wt_rlmn");

    size_t samp_size = 1000000;
    std::vector<std::pair<uint64_t, uint8_t>> tests = compute_random_rank_queries(wt_rlmn.size(), wt_rlmn.sigma, samp_size);

    double acc_time=0;
    std::vector<int64_t> my_dt_ans(samp_size);
    for(size_t j=0;j<tests.size();j++){
        //std::cout<<"rank i:"<<tests[j].first<<" c:"<<my_dt.eff2byte(tests[j].second)<<std::endl;
        MEASURE(my_dt.rank(tests[j].first, my_dt.eff2byte(tests[j].second)), acc_time, my_dt_ans[j], std::chrono::nanoseconds);
    }

    std::cout<<"\t"<<my_dt_name<<": "<<acc_time/double(samp_size)<<", ";

    acc_time=0;
    std::vector<int64_t> wt_rlmn_ans(samp_size);
    for(size_t j=0;j<tests.size();j++){
        MEASURE(wt_rlmn.rank(tests[j].first, my_dt.eff2byte(tests[j].second)), acc_time, wt_rlmn_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"wt_rlmn: "<<acc_time/double(samp_size)<<std::endl;

    for(size_t j=0;j<tests.size();j++){
        if(my_dt_ans[j]<0) continue;
        if(wt_rlmn_ans[j]!=uint64_t(my_dt_ans[j])){
            std::cout<<"query:  idx:"<<tests[j].first<<", sym:"<<int(my_dt.eff2byte(tests[j].second))<<", test_id:"<<j<<std::endl;
            std::cout<<"wt_huff rank answer: "<<wt_rlmn_ans[j]<<std::endl;
            std::cout<<"my_dt   rank answer: "<<my_dt_ans[j]<<"\n"<<std::endl;
        }
        assert(wt_rlmn_ans[j]==uint64_t(my_dt_ans[j]));
    }
}

template<class bwt_type>
void test_inverse_select(bwt_type& my_dt, std::string& input_file, std::string my_dt_name){

    std::cout<<"Testing inverse select (avg_time in nanoseconds)"<<std::endl;
    sdsl::custom_wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_file+".wt_rlmn");

    size_t samp_size = 10000000;
    std::vector<uint64_t> samples = sample_unique(wt_rlmn.size(), samp_size);

    double acc_time=0;
    std::vector<std::pair<uint8_t, uint64_t>> my_dt_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(my_dt.inverse_select(samples[j]), acc_time, my_dt_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"\t"<<my_dt_name<<": "<<acc_time/double(samp_size)<<", ";

    acc_time=0;
    std::vector<std::pair<uint8_t, uint64_t>> wt_rlmn_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(wt_rlmn.inverse_select(samples[j]), acc_time, wt_rlmn_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"wt_rlmn: "<<acc_time/double(samp_size)<<std::endl;

    for(size_t j=0;j<samples.size();j++){
        if(wt_rlmn_ans[j].first!=my_dt_ans[j].first ||
           wt_rlmn_ans[j].second!=my_dt.eff2byte(my_dt_ans[j].second)){
            std::cout<<"wt_huff idx: "<<samples[j]<<" |\t sym: "<<int(wt_rlmn_ans[j].second)<<" rank: "<<int(wt_rlmn_ans[j].first)<<std::endl;
            std::cout<<"my_dt   idx: "<<samples[j]<<" |\t sym: "<<int(my_dt.eff2byte(my_dt_ans[j].second))<<" rank: "<<int(my_dt_ans[j].first)<<"\n"<<std::endl;
        }
        assert(wt_rlmn_ans[j].first==my_dt_ans[j].first &&
               wt_rlmn_ans[j].second==my_dt.eff2byte(my_dt_ans[j].second));
    }
}

void test_bwt(std::string& input_prefix, std::string& output_prefix){

    std::cout<<"Testing VLBT BWT"<<std::endl;
    std::string input_bwt = input_prefix+".ebwt";
    using bwt_type = vlbt_bwt<65536, 64, 4>;
    bwt_type bwt_dt;

    build_bwt<bwt_type>(bwt_dt, input_bwt, BWT_FORMAT::GRL_BWT);
    std::string output_file = output_prefix+".vlbt_bwt";
    size_t written_bytes = store_to_file(output_file, bwt_dt);
    std::cout<<"We store "<<written_bytes<<" in "<<output_file<<"\n"<<std::endl;

    test_rank(bwt_dt, input_prefix, "vlbt_bwt");
    test_count(bwt_dt, input_prefix, "vlbt_bwt");
    test_inverse_select(bwt_dt, input_prefix, "vlbt_bwt");
    test_access(bwt_dt, input_prefix, "vlbt_bwt");
}

template<class sa_samp_type>
void test_bwt_th(std::string& input_prefix, size_t subsamp_val, std::string& output_prefix){

    std::cout<<"Testing VLBT BWT with toeholds"<<std::endl;
    std::string bwt_file = input_prefix+".ebwt";
    using bwt_th_type = vlbt_bwt_th<65536, 64, 4>;
    bwt_th_type bwt_dt;

    std::string samp_sa_file = input_prefix+".sa_samples";
    std::string str_ranges_file = input_prefix+".str_ranges";
    std::string out_ssamp_heads_file = output_prefix+".ssamp_heads";
    std::string out_ssamp_tails_file = output_prefix+".ssamp_tails";
    subsample_sa_samples<sa_samp_type>(samp_sa_file, str_ranges_file, subsamp_val, out_ssamp_heads_file, out_ssamp_tails_file);

    build_bwt_th<bwt_th_type, sa_samp_type>(bwt_dt, bwt_file, BWT_FORMAT::GRL_BWT, out_ssamp_heads_file);
    std::string output_file = output_prefix+".vlbt_bwt_th";
    size_t written_bytes = store_to_file(output_file, bwt_dt);
    std::cout<<"We store "<<written_bytes<<" in "<<output_file<<"\n"<<std::endl;

    test_inverse_select(bwt_dt, input_prefix, "vlbt_bwt_th");
    test_access(bwt_dt, input_prefix, "vlbt_bwt_th");
    test_rank(bwt_dt, input_prefix, "vlbt_bwt_th");
    test_count(bwt_dt, input_prefix, "vlbt_bwt_th");
    test_locate(bwt_dt, input_prefix, "vlbt_bwt_th");
}

template<class size_type>
void test_phi(std::string& input_prefix, size_t ssamp_val, std::string& output_prefix){

    std::string samp_sa_file = input_prefix+".sa_samples";
    std::string str_ranges_file = input_prefix+".str_ranges";

    std::string ssamp_phi_file = output_prefix+".ssamps_phi";
    std::string ssamp_th_file = output_prefix+".ssamps_th";

    subsample_sa_samples<size_type>(samp_sa_file, str_ranges_file, ssamp_val, ssamp_phi_file, ssamp_th_file);
    using phi_type = vlbt_phi<262144, 64, 4>;
    phi_type phi_dt;
    build_phi<phi_type, uint64_t>(phi_dt, ssamp_phi_file);

    std::string output_file = output_prefix+".vlbt_phi";
    size_t written_bytes = store_to_file(output_file, phi_dt);
    std::cout<<"We store "<<written_bytes<<" in "<<output_file<<std::endl;
}

template<class size_type>
void test_sr_index(std::string& input_prefix, size_t ssamp_val, std::string& output_prefix){

    using bwt_th_type = vlbt_bwt_th<65536, 64, 4>;
    using phi_type = vlbt_phi<65536, 64, 4>;
    using sr_index_type = vlbt_sr_index<bwt_th_type, phi_type>;
    sr_index_type sr_index;
    build_sr_index<sr_index_type , size_type>(sr_index, input_prefix, ssamp_val, output_prefix);

    std::string output_sr_index_file = output_prefix+".sr_index";
    size_t written_bytes = store_to_file(output_sr_index_file, sr_index);
    std::cout<<"Final sr-index uses "<<written_bytes<<" bytes ("<< double(written_bytes*8)/double(sr_index.size())<<" bps)"<<std::endl;

    //test_count(sr_index, input_prefix, "sr_index");
    //test_inverse_select(sr_index.bwt_with_th, input_prefix, "sr_index");
    //test_rank(sr_index.bwt_with_th, input_prefix, "sr_index");
    //test_access(sr_index.bwt_with_th, input_prefix, "sr_index");
}

int main(int argc, char** argv){

    if(argc!=4){
        std::cout<<"usage: ./test_vlbt input_prefix build_other_dts=0|1 output_prefix"<<std::endl;
        exit(1);
    }

    std::string input_prefix = std::string(argv[1]);
    char *pend;
    long int other_dts = strtol(argv[2], &pend, 10);
    assert(other_dts>=0 && other_dts<=1);
    std::string output_prefix = std::string(argv[3]);

    if(other_dts){
        //std::string plain_input_file = "tmp_plain.txt";

        std::string bwt_file = input_prefix+".ebwt";
        std::string wt_file = input_prefix+".wt_rlmn";
        std::cout<<"Creating wavelet trees for "<<bwt_file<<std::endl;
        sdsl::custom_wt_rlmn<> wt(bwt_file);
        sdsl::store_to_file(wt, wt_file);
        //rl2plain(bwt_file, plain_input_file);
        //TESTED_DTS
    }
    test_bwt(input_prefix, output_prefix);
    //test_bwt_th<uint64_t>(input_prefix, 4, output_prefix);
    //test_phi<uint64_t>(input_prefix, 4, output_prefix);
    //test_sr_index<uint64_t>(input_prefix, 4, output_prefix);
}