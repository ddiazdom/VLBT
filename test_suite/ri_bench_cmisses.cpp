//
// Created by Diaz, Diego on 27.11.2025.
//
#include <iostream>
#include "r-index/internal/r_index.hpp"
#include "../scripts/utils.h"

//this function discards any data in the LD1 and LD2 caches
void flush_cache() {
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<uint8_t> dist_sym(0, std::numeric_limits<uint8_t>::max());
    static constexpr size_t SIZE = 128*1024*1024; // 64 MB
    volatile uint8_t* buf = new uint8_t[SIZE];
    for (size_t i = 0; i < SIZE; i++) buf[i] = dist_sym(rng);
    delete[] buf;
}

//isolate the function to count the number of LD1 and LD2 cache misses
__attribute__((noinline)) void bench_iso_locate(ri::r_index<>& dt, std::vector<std::string>& queries, size_t& dummy) {
    for(const auto &query : queries) {
        std::vector<ulint> occ = dt.locate_all(query);
        dummy += occ.size();
    }
}
void benchmark_locate(ri::r_index<>&dt, const std::string& pat_file) {
    ulint n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);
    //std::cout<<"\tSearching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    size_t dummy = 0;

    //warm up the data structure to avoid page faults
    size_t warmup = std::min<size_t>(200, n_pats);
    for (int k = 0; k < warmup; k++) {
        std::vector<ulint> occ = dt.locate_all(pat_list[k]);
        dummy+=occ.size();
    }
    //

    //shuffle the list of patterns to prevent bias
    auto rd = std::random_device {};
    auto rng = std::default_random_engine {rd()};
    std::shuffle(pat_list.begin(), pat_list.end(), rng);
    //

    //pollute LD1 and LD2 to evict the dt and start the analysis in cold
    flush_cache();

    //perform the benchmark
    bench_iso_locate(dt, pat_list, dummy);

    std::cout<<"count dummy: "<<dummy<<std::endl;//print it to avoid optimizations
}

//isolate the function to count the number of LD1 and LD2 cache misses
__attribute__((noinline)) void bench_iso_count(ri::r_index<>& dt, std::vector<std::string>& queries, size_t& dummy) {
    for(const auto &query : queries) {
        auto range = dt.count(query);
        dummy += range.second-range.first+1;
    }
}
void benchmark_count(ri::r_index<>&dt, const std::string& pat_file) {
    ulint n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);
    //std::cout<<"\tSearching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    size_t dummy = 0;

    //warm up the data structure to avoid page faults
    size_t warmup = std::min<size_t>(200, n_pats);
    for (int k = 0; k < warmup; k++) {
        auto range = dt.count(pat_list[k]);
        dummy+=range.second-range.first+1;
    }
    //

    //shuffle the list of patterns to prevent bias
    auto rd = std::random_device {};
    auto rng = std::default_random_engine {rd()};
    std::shuffle(pat_list.begin(), pat_list.end(), rng);
    //

    //pollute LD1 and LD2 to evict the dt and start the analysis in cold
    flush_cache();

    //perform the benchmark
    bench_iso_count(dt, pat_list, dummy);

    std::cout<<"count dummy: "<<dummy<<std::endl;//print it to avoid optimizations
}

int main(int argc, char** argv) {

    if(argc!=5){
        std::cout<<"usage: ./ri_bench_cmisses <input_dt> <pat_file>"<<std::endl;
        exit(1);
    }

    const auto input_index = std::string(argv[1]);
    const auto pat_file = std::string(argv[2]);

    bool fast=false;
    std::ifstream in(input_index);
    in.read((char*)&fast,sizeof(fast));
    ri::r_index<> ri;
    ri.load(in);
    in.close();

    benchmark_count(ri, pat_file);
    benchmark_locate(ri, pat_file);
}