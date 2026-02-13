//
// Created by Diaz, Diego on 27.11.2025.
//
#include <iostream>
#include <random>

#include "sr-index/include/sr-index/sr_index.h"
#include "sr-index/sri_cli_utils.h"

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
template<class dt_type>
__attribute__((noinline)) void bench_iso_locate(dt_type& dt, std::vector<std::string>& queries, size_t& dummy) {
    for(const auto &query : queries) {
        std::vector<size_t> occ = dt.Locate(query);
        dummy += occ.size();
    }
}

template<class dt_type>
void benchmark_locate(dt_type&dt, const std::string& pat_file) {
    ulint n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);
    //std::cout<<"\tSearching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    size_t dummy = 0;

    //warm up the data structure to avoid page faults
    size_t warmup = std::min<size_t>(200, n_pats);
    for (int k = 0; k < warmup; k++) {
        std::vector<size_t> occ = dt.Locate(pat_list[k]);
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
    std::cout<<"locate dummy: "<<dummy<<std::endl;//print it to avoid optimizations
}

//isolate the function to count the number of LD1 and LD2 cache misses
template<class dt_type>
__attribute__((noinline)) void bench_iso_count(dt_type& dt, std::vector<std::string>& queries, size_t& dummy) {
    for(const auto &query : queries) {
        auto range = dt.Count(query);
        dummy += range.second-range.first+1;
    }
}

template<class dt_type>
void benchmark_count(dt_type&dt, const std::string& pat_file) {
    ulint n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);
    //std::cout<<"\tSearching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    size_t dummy = 0;

    //warm up the data structure to avoid page faults
    size_t warmup = std::min<size_t>(200, n_pats);
    for (int k = 0; k < warmup; k++) {
        auto range = dt.Count(pat_list[k]);
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

    if(argc!=3){
        std::cout<<"usage: ./sri_bench_cmisses <input_dt> <pat_file>"<<std::endl;
        exit(1);
    }

    const auto index_file = std::string(argv[1]);
    const auto pat_file = std::string(argv[2]);
    std::string ext = index_file.substr(index_file.find_last_of(".")+1);

    if(ext=="sri") {
        sri::SrIndex<> sri;
        load_from_file(sri, index_file);
        benchmark_count(sri, pat_file);
        benchmark_locate(sri, pat_file);
    } else if (ext=="sri_vm") {
        sri::SrIndexValidMark<> sri;
        load_from_file(sri, index_file);
        benchmark_count(sri, pat_file);
        benchmark_locate(sri, pat_file);
    } else if (ext=="sri_va") {
        sri::SrIndexValidArea<> sri;
        load_from_file(sri, index_file);
        benchmark_count(sri, pat_file);
        benchmark_locate(sri, pat_file);
    }else {
        std::cout<<"unsupported file format"<<std::endl;
        exit(1);
    }
}