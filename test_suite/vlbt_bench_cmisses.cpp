//
// Created by Diaz, Diego on 25.11.2025.
//
#include <iostream>
#include <ostream>
#include <algorithm>

#include "../include/vlbt_bwt.h"
#include "../include/vlbt_sr_index.h"
#include "../scripts/utils.h"

#include <unordered_set>
#include <random>

//this function discards any data in the LD1 and LD2 caches
void flush_cache() {
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<uint8_t> dist_sym(0, std::numeric_limits<uint8_t>::max());
    static constexpr size_t SIZE = 64*1024*1024; // 64 MB
    volatile uint8_t* buf = new uint8_t[SIZE];
    for (size_t i = 0; i < SIZE; i++) buf[i] = dist_sym(rng);
    delete[] buf;
}

void print_histogram(std::vector<double>& times)  {

    std::sort(times.begin(), times.end());

    const auto q1 = times[times.size() / 4];
    const auto q3 = times[(3 * times.size()) / 4];
    const double iqr = q3 - q1;
    const double upper = q3 + 1.5 * iqr;

    std::vector<double> filtered;
    for (double t : times) {
        if (t <= upper) {
            filtered.push_back(t);
        }
    }

    // Find range
    const double min_val = filtered[0];
    const double max_val = filtered.back();

    // Create histogram bins
    constexpr int num_bins = 20;
    std::vector bins(num_bins, 0);
    const double bin_width = (max_val - min_val) / num_bins;

    for (const double t : filtered) {
        const int bin = std::min(int((t - min_val) / bin_width), num_bins - 1);
        bins[bin]++;
    }

    // Print ASCII histogram
    std::cout << "Histogram of runtimes (microseconds):\n";
    for (int i = 0; i < num_bins; ++i) {
        const double bin_start = min_val + i * bin_width;
        const double bin_end = bin_start + bin_width;
        std::cout << "[" << bin_start << ", " << bin_end << "): ";

        const int count = bins[i];
        for (int j = 0; j < static_cast<int>(count * 50 / filtered.size()); ++j) { // scale to max width 50
            std::cout << '#';
        }
        std::cout << " (" << count << ")\n";
    }
}

std::vector<uint64_t> sample_unique(const uint64_t n, uint64_t x) {
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

template<class dt_type>
std::vector<std::pair<uint64_t, uint8_t>> compute_random_rank_queries(const dt_type& dt, uint64_t n_samps) {

    const size_t text_size = dt.size();
    const size_t alphabet_size = dt.alphabet_size();

    if(n_samps > text_size*alphabet_size) throw std::invalid_argument("invalid sample size");

    std::unordered_set<std::pair<uint64_t, uint8_t>> seen;
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist_idx(0, text_size - 1);
    std::uniform_int_distribution<uint64_t> dist_sym(0, alphabet_size- 1);

    while (seen.size() < n_samps) {
        uint64_t idx = dist_idx(rng);
        uint8_t symbol = dt.eff2byte(dist_sym(rng));
        if (dt.rank(idx, symbol)<0)  continue;
        seen.insert({idx, symbol});
    }
    return {seen.begin(), seen.end()};
}

//isolate the function to count the number of LD1 and LD2 cache misses
template<class dt_type>
__attribute__((noinline)) void bench_iso_locate(dt_type& dt, std::vector<std::string>& queries, size_t& dummy) {
    for(const auto &query : queries) {
        std::vector<uint64_t> occ = dt.locate(query);
        dummy += occ.size();
    }
}
template<class dt_type>
void benchmark_locate(dt_type &dt, const std::string& pat_file) {
    ulint n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);
    size_t dummy = 0;

    //warm up the data structure to avoid page faults
    size_t warmup = std::min<size_t>(200, n_pats);
    for (int k = 0; k < warmup; k++) {
        std::vector<uint64_t> occ = dt.locate(pat_list[k]);
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
    bench_iso_count(dt, pat_list, dummy);
    std::cout<<"locate dummy: "<<dummy<<std::endl;//print it to avoid optimizations
}

//isolate the function to count the number of LD1 and LD2 cache misses
template<class dt_type>
__attribute__((noinline)) void bench_iso_count(dt_type& dt, std::vector<std::string>& queries, size_t& dummy) {
    for(const auto &query : queries) {
        auto range = dt.count(query);
        dummy += range.second-range.first+1;
    }
}
template<class dt_type>
void benchmark_count(dt_type &dt, const std::string& pat_file) {
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

//isolate the function to count the number of LD1 and LD2 cache misses
template<class dt_type>
__attribute__((noinline)) void bench_iso_rank(dt_type& dt, const std::vector<std::pair<uint64_t, uint8_t>>& queries, size_t& dummy) {
    for(const auto &[pos ,sym] : queries) {
        dummy += dt.rank(pos, sym);
    }
}
template<class dt_type>
void benchmark_rank(dt_type &dt, const size_t n) {

    std::vector<std::pair<uint64_t, uint8_t>> queries = compute_random_rank_queries(dt, n);
    size_t dummy = 0;
    //NOTE: no need to warmup because the computation of random queries calls ``dt.rank'' to check the answer is >0
    //pollute the cache so the dt starts cold
    flush_cache();

    //perform the benchmarks
    bench_iso_rank(dt, queries, dummy);
    std::cout<<"rank dummy: "<<dummy<<std::endl;//print it to avoid optimizations
}

//isolate the function to count the number of LD1 and LD2 cache misses
template<class dt_type>
__attribute__((noinline)) void bench_iso_access(dt_type& dt, const std::vector<uint64_t>& queries, size_t& dummy) {
    for(unsigned long long query : queries) {
        dummy += dt[query];
    }
}
template<class dt_type>
void benchmark_access(dt_type &dt, const size_t n) {

    //fill in with random positions
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist(0, n-1);
    std::vector<uint64_t> query_pos(n);
    size_t dummy = 0;
    for (int k = 0; k < n; k++) {
        query_pos[k] = dist(rng);
    }
    //

    //warmup
    size_t warmup = std::min<size_t>(5000, n);
    for (int k = 0; k < warmup; k++) {
        dummy+=dt[dist(rng)];
    }
    //

    //pollute the cache so the dt starts cold
    flush_cache();
    bench_iso_access(dt, query_pos, dummy);
    std::cout<<"access dummy: "<<dummy<<std::endl;//print it to avoid optimizations
}

template<class dt_type>
void bench_rlbwt_int2(const std::string& input_index, const size_t n_samp, const std::string& pat_file) {
    dt_type bwt_th_dt;
    load_from_file(input_index, bwt_th_dt);
    benchmark_access(bwt_th_dt, std::min<size_t>(n_samp, bwt_th_dt.size()));
    benchmark_rank(bwt_th_dt, std::min<size_t>(n_samp, bwt_th_dt.size()));
    benchmark_count(bwt_th_dt, pat_file);
}

template<class dt_type>
void bench_sri_int2(const std::string& input_index, const std::string& pat_file) {
    dt_type bwt_th_dt;
    load_from_file(input_index, bwt_th_dt);
    benchmark_count(bwt_th_dt, pat_file);
    benchmark_locate(bwt_th_dt, pat_file);
}

template<size_t b_size>
void bench_int(const std::string& input_index, const size_t n_samp, VLBT_TYPE& tag, const std::string& pat_file) {

    std::cout<<"Using "<<n_samp<<" samples"<<std::endl;

    switch (tag) {
        case RLBWT:
            std::cout<<"Testing VLBT RLBWT with block size "<<b_size<<std::endl;
            bench_rlbwt_int2<vlbt_rlbwt<b_size>>(input_index, n_samp, pat_file);
            break;
        case RLBWT_WITH_TOEHOLDS:
            std::cout<<"Testing VLBT RLBWT with toeholds and block size "<<b_size<<std::endl;
            bench_rlbwt_int2<vlbt_rlbwt_th<b_size>>(input_index, n_samp, pat_file);
            break;
        case SRI_VALID_AREA:
            std::cout<<"Testing VLBT sr-index with valid area and block size "<<b_size<<std::endl;
            bench_sri_int2<vlbt_sri_va<b_size, b_size>>(input_index, pat_file);
            break;
        default:
            std::cerr<<"Unknown index_type"<<std::endl;
    }
}

int main(int argc, char** argv) {

    if(argc!=4){
        std::cout<<"usage: ./bench_cmiss_vlbt <input_dt> <n_samples> <pat_file>"<<std::endl;
        exit(1);
    }

    const auto input_index = std::string(argv[1]);
    char *pend;
    long int n_samp = strtol(argv[2], &pend, 10);
    const auto pat_file = std::string(argv[3]);

    temp_param_t tp = read_template_param(input_index);
    switch (tp.b_size) {
        case 1024:
            bench_int<1024>(input_index, n_samp, tp.tag, pat_file);
            break;
        case 4096:
            bench_int<4096>(input_index, n_samp, tp.tag, pat_file);
            break;
        case 16384:
            bench_int<16384>(input_index, n_samp, tp.tag, pat_file);
            break;
        case 65536:
            bench_int<65536>(input_index, n_samp, tp.tag, pat_file);
            break;
        case 262144:
            bench_int<262144>(input_index, n_samp, tp.tag, pat_file);
            break;
        case 1048576:
            bench_int<1048576>(input_index, n_samp, tp.tag, pat_file);
            break;
        default:
            exit(1);
    }
}