//
// Created by Diaz, Diego on 17.10.2022.
//

#include <iostream>
#include <ostream>
#include <cxxabi.h>
#include <memory>

/*#include <sdsl/wt_huff.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/wt_blcd.hpp>
#include <sdsl/wt_int.hpp>
#include <sdsl/wt_rlmn.hpp>
#include "fb_wt/wt-fbb-0.1.0/wt_fbb.hpp"*/

#include "../rlbwt_small_alpha.h"
#include "../rlbwt_vlb.h"

#include "performancecounters/benchmarker.h"
#include <cassert>

void pretty_print(size_t volume, size_t bytes, std::string name,
                  event_aggregate agg) {
    printf("%-40s : ", name.c_str());
    printf(" %5.2f GB/s ", bytes / agg.fastest_elapsed_ns());
    printf(" %5.1f Ma/s ", volume * 1000.0 / agg.fastest_elapsed_ns());
    printf(" %5.2f ns/d ", agg.fastest_elapsed_ns() / volume);
    if (collector.has_events()) {
        printf(" %5.2f GHz ", agg.fastest_cycles() / agg.fastest_elapsed_ns());
        printf(" %5.2f c/d ", agg.fastest_cycles() / volume);
        printf(" %5.2f i/d ", agg.fastest_instructions() / volume);
        printf(" %5.2f c/b ", agg.fastest_cycles() / bytes);
        printf(" %5.2f i/b ", agg.fastest_instructions() / bytes);
        printf(" %5.2f i/c ", agg.fastest_instructions() / agg.fastest_cycles());
        printf(" %5.2f br_miss ", agg.branch_misses());
    }
    printf("\n");
}

/*template<class time_t>
std::string report_time(time_t start, time_t end, size_t padding){
    auto dur = end - start;
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);
    std::stringstream time;

    for(size_t i=0;i<padding;i++) std::cout<<" ";
    if(h.count()>0){
        time <<"build_time (hh:mm:ss.ms): "<<std::setfill('0')<<std::setw(2)<<h.count()<<":"<<std::setfill('0')<<std::setw(2)<<m.count()<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count();
    }else if(m.count()>0){
        time <<"build_time (mm:ss.ms): "<<std::setfill('0')<<std::setw(2)<<size_t(m.count())<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count();
    }else if(s.count()>0){
        time <<"build_time (ss.ms): "<<std::setfill('0')<<std::setw(2)<<size_t(s.count())<<"."<<ms.count();
    }else{
        time <<"build_time (ms): "<<ms.count();
    }

    return time.str();
}*/

#define build_dt(dt, suffix) \
{                            \
dt instance;                 \
auto t1 = std::chrono::high_resolution_clock::now();\
sdsl::construct(instance, plain_input_file, 1);\
auto t2 = std::chrono::high_resolution_clock::now();\
sdsl::store_to_file(instance, output_file+"."+suffix); \
std::cout<<suffix<<" "<<report_time(t1, t2, 0)<<",  space_usage:"<<float(sdsl::size_in_bytes(instance)*8)/float(instance.size())<<" bps"<<std::endl;\
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

/*void rl2plain(std::string& rl_file, std::string& output_plain_file){

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

template<class bwt_type>
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
}*/
#include <numeric>
template<typename func>
int measure_avg_time(func function, std::vector<std::pair<uint8_t, uint64_t>>& tests, const uint16_t* vec, std::string name) {

    std::vector<long long> durations;

    for (auto & test : tests) {
        auto start = std::chrono::high_resolution_clock::now();
        function(vec, test.first, test.second);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        durations.push_back(duration);
    }

    // Compute average duration
    long long totalDuration = std::accumulate(durations.begin(), durations.end(), 0LL);
    double averageDuration = static_cast<double>(totalDuration) / static_cast<double>(tests.size());
    std::cout << "Average running time "<<name<<": "<< averageDuration << " nanoseconds" << std::endl;
    return 0;
}

static inline uint64_t count_scalar(const uint16_t* stream, uint8_t sym, uint64_t idx){

    uint8_t alpha_bits = 4;
    uint8_t alpha_mask = 15;

    uint64_t r_len[2]={0};
    r_len[1] = stream[0]>>alpha_bits;
    size_t acc = r_len[1], i=0, rank=0;
    while(acc<idx){
        rank+= r_len[(stream[i] & alpha_mask)==sym];
        r_len[1] = stream[++i]>>alpha_bits;
        acc+= r_len[1];
    }
    return rank + (idx-(acc-r_len[1]))*((stream[i]&15)==sym);
}

#ifdef __ARM_NEON__
#include <arm_neon.h>
static inline uint64_t count_neon(const uint16_t * stream, uint8_t sym, uint64_t idx){

    const uint16x8_t sym_vec = vdupq_n_u16(sym);
    //const uint16x8_t zeroes = vdupq_n_u16(0);
    const uint16x8_t alpha_mask = vdupq_n_u16(15);
    //uint16x8_t rank_acc = vdupq_n_u16(0);
    size_t i=0;
    /*do{
        uint16x8_t block = vld1q_u16(&stream[i]);
        uint16x8_t bk_lengths = vshrq_n_u16(block, 4);
        uint16x8_t sym_mask = vceqq_u16(vandq_u16(block, alpha_mask), sym_vec);

        //prefix sum
        uint16x8_t pf_sum = vaddq_u16(bk_lengths, vextq_u16(zeroes, bk_lengths, 7));
        pf_sum = vaddq_u16(pf_sum, vextq_u16(zeroes, pf_sum, 6));
        pf_sum = vaddq_u16(pf_sum, vextq_u16(zeroes, pf_sum, 4));

        uint16x8_t counter = vcltq_u16(pf_sum, vdupq_n_u16(idx));
        sym_mask = vandq_u16(sym_mask, counter);

        rank_acc = vaddq_u16(rank_acc, vandq_u16(bk_lengths, sym_mask));
        uint64_t last_pfs = vgetq_lane_u16(pf_sum, 7);

        if(last_pfs>=idx){
            const uint8x8_t res = vshrn_n_u16(counter, 4);
            const uint64_t cnt_mask = vget_lane_u64(vreinterpret_u64_u8(res), 0);
            uint8_t last = 8 - (__builtin_clzll(cnt_mask)>>3);
            last_pfs = vmaxvq_u16(vandq_u16(pf_sum, counter));

            //horizontal add of rank_acc
            uint64x2_t tmp   = vpaddlq_u32(vpaddlq_u16(rank_acc));
            uint64x1_t rank = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp));
            return rank[0] + ((idx-last_pfs)*((stream[i+last] & 15)==sym));
        }
        idx-=last_pfs;
        i+=8;
    } while(true);*/
    uint16x8_t block = vld1q_u16(&stream[0]);
    uint16x8_t bk_lengths = vshrq_n_u16(block, 4);
    uint64x2_t tmp   = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0], rank=0;
    while(acc<idx){
        bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));
        tmp   = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        i+=8;
        block = vld1q_u16(&stream[i]);
        bk_lengths = vshrq_n_u16(block, 4);
        tmp   = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }
    return rank;
}
#endif

#ifdef __AVX2__
#include <immintrin.h>

#define _mm256_cmpge_epu16(a, b) \
        _mm256_cmpeq_epi16(_mm256_max_epu16(a, b), a)

#define _mm256_cmple_epu16(a, b) _mm256_cmpge_epu16(b, a)

static inline uint64_t count_avx2(const uint16_t * stream, uint8_t sym, uint64_t idx){


    /*const __m256i sym_vec = _mm256_set1_epi16(sym);
    const __m256i alpha_mask = _mm256_set1_epi16(15);
    size_t i=0;

    __m256i block = _mm256_loadu_si256((const __m256i*)&stream[0]);
    __m256i bk_lengths = _mm256_srli_epi16(block, 4);

    __m256i hadd = _mm256_hadd_epi16(bk_lengths, bk_lengths);
    __m128i sum128 = _mm_add_epi16(_mm256_extracti128_si256(hadd, 0), _mm256_extracti128_si256(hadd, 1));
    sum128 = _mm_hadd_epi16(sum128, sum128);
    sum128 = _mm_hadd_epi16(sum128, sum128);
    uint64_t acc = _mm_extract_epi16(sum128, 0);

    size_t rank=0;

    while(acc<idx){
        bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi16(_mm256_and_si256(block, alpha_mask), sym_vec));

        hadd = _mm256_hadd_epi16(bk_lengths, bk_lengths);
        sum128 = _mm_add_epi16(_mm256_extracti128_si256(hadd, 0), _mm256_extracti128_si256(hadd, 1));
        sum128 = _mm_hadd_epi16(sum128, sum128);
        sum128 = _mm_hadd_epi16(sum128, sum128);
        rank += _mm_extract_epi16(sum128, 0);

        i+=8;
        block = _mm256_load_si256((const __m256i*)stream);
        bk_lengths = _mm256_srli_epi16(block, 4);

        hadd = _mm256_hadd_epi16(bk_lengths, bk_lengths);
        sum128 = _mm_add_epi16(_mm256_extracti128_si256(hadd, 0), _mm256_extracti128_si256(hadd, 1));
        sum128 = _mm_hadd_epi16(sum128, sum128);
        sum128 = _mm_hadd_epi16(sum128, sum128);
        acc += _mm_extract_epi16(sum128, 0);
    }*/

    /*size_t tmp=0;
    for(size_t j=0;j<64;j++){
        tmp+=stream[j]>>4;
        std::cout<<j<<" -> sym:"<<(stream[j] & 15)<<" len:"<<(stream[j]>>4)<<" acc:"<<tmp<<std::endl;
    }*/

    const __m256i sym_vec = _mm256_set1_epi16(sym);
    const __m256i alpha_mask = _mm256_set1_epi16(15);
    __m256i rank_acc = _mm256_set1_epi16(0);
    size_t i=0;
    uint16_t acc[16]={0};

    do{
        __m256i block = _mm256_loadu_si256((const __m256i*)&stream[i]);
        __m256i bk_lengths = _mm256_srli_epi16(block, 4);
        __m256i sym_mask = _mm256_cmpeq_epi16(_mm256_and_si256(block, alpha_mask), sym_vec);

        _mm256_storeu_si256((__m256i *)&acc, bk_lengths);
        /*for(size_t j=0;j<16;j++){
            std::cout<<i<<" = "<<idx<<" -> "<<j<<" / "<<acc[j]<<std::endl;
        }*/

        __m256i pf_sum = _mm256_add_epi16(bk_lengths, _mm256_slli_si256(bk_lengths, 2));
        pf_sum = _mm256_add_epi16(pf_sum, _mm256_slli_si256(pf_sum, 4));
        pf_sum = _mm256_add_epi16(pf_sum, _mm256_slli_si256(pf_sum, 8));
        _mm256_storeu_si256((__m256i *)&acc, pf_sum);
        uint16_t last_pf = acc[7]+acc[15];
        pf_sum = _mm256_add_epi16(pf_sum, _mm256_set_m128i( _mm_set1_epi16(acc[7]), _mm_set1_epi16(0)));

        __m256i counter =_mm256_cmpgt_epi16(_mm256_set1_epi16(idx),pf_sum);
        sym_mask = _mm256_and_si256(sym_mask, counter);
        rank_acc = _mm256_add_epi16(rank_acc, _mm256_and_si256(bk_lengths, sym_mask));

        //TODO
        /*_mm256_storeu_si256((__m256i *)&acc, pf_sum);
        for(size_t j=0;j<16;j++){
            std::cout<<idx<<" -> "<<j<<" / "<<acc[j]<<std::endl;
        }
        _mm256_storeu_si256((__m256i *)&acc, counter);
        for(size_t j=0;j<16;j++){
            std::cout<<j<<" | "<<acc[j]<<std::endl;
        }*/
        //

        if(idx<=last_pf){
            uint32_t cnt_mask = _mm256_movemask_epi8(counter);
            uint64_t last = (_mm_popcnt_u32(cnt_mask)>>1);

            __m256i hadd = _mm256_hadd_epi16(rank_acc, rank_acc);
            __m128i sum128 = _mm_add_epi16(_mm256_extracti128_si256(hadd, 0), _mm256_extracti128_si256(hadd, 1));
            sum128 = _mm_hadd_epi16(sum128, sum128);
            sum128 = _mm_hadd_epi16(sum128, sum128);
            uint64_t rank = _mm_extract_epi16(sum128, 0);
            return rank + ((idx-last_pf)*((stream[i+last] & 15)==sym));
        }
        i+=16;
        idx -= last_pf;
    } while(true);
}
#endif

#include <algorithm>
#include <random>
#include <set>
// Function to generate a strictly increasing sequence of random numbers within a specified range
std::vector<uint16_t> generate_strictly_increasing_sequence(int count, int min_value, int max_value) {
    if (count > (max_value - min_value + 1)) {
        throw std::invalid_argument("Count is too large for the specified range.");
    }

    std::vector<uint16_t> sequence;
    //std::random_device rd;
    //std::mt19937 gen(rd());

    std::uniform_int_distribution<> dis(min_value, max_value);
    std::set<uint16_t> idxs;

    while (idxs.size() < count) {
        uint16_t num = rand() % max_value;
        idxs.insert(num);
    }

    for(auto &pos : idxs){
        sequence.push_back(pos);
    }
    return sequence;
}

void test_scan(){
    /*uint16_t vect[64] = {
            0x1A2B, 0x3C48, 0x5E6F, 0x7F80, 0x9A9B, 0xBCBC, 0xDED8, 0xF0F8,
            0x1122, 0x3344, 0x5566, 0x7788, 0x99AA, 0xBBCC, 0xDDEE, 0xFFE0,
            0x1234, 0x5678, 0x9ABC, 0xDEF0, 0x1357, 0x9BDF, 0x2468, 0xACEF,
            0x1987, 0x7654, 0x3210, 0xFEDC, 0xBA98, 0x7654, 0x3210, 0x1234,
            0x5678, 0x9ABC, 0xDEF1, 0x1357, 0x9BDF, 0x2468, 0xACEF, 0x1987,
            0x7654, 0x3210, 0xFEDC, 0xBA98, 0x7654, 0x3210, 0x1234, 0x5678,
            0x9ABC, 0xDEF1, 0x1357, 0x9BDF, 0x2468, 0xACEF, 0x1987, 0x7654,
            0x3210, 0xFEDC, 0xBA98, 0x7654, 0x3210, 0x1234, 0x5678, 0x9ABC
    };*/

    std::vector<uint16_t> vec(64);
    std::vector<uint16_t> sequence = generate_strictly_increasing_sequence(65, 0, 32767);
    size_t acc=0;
    for(size_t i=65;i-->1;){
        sequence[i]=sequence[i]-sequence[i-1];
        acc+=sequence[i];
        uint16_t sym = rand()%16;
        sequence[i] = sequence[i]<<4 | sym;
    }
    std::cout<<acc<<std::endl;
    sequence.pop_back();
    assert(acc<32767);

    std::vector<std::pair<uint8_t, uint64_t>> tests(20000);
    for(auto & test : tests){
        test.first = rand() % 16;
        test.second = rand() % acc;
        //std::cout<<tests[i].second<<std::endl;
    }

    for(auto & test : tests){
        //std::cout<<int(test.first)<<" "<<test.second<<"->  avx2:"<<count_avx2(sequence.data(), test.first, test.second)<<" scalar:"<<count_scalar(sequence.data(), test.first, test.second)<<std::endl;
        //assert(count_neon(vect, test.first, test.second)==count_scalar(vect, test.first, test.second));
    }
    //std::cout<<"whut? "<<small<<" "<<tests.size()<<std::endl;

#ifdef __ARM_NEON__
    measure_avg_time(count_neon, tests, sequence.data(), "neon");
#endif

#ifdef __AVX2__
    measure_avg_time(count_avx2, tests, sequence.data(), "avx2");
#endif

    measure_avg_time(count_scalar, tests, sequence.data(), "scalar");

    /*{
        std::string name = "count neon";
        pretty_print(tests.size(), tests.size(), name,
                     bench([&tests, &vect]() {
                         for(auto & test : tests){
                             count_neon(vect, test.first, test.second);
                         }
                     }));
    }

    {
        std::string name = "count scalar";
        pretty_print(tests.size(), tests.size(), name,
                     bench([&tests, &vect]() {
                         for(auto & test : tests){
                             count_neon(vect, test.first, test.second);
                         }
                     }));
    }*/
}

int main(int argc, char** argv){

    if(argc!=4){
        std::cout<<"usage: ./build_bwt_dts plain_bwt.rl_bwt alphabet output_file"<<std::endl;
        exit(1);
    }

    std::string input_file = std::string(argv[1]);
    char *pend;
    long int alphabet = strtol(argv[2], &pend, 10);
    //assert(alphabet>2 && alphabet<=16);
    std::string output_file = std::string(argv[3]);

    using bwt_type = rlbwt_vlb<4096, 64, 4>;
    bwt_type bwt;
    build_rlbwt_vlb<bwt_type>(bwt, input_file, INPUT_FORMAT::GRL_BWT);

    //build_my_bwt<rlbwt_small_alpha<16>>(input_file, output_file);
    /*std::string plain_input_file = "tmp_plain.txt";
    rl2plain(input_file, plain_input_file);
    std::cout<<"Creating wavelet trees for "<<input_file<<std::endl;
    TESTED_DTS

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
    }*/
}