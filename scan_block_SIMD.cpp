//
// Created by Diaz, Diego on 2.4.2025.
//
#include <iostream>
#include <cassert>
#include "benchmarks/perf_utils.h"

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

int main(){
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