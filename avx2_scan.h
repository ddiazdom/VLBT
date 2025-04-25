//
// Created by Diaz, Diego on 25.4.2025.
//

#ifndef BWT_DTS_BENCHMARKS_AVX_SCAN_H
#define BWT_DTS_BENCHMARKS_AVX_SCAN_H

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
#endif //BWT_DTS_BENCHMARKS_AVX_SCAN_H
