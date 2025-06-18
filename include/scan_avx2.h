//
// Created by Diaz, Diego on 25.4.2025.
//

#ifndef VLBT_SCAN_AVX2_H
#define VLBT_SCAN_AVX2_H

#include <immintrin.h>
#include "simd_tables.h"
#include "utils.h"

#define _mm256_cmpge_epu8(a, b) \
        _mm256_cmpeq_epi8(_mm256_max_epu8(a, b), a)
#define _mm256_cmple_epu8(a, b) _mm256_cmpge_epu8(b, a)

#define _mm256_cmpge_epu16(a, b) \
        _mm256_cmpeq_epi16(_mm256_max_epu16(a, b), a)
#define _mm256_cmple_epu16(a, b) _mm256_cmpge_epu16(b, a)

#define _mm256_cmpge_epu32(a, b) \
        _mm256_cmpeq_epi32(_mm256_max_epu32(a, b), a)
#define _mm256_cmple_epu32(a, b) _mm256_cmpge_epu32(b, a)

#define _mm256_cmpge_epu64(a, b) \
        _mm256_cmpeq_epi64(_mm256_max_epu64(a, b), a)
#define _mm256_cmple_epu64(a, b) _mm256_cmpge_epu64(b, a)

static inline void print8x32(__m256i vec) {
    uint8_t values[32];
    _mm256_storeu_si256((__m256i*)values, vec); // Unaligned store
    for (int i = 0; i < 32; ++i) {
        std::cout<<int(values[i])<<" ";
    }
    std::cout<<""<<std::endl;
}

static inline void print16x16(__m256i vec) {
    uint16_t values[16];
    _mm256_storeu_si256((__m256i*)values, vec); // Unaligned store
    for (int i = 0; i < 16; ++i) {
        std::cout<<int(values[i])<<" ";
    }
    std::cout<<""<<std::endl;
}

static inline void print32x8(__m256i vec) {
    uint32_t values[8];
    _mm256_storeu_si256((__m256i*)values, vec); // Unaligned store
    for (int i = 0; i < 8; ++i) {
        std::cout<<values[i]<<" ";
    }
    std::cout<<""<<std::endl;
}

static inline void print64x4(__m256i vec) {
    uint64_t values[4];
    _mm256_storeu_si256((__m256i*)values, vec); // Unaligned store
    for (int i = 0; i < 4; ++i) {
        std::cout<<values[i]<<" ";
    }
    std::cout<<""<<std::endl;
}

static inline __m256i _m256_shift_right_epi8(__m256i input, uint8_t shift) {
    //from https://wunkolo.github.io/post/2020/11/gf2p8affineqb-int8-shifting/
    switch (shift) {
        case 0: return input;
        case 1: return _mm256_and_si256(_mm256_srli_epi64(input, 1), _mm256_set1_epi8(127));
        case 2: return _mm256_and_si256(_mm256_srli_epi64(input, 2), _mm256_set1_epi8(63));
        case 3: return _mm256_and_si256(_mm256_srli_epi64(input, 3), _mm256_set1_epi8(31));
        case 4: return _mm256_and_si256(_mm256_srli_epi64(input, 4), _mm256_set1_epi8(15));
        case 5: return _mm256_and_si256(_mm256_srli_epi64(input, 5), _mm256_set1_epi8(7));
        case 6: return _mm256_and_si256(_mm256_srli_epi64(input, 6), _mm256_set1_epi8(3));
        case 7: return _mm256_and_si256(_mm256_srli_epi64(input, 7), _mm256_set1_epi8(1));
        default: return _mm256_setzero_si256(); // for shift >= 8
    }
}

static inline __m256i _m256_shift_right_epi16(__m256i input, int shift) {
    switch (shift) {
        case 0: return input;
        case 1: return _mm256_srli_epi16(input, 1);
        case 2: return _mm256_srli_epi16(input, 2);
        case 3: return _mm256_srli_epi16(input, 3);
        case 4: return _mm256_srli_epi16(input, 4);
        case 5: return _mm256_srli_epi16(input, 5);
        case 6: return _mm256_srli_epi16(input, 6);
        case 7: return _mm256_srli_epi16(input, 7);
        default: return _mm256_setzero_si256(); // for shift >= 8
    }
}

static inline __m256i _m256_shift_right_epi32(__m256i input, int shift) {
    switch (shift) {
        case 0: return input;
        case 1: return _mm256_srli_epi32(input, 1);
        case 2: return _mm256_srli_epi32(input, 2);
        case 3: return _mm256_srli_epi32(input, 3);
        case 4: return _mm256_srli_epi32(input, 4);
        case 5: return _mm256_srli_epi32(input, 5);
        case 6: return _mm256_srli_epi32(input, 6);
        case 7: return _mm256_srli_epi32(input, 7);
        default: return _mm256_setzero_si256(); // for shift >= 8
    }
}

static inline __m256i m256_shift_right_epi64(__m256i input, int shift) {
    switch (shift) {
        case 0: return input;
        case 1: return _mm256_srli_epi64(input, 1);
        case 2: return _mm256_srli_epi64(input, 2);
        case 3: return _mm256_srli_epi64(input, 3);
        case 4: return _mm256_srli_epi64(input, 4);
        case 5: return _mm256_srli_epi64(input, 5);
        case 6: return _mm256_srli_epi64(input, 6);
        case 7: return _mm256_srli_epi64(input, 7);
        default: return _mm256_setzero_si256(); // for shift >= 8
    }
}

static inline void m256_psum_epi8_ovf(__m256i input, uint32_t idx, uint32_t& pf_sum, uint32_t& idx_run) {

    __m256i halves[2];

    uint16_t tmp[16];
    const __m256i idx_vec = _mm256_set1_epi16(idx);
    __m128i low_bytes  = _mm256_extracti128_si256(input, 0);
    __m128i high_bytes = _mm256_extracti128_si256(input, 1);

    halves[0] = _mm256_cvtepu8_epi16(low_bytes);
    halves[1] = _mm256_cvtepu8_epi16(high_bytes);

    halves[0] = _mm256_add_epi16(halves[0], _mm256_slli_si256(halves[0], 2));
    halves[0] = _mm256_add_epi16(halves[0], _mm256_slli_si256(halves[0], 4));
    halves[0] = _mm256_add_epi16(halves[0], _mm256_slli_si256(halves[0], 8));
    _mm256_storeu_si256((__m256i *)&tmp, halves[0]);
    halves[0] = _mm256_add_epi16(halves[0], _mm256_set_m128i(_mm_set1_epi16(tmp[7]), _mm_set1_epi16(0)));

    //print16x16(halves[0]);

    __m256i idx_mask = _mm256_cmple_epu16(halves[0], idx_vec);//mask for <=idx
    uint64_t lt_low = (uint32_t) _mm256_movemask_epi8(idx_mask);
    //print16x8(idx_mask);

    halves[1] = _mm256_add_epi16(halves[1], _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(tmp[7]+tmp[15])));
    halves[1] = _mm256_add_epi16(halves[1], _mm256_slli_si256(halves[1], 2));
    halves[1] = _mm256_add_epi16(halves[1], _mm256_slli_si256(halves[1], 4));
    halves[1] = _mm256_add_epi16(halves[1], _mm256_slli_si256(halves[1], 8));
    _mm256_storeu_si256((__m256i *)&tmp, halves[1]);
    halves[1] = _mm256_add_epi16(halves[1], _mm256_set_m128i(_mm_set1_epi16(tmp[7]), _mm_set1_epi16(0)));

    //print16x16(halves[1]);

    idx_mask = _mm256_cmple_epu16(halves[1], idx_vec);//mask for <=idx
    //print16x16(idx_mask);

    uint64_t lt_high = (uint32_t) _mm256_movemask_epi8(idx_mask);
    uint64_t lt = lt_low | (lt_high<<16);//two bits per element
    idx_run = __builtin_ctzll(~lt)>>1;

    _mm256_storeu_si256((__m256i*)tmp, halves[idx_run>>3]);
    pf_sum = tmp[idx_run & 7];
}

static inline uint32_t m256_hsum_epi8_ovf(__m256i input) {

    __m128i low_bytes  = _mm256_extracti128_si256(input, 0);
    __m128i high_bytes = _mm256_extracti128_si256(input, 1);
    __m256i low = _mm256_cvtepu8_epi16(low_bytes);
    __m256i high = _mm256_cvtepu8_epi16(high_bytes);

    const __m256i sum1 = _mm256_add_epi16(low, high);
    //print16x16(sum1);
    const __m256i sum2 = _mm256_add_epi16(sum1, _mm256_srli_si256(sum1, 2));
    //print16x16(sum2);
    const __m256i sum3 = _mm256_add_epi16(sum2, _mm256_srli_si256(sum2, 4));
    //print16x16(sum3);
    const __m256i sum4 = _mm256_add_epi16(sum3, _mm256_srli_si256(sum3, 8));
    //print16x16(sum4);

    uint16_t tmp[16];
    _mm256_storeu_si256((__m256i*)tmp, sum4);
    return tmp[0]+tmp[8];
}

static inline uint32_t m256_hsum_epi8(__m256i input) {
    __m256i sad = _mm256_sad_epu8(input, _mm256_setzero_si256());
    sad = _mm256_add_epi64(sad, _mm256_srli_si256(sad, 8));

    __m256i v_shifted = _mm256_permute4x64_epi64(sad, _MM_SHUFFLE(0, 0, 0, 2));
    sad = _mm256_add_epi64(sad, v_shifted);

    return _mm256_cvtsi256_si32(sad);
}

static inline void m256_psum_epi16_ovf(__m256i input, uint32_t idx, uint32_t& pf_sum, uint32_t& idx_run) {

    __m256i halves[2];

    uint32_t tmp[8];
    const __m256i idx_vec =  _mm256_set1_epi32(idx);

    __m128i low_bytes  = _mm256_extracti128_si256(input, 0);
    __m128i high_bytes = _mm256_extracti128_si256(input, 1);

    halves[0] = _mm256_cvtepu16_epi32(low_bytes);
    halves[1] = _mm256_cvtepu16_epi32(high_bytes);

    halves[0] = _mm256_add_epi32(halves[0], _mm256_slli_si256(halves[0], 4));
    halves[0] = _mm256_add_epi32(halves[0], _mm256_slli_si256(halves[0], 8));
    _mm256_storeu_si256((__m256i *)&tmp, halves[0]);
    halves[0] = _mm256_add_epi32(halves[0], _mm256_set_m128i(_mm_set1_epi16(tmp[3]), _mm_set1_epi32(0)));

    //print32x8(halves[0]);

    __m256i idx_mask = _mm256_cmple_epu32(halves[0], idx_vec);//mask for >idx
    uint64_t lt_low = (uint32_t)_mm256_movemask_epi8(idx_mask);//4 bits represent one element

    halves[1] = _mm256_add_epi32(halves[1], _mm256_set_epi32(0,0,0,0,0,0,0,(tmp[3]+tmp[7])));
    halves[1] = _mm256_add_epi32(halves[1], _mm256_slli_si256(halves[1], 4));
    halves[1] = _mm256_add_epi32(halves[1], _mm256_slli_si256(halves[1], 8));
    _mm256_storeu_si256((__m256i *)&tmp, halves[1]);
    halves[1] = _mm256_add_epi16(halves[1], _mm256_set_m128i(_mm_set1_epi32(tmp[3]), _mm_set1_epi32(0)));

    //print32x8(halves[1]);
    //print32x8(idx_vec);

    idx_mask = _mm256_cmple_epu32(halves[1], idx_vec);//4 bits represent one element
    uint64_t lt_high = (uint32_t)_mm256_movemask_epi8(idx_mask);//4 bits represent one element

    lt_low |= lt_high<<16;
    //print32x8(idx_mask);
    idx_run = __builtin_ctzll(~lt_low)>>2;//

    uint32_t pf_sum_vec[8];
    _mm256_storeu_si256((__m256i*)&pf_sum_vec, halves[idx_run>>2]);

    pf_sum = pf_sum_vec[idx_run & 3];
}

static inline uint32_t m256_hsum_epi16_ovf(__m256i input) {

    const __m256i sum = _mm256_add_epi32(_mm256_shuffle_epi8(input, _mm256_set_epi8(-1,-1,7,6, -1,-1,5,4, -1,-1,3,2, -1,-1,1,0,
                                                                                    -1,-1,7,6, -1,-1,5,4, -1,-1,3,2, -1,-1,1,0)),
                                         _mm256_shuffle_epi8(input, _mm256_set_epi8(-1,-1,15,14, -1,-1,13,12, -1,-1,11,10, -1,-1,9,8,
                                                                                    -1,-1,15,14, -1,-1,13,12, -1,-1,11,10, -1,-1,9,8)));

    const __m256i sum2 = _mm256_hadd_epi32(sum, sum);

    __m128i sum_high = _mm256_extractf128_si256(sum2, 1);
    __m128i result = _mm_add_epi64(sum_high, _mm256_castsi256_si128(sum));
    return _mm_cvtsi128_si32(result);
}

static inline uint32_t hsum_epi16(__m128i input) {
    const __m128i sum1 = _mm_add_epi16(input, _mm_srli_si128(input, 2));
    const __m128i sum2 = _mm_add_epi16(sum1, _mm_srli_si128(sum1, 4));
    const __m128i sum3 = _mm_add_epi16(sum2, _mm_srli_si128(sum2, 8));
    return _mm_cvtsi128_si32(sum3);
}

static inline uint32_t hsum_epi32(__m128i input) {
    const __m128i sum1 = _mm_add_epi32(input, _mm_srli_si128(input, 4));
    const __m128i sum2 = _mm_add_epi32(sum1, _mm_srli_si128(sum1, 8));
    return _mm_cvtsi128_si32(sum2);
}

static inline uint64_t hsum_sum_epi64(__m128i input) {
    const __m128i sum = _mm_add_epi32(input, _mm_srli_si128(input, 8));
    return _mm_cvtsi128_si64(sum);
}

template<bool vbyte_compressed, uint8_t ctr_width, uint8_t bytes_per_run>
static inline __m128i decode_block_sse42(const uint8_t **stream) {

    if constexpr (vbyte_compressed) {
        uint8_t *pshuf, ctrl_bits = **stream;
        uint8_t len;
        if constexpr (ctr_width == 1) {
            pshuf = (uint8_t *) &dec_table_16x8[ctrl_bits];
            len = pshuf[14 + (ctrl_bits >> 7)] + 1;
        } else if constexpr (ctr_width == 2) {
            pshuf = (uint8_t *) &dec_table_32x4[ctrl_bits];
            len = pshuf[12 + (ctrl_bits >> 6)] + 1;
        } else if constexpr (ctr_width == 3) {
            pshuf = (uint8_t *) &dec_table_64x2[ctrl_bits];
            len = pshuf[8 + ((ctrl_bits >> 3) & 7)] + 1;
        } else {
            exit(1);
        }
        __m128i compressed =  _mm_loadu_si128((const __m128i*)(*stream+1));
        __m128i dec_shuffle = _mm_loadu_si128((const __m128i*)pshuf);
        __m128i data = _mm_shuffle_epi8(compressed, dec_shuffle);
        *stream += len + 1;

        return data;

    } else {
        __m128i data = _mm_loadu_si128((const __m128i*)(*stream));

        //align the bytes
        if constexpr (bytes_per_run==3){
            const __m128i dec_shuff = _mm_set_epi8(-1,11,10,9, -1,8,7,6, -1,5,4,3, -1,2,1,0);
            data = _mm_shuffle_epi8(data, dec_shuff);
            *stream+=12;
        } else if constexpr (bytes_per_run==5){
            const __m128i dec_shuff = _mm_set_epi8(-1,-1,-1,9,8,7,6,5, -1,-1,-1,4,3,2,1,0);
            data = _mm_shuffle_epi8(data, dec_shuff);
            *stream+=10;
        } else if constexpr (bytes_per_run==6){
            const __m128i dec_shuff = _mm_set_epi8(-1,-1,11,10,9,8,7,6, -1,-1,5,4,3,2,1,0);
            data = _mm_shuffle_epi8(data, dec_shuff);
            *stream+=12;
        } else if constexpr (bytes_per_run==7){
            const __m128i dec_shuff = _mm_set_epi8(-1,13,12,11,10,9,8,7, -1,6,5,4,3,2,1,0);
            data = _mm_shuffle_epi8(data, dec_shuff);
            *stream+=14;
        } else {
            *stream+=16;
        }
        return data;
    }
}

template<bool overflow16, bool overflow32=false>
static inline std::pair<uint64_t, uint8_t> inv_select_sse42_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx) {

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    __m128i block = _mm_loadu_si128((const __m128i*)*stream);
    *stream+=16;
    __m128i bk_lengths =  shift_right_epi8(block, sigma_bits);
    uint32_t prev_acc = 0, acc;

    //print8x16(bk_lengths);

    //sometimes the back of the block has garbage, so I have to assume overflow
    acc = hsum_epi8_ovf(bk_lengths);

    while(acc<=idx){
        block = _mm_loadu_si128((const __m128i*)*stream);
        *stream+=16;
        bk_lengths =  shift_right_epi8(block, sigma_bits);

        //print8x16(bk_lengths);
        prev_acc = acc;
        acc += hsum_epi8_ovf(bk_lengths);
        l++;
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    uint8_t alpha_m = (1UL << sigma_bits)-1, run, sym;
    const __m128i alpha_mask = _mm_set1_epi8(alpha_m);
    __m128i sym_vec;

    if constexpr(overflow16){
        psum_epi8_ovf(bk_lengths, idx, pf_sum, idx_run);
        run = stream_start[l*16+idx_run];
        sym = run & alpha_m;
        sym_vec = _mm_set1_epi8(sym);
    }else{
        //prefix sum without overflow
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 1));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 2));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 4));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 8));
        //print8x16(bk_lengths);
        const __m128i idx_mask = _mm_cmple_epu8(bk_lengths, _mm_set1_epi8(idx));//mask for >idx
        uint16_t less_than = _mm_movemask_epi8(idx_mask);
        idx_run = __builtin_ctzll(~less_than);

        //print8x16(idx_mask);
        const __m128i shuff = _mm_set1_epi8(idx_run);
        const __m128i run_vec = _mm_shuffle_epi8(block, shuff);
        //print8x16(shuff);
        run = ((uint8_t*)&run_vec)[0];

        sym = run & alpha_m;
        sym_vec = _mm_and_si128(run_vec, alpha_mask);

        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
        pf_sum = ((uint8_t*)&pf_sum_vec)[0];
    }

    pf_sum-= run>>sigma_bits;

    *stream = stream_start;
    //compute the rank answer
    uint64_t rank = 0;
    for(size_t p=0;p<l;p++){
        block = _mm_loadu_si128((const __m128i*)*stream);
        *stream+=16;
        bk_lengths =  shift_right_epi8(block, sigma_bits);
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank += hsum_epi8_ovf(bk_lengths);
        } else {
            rank += hsum_epi8(bk_lengths);
        }
    }

    block = _mm_loadu_si128((const __m128i*)*stream);
    bk_lengths =  shift_right_epi8(block, sigma_bits);

    bk_lengths = _mm_and_si128(bk_lengths, _mm_loadu_si128((const __m128i*)mask8x16[idx_run]));
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));

    if constexpr (overflow16){
        rank += hsum_epi8_ovf(bk_lengths);
    }else{
        rank += hsum_epi8(bk_lengths);
    }

    rank +=idx-pf_sum;
    return {rank, sym};
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline std::pair<uint64_t, uint8_t> inv_select_sse42_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx) {

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    __m128i block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    __m128i bk_lengths =  shift_right_epi16(block, sigma_bits);
    uint32_t prev_acc = 0, acc;

    //some time the block has some garbage, so we have to assume overflow at the end
    acc = hsum_epi16_ovf(bk_lengths);
    //print16x8(bk_lengths);

    while(acc<=idx){
        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        bk_lengths =  shift_right_epi16(block, sigma_bits);
        //print16x8(bk_lengths);

        prev_acc = acc;
        acc += hsum_epi16_ovf(bk_lengths);
        l++;
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi16(alpha_m);

    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx, pf_sum, idx_run);
    }else{
        //prefix sum without overflow
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 2));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 4));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 8));
        //print16x8(bk_lengths);

        const __m128i idx_mask = _mm_cmple_epu16(bk_lengths, _mm_set1_epi16(idx));//mask for <=idx
        uint8_t less_than = (uint8_t)_mm_movemask_epi8(_mm_packs_epi16(idx_mask, _mm_setzero_si128()));
        idx_run = __builtin_ctzll(~less_than);
        //print16x8(idx_mask);
    }

    const __m128i shuff = _mm_add_epi8(_mm_set_epi8(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0),
                                       _mm_set1_epi8(idx_run<<1));
    //print8x16(shuff);

    const __m128i run_vec = _mm_shuffle_epi8(block, shuff);
    uint16_t run = (uint16_t)_mm_cvtsi128_si32(run_vec);
    uint8_t sym = run & alpha_m;
    const __m128i sym_vec = _mm_and_si128(run_vec, alpha_mask);

    if constexpr (!overflow8){
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
        pf_sum = (uint16_t)_mm_cvtsi128_si32(pf_sum_vec);
    }

    pf_sum-= run>>sigma_bits;

    *stream = stream_start;
    //compute the rank answer
    uint64_t rank = 0;
    for(size_t p=0;p<l;p++){
        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        bk_lengths = shift_right_epi16(block, sigma_bits);
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank += hsum_epi16_ovf(bk_lengths);
        }else{
            rank += hsum_epi16(bk_lengths);
        }
    }

    block = decode_block_sse42<vbyte_compressed,1,2>(stream);
    bk_lengths =  shift_right_epi16(block, sigma_bits);

    bk_lengths = _mm_and_si128(bk_lengths, _mm_loadu_si128((const __m128i*)mask16x8[idx_run]));
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));

    if constexpr (overflow8){
        rank += hsum_epi16_ovf(bk_lengths);
    }else{
        rank += hsum_epi16(bk_lengths);
    }

    rank +=idx-pf_sum;
    return {rank, sym};
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint8_t> inv_select_sse42_32x4(const uint8_t ** stream, uint8_t sigma, uint64_t idx) {

    //TODO: assert idx fits 2 bytes
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    __m128i block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    __m128i bk_lengths =  shift_right_epi32(block, sigma_bits);
    uint32_t prev_acc=0, acc;

    acc= hsum_epi32(bk_lengths);
    //print32x4(bk_lengths);

    while(acc<=idx){
        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        bk_lengths =  shift_right_epi32(block, sigma_bits);

        //print32x4(bk_lengths);

        prev_acc = acc;
        acc += hsum_epi32(bk_lengths);
        l++;
    }

    idx-=prev_acc;

    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 4));
    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 8));

    //print32x4(bk_lengths);

    __m128i idx_mask = _mm_cmple_epu32(bk_lengths, _mm_set1_epi32(idx));//mask for >idx

    //print32x4(idx_mask);

    uint64_t less_than = _mm_cvtsi128_si64(_mm_packs_epi32(idx_mask, _mm_setzero_si128()));
    uint8_t idx_run = __builtin_ctzll(~less_than)>>4;

    uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi32(alpha_m);

    const __m128i shuff_idxs = _mm_set_epi8(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
    const __m128i shuff = _mm_add_epi8(shuff_idxs, _mm_set1_epi8(idx_run<<2));

    //print8x16(shuff);

    const __m128i run_vec = _mm_shuffle_epi8(block, shuff);
    uint32_t run = _mm_cvtsi128_si32(run_vec);

    uint8_t sym = run & alpha_m;
    const __m128i sym_vec = _mm_and_si128(run_vec, alpha_mask);

    const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
    uint32_t pf_sum = _mm_cvtsi128_si32(pf_sum_vec);
    pf_sum-= run>>sigma_bits;

    *stream = stream_start;

    //compute the rank answer
    uint64_t rank = 0;
    for(size_t p=0;p<l;p++){
        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        bk_lengths = shift_right_epi32(block, sigma_bits);
        //print32x4(_mm_and_si128(block, alpha_mask));

        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));
        //print32x4(bk_lengths);
        rank += hsum_epi32(bk_lengths);
    }

    block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    bk_lengths =  shift_right_epi32(block, sigma_bits);

    bk_lengths = _mm_and_si128(bk_lengths, _mm_loadu_si128((const __m128i*)mask32x4[idx_run]));
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));
    //print32x4(bk_lengths);

    rank += hsum_epi32(bk_lengths);
    rank +=idx-pf_sum;

    return {rank, sym};
}

//byte compressed by default
template<uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint8_t> inv_select_sse42_64x2(const uint8_t ** stream, uint8_t sigma, uint64_t idx){
    return {0,0};
}

template<bool overflow16, bool overflow32=false>
static inline uint8_t access_sse42_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx){
    return 0;
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline uint8_t access_sse42_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx){
    return 0;
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline uint8_t access_sse42_32x4(const uint8_t **stream, uint8_t sigma, uint64_t idx){
    return 0;
}

template<uint8_t bytes_per_run>
static inline uint8_t access_sse42_64x2(const uint8_t **stream, uint8_t sigma, uint64_t idx){
    return 0;
}
#endif //VLBT_SCAN_AVX2_H


/*static inline uint64_t count_avx2(const uint16_t * stream, uint8_t sym, uint64_t idx){
    / *const __m256i sym_vec = _mm256_set1_epi16(sym);
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
    }* /

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
        }* /

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
        /* _mm256_storeu_si256((__m256i *)&acc, pf_sum);
        for(size_t j=0;j<16;j++){
            std::cout<<idx<<" -> "<<j<<" / "<<acc[j]<<std::endl;
        }
        _mm256_storeu_si256((__m256i *)&acc, counter);
        for(size_t j=0;j<16;j++){
            std::cout<<j<<" | "<<acc[j]<<std::endl;
        }* /
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
}*/
