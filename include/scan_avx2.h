//
// Created by Diaz, Diego on 25.4.2025.
//

#ifndef VLBT_SCAN_AVX2_H
#define VLBT_SCAN_AVX2_H

#include <immintrin.h>
#include "simd_tables.h"
#include "utils.h"

#define _mm256_set_m128i(v0, v1) \
    _mm256_insertf128_si256(_mm256_castsi128_si256(v1), (v0), 1)

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

static inline __m256i _m256_shift_right_epi8(__m256i input, const uint8_t shift) {
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

static inline __m256i _m256_shift_right_epi16(__m256i input, const uint8_t shift) {
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

static inline __m256i _m256_shift_right_epi32(__m256i input, const uint8_t shift) {
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

static inline __m256i m256_shift_right_epi64(__m256i input, const uint8_t shift) {
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

static inline void _m256_psum_epi8_ovf(const __m256i& input, const uint32_t idx, uint32_t& pf_sum, uint32_t& idx_run) {

    const __m256i idx_vec = _mm256_set1_epi16(idx);
    const __m128i low_bytes  = _mm256_extracti128_si256(input, 0);
    const __m128i high_bytes = _mm256_extracti128_si256(input, 1);

    __m256i halves[2];
    halves[0] = _mm256_cvtepu8_epi16(low_bytes);
    halves[1] = _mm256_cvtepu8_epi16(high_bytes);

    uint16_t tmp[16];
    halves[0] = _mm256_add_epi16(halves[0], _mm256_slli_si256(halves[0], 2));
    halves[0] = _mm256_add_epi16(halves[0], _mm256_slli_si256(halves[0], 4));
    halves[0] = _mm256_add_epi16(halves[0], _mm256_slli_si256(halves[0], 8));
    _mm256_storeu_si256((__m256i *)&tmp, halves[0]);
    halves[0] = _mm256_add_epi16(halves[0], _mm256_set_m128i(_mm_set1_epi16(tmp[7]), _mm_set1_epi16(0)));


    __m256i idx_mask = _mm256_cmple_epu16(halves[0], idx_vec);//mask for <=idx
    const uint64_t lt_low = static_cast<uint32_t>(_mm256_movemask_epi8(idx_mask));//the cast is necessary to avoid incorrect casting
    //print16x8(idx_mask);

    halves[1] = _mm256_add_epi16(halves[1], _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,tmp[7]+tmp[15]));

    halves[1] = _mm256_add_epi16(halves[1], _mm256_slli_si256(halves[1], 2));
    halves[1] = _mm256_add_epi16(halves[1], _mm256_slli_si256(halves[1], 4));
    halves[1] = _mm256_add_epi16(halves[1], _mm256_slli_si256(halves[1], 8));
    _mm256_storeu_si256((__m256i *)&tmp, halves[1]);
    halves[1] = _mm256_add_epi16(halves[1], _mm256_set_m128i(_mm_set1_epi16(tmp[7]), _mm_set1_epi16(0)));

    //print16x16(halves[0]);
    //print16x16(halves[1]);

    idx_mask = _mm256_cmple_epu16(halves[1], idx_vec);//mask for <=idx
    //print16x16(idx_mask);

    const uint64_t lt_high = static_cast<uint32_t>(_mm256_movemask_epi8(idx_mask));//the cast is necessary to avoid incorrect casting
    const uint64_t lt = lt_low | (lt_high<<32);//two bits per element
    idx_run = __builtin_ctzll(~lt)>>1;

    _mm256_storeu_si256((__m256i*)tmp, halves[idx_run>15]);
    pf_sum = tmp[idx_run & 15];
}

static inline uint32_t _m256_hsum_epi8_ovf(__m256i input) {

    const __m128i low_bytes  = _mm256_extracti128_si256(input, 0);
    const __m128i high_bytes = _mm256_extracti128_si256(input, 1);
    const __m256i low = _mm256_cvtepu8_epi16(low_bytes);
    const __m256i high = _mm256_cvtepu8_epi16(high_bytes);

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

static inline uint32_t _m256_hsum_epi8(const __m256i& input) {
    //__m256i sad = _mm256_sad_epu8(input, _mm256_setzero_si256());
    //sad = _mm256_add_epi64(sad, _mm256_srli_si256(sad, 8));
    //const auto v_shifted = _mm256_permute4x64_epi64(sad, _MM_SHUFFLE(0, 0, 0, 2));
    //sad = _mm256_add_epi64(sad, v_shifted);
    //return _mm_cvtsi128_si32(_mm256_castsi256_si128(sad));
    const __m256i sum1 = _mm256_add_epi8(input, _mm256_srli_si256(input, 1));
    const __m256i sum2 = _mm256_add_epi8(sum1, _mm256_srli_si256(sum1, 2));
    const __m256i sum3 = _mm256_add_epi8(sum2, _mm256_srli_si256(sum2, 4));
    const __m256i sum4 = _mm256_add_epi8(sum3, _mm256_srli_si256(sum3, 8));

    uint8_t tmp[32];
    _mm256_storeu_si256((__m256i *)&tmp, sum4);
    return tmp[0]+tmp[16];
}

static inline void _m256_psum_epi16_ovf(const __m256i& input, const uint32_t idx, uint32_t& pf_sum, uint32_t& idx_run) {

    __m256i halves[2];

    uint32_t tmp[8];
    const __m256i idx_vec =  _mm256_set1_epi32(idx);

    const auto low_bytes  = _mm256_extracti128_si256(input, 0);
    const auto high_bytes = _mm256_extracti128_si256(input, 1);

    halves[0] = _mm256_cvtepu16_epi32(low_bytes);
    halves[1] = _mm256_cvtepu16_epi32(high_bytes);

    //print32x8(halves[0]);
    //print32x8(halves[1]);

    halves[0] = _mm256_add_epi32(halves[0], _mm256_slli_si256(halves[0], 4));
    halves[0] = _mm256_add_epi32(halves[0], _mm256_slli_si256(halves[0], 8));
    _mm256_storeu_si256((__m256i *)&tmp, halves[0]);
    halves[0] = _mm256_add_epi32(halves[0], _mm256_set_m128i(_mm_set1_epi32(tmp[3]), _mm_set1_epi32(0)));

    __m256i idx_mask = _mm256_cmple_epu32(halves[0], idx_vec);//mask for >idx
    uint64_t lt_low = static_cast<uint32_t>(_mm256_movemask_epi8(idx_mask));//4 bits represent one element

    halves[1] = _mm256_add_epi32(halves[1], _mm256_set_epi32(0,0,0,0,0,0,0,(tmp[3]+tmp[7])));

    halves[1] = _mm256_add_epi32(halves[1], _mm256_slli_si256(halves[1], 4));
    halves[1] = _mm256_add_epi32(halves[1], _mm256_slli_si256(halves[1], 8));
    _mm256_storeu_si256((__m256i *)&tmp, halves[1]);
    halves[1] = _mm256_add_epi32(halves[1], _mm256_set_m128i(_mm_set1_epi32(tmp[3]), _mm_set1_epi32(0)));

    //print32x8(halves[0]);
    //print32x8(halves[1]);
    //print32x8(idx_vec);

    idx_mask = _mm256_cmple_epu32(halves[1], idx_vec);//4 bits represent one element
    const uint64_t lt_high = static_cast<uint32_t>(_mm256_movemask_epi8(idx_mask));//4 bits represent one element

    lt_low |= lt_high<<32;
    //print32x8(idx_mask);
    idx_run = __builtin_ctzll(~lt_low)>>2;//

    uint32_t pf_sum_vec[8];
    _mm256_storeu_si256((__m256i*)&pf_sum_vec, halves[idx_run>7]);

    pf_sum = pf_sum_vec[idx_run & 7];
}

static inline uint32_t _m256_hsum_epi16_ovf(const __m256i& input) {

    const __m256i sum = _mm256_add_epi32(_mm256_shuffle_epi8(input, _mm256_set_epi8(-1,-1,7,6, -1,-1,5,4, -1,-1,3,2, -1,-1,1,0,
                                                                                             -1,-1,7,6, -1,-1,5,4, -1,-1,3,2, -1,-1,1,0)),
                                         _mm256_shuffle_epi8(input, _mm256_set_epi8(-1,-1,15,14, -1,-1,13,12, -1,-1,11,10, -1,-1,9,8,
                                                                                             -1,-1,15,14, -1,-1,13,12, -1,-1,11,10, -1,-1,9,8)));
    //print32x8(sum);
    const __m256i sum1 = _mm256_add_epi32(sum, _mm256_srli_si256(sum, 4));
    const __m256i sum2 = _mm256_add_epi32(sum1, _mm256_srli_si256(sum1, 8));
    //print32x8(sum2);
    uint32_t tmp[8];
    _mm256_storeu_si256((__m256i *)&tmp, sum2);
    return tmp[0]+tmp[4];
}

static inline uint16_t _m256_hsum_epi16(const __m256i& input) {
    const __m256i sum1 = _mm256_add_epi16(input, _mm256_srli_si256(input, 2));
    const __m256i sum2 = _mm256_add_epi16(sum1, _mm256_srli_si256(sum1, 4));
    const __m256i sum3 = _mm256_add_epi16(sum2, _mm256_srli_si256(sum2, 8));
    uint16_t tmp[16];
    _mm256_storeu_si256((__m256i *)&tmp, sum3);
    return tmp[0]+tmp[8];
}

static inline uint32_t _m256_hsum_epi32(const __m256i& input) {
    const __m256i sum1 = _mm256_add_epi32(input, _mm256_srli_si256(input, 4));
    const __m256i sum2 = _mm256_add_epi32(sum1, _mm256_srli_si256(sum1, 8));
    uint32_t tmp[8];
    _mm256_storeu_si256((__m256i *)&tmp, sum2);
    return tmp[0]+tmp[4];
}

static inline uint64_t _m256_hsum_epi64(const __m256i& input) {
    const __m256i sum = _mm256_add_epi32(input, _mm256_srli_si256(input, 8));
    uint64_t tmp[4];
    _mm256_storeu_si256((__m256i *)&tmp, sum);
    return tmp[0]+tmp[2];
}

template<bool vbyte_compressed, uint8_t ctr_width, uint8_t bytes_per_run>
static inline __m256i decode_block_avx2(const uint8_t **stream) {

    if constexpr (vbyte_compressed) {

        uint8_t *pshuf_low, *pshuf_high;
        const uint8_t *ctrl_bits = *stream;
        uint8_t len1, len2;

        if constexpr (ctr_width == 1) {
            pshuf_low = (uint8_t *) &dec_table_16x8[*ctrl_bits];
            len1 = pshuf_low[14 + (*ctrl_bits >> 7)] + 1;
            ctrl_bits+=len1+1;
            pshuf_high = (uint8_t *) &dec_table_16x8[*ctrl_bits];
            len2 = pshuf_high[14 + (*ctrl_bits >> 7)] + 1;
        } else if constexpr (ctr_width == 2) {
            pshuf_low = (uint8_t *) &dec_table_32x4[*ctrl_bits];
            len1 = pshuf_low[12 + (*ctrl_bits >> 6)] + 1;
            ctrl_bits+=len1+1;
            pshuf_high = (uint8_t *) &dec_table_32x4[*ctrl_bits];
            len2 = pshuf_high[12 + (*ctrl_bits >> 6)] + 1;
        } else if constexpr (ctr_width == 3) {
            pshuf_low = (uint8_t *) &dec_table_64x2[*ctrl_bits];
            len1 = pshuf_low[8 + ((*ctrl_bits >> 3) & 7)] + 1;
            ctrl_bits+=len1+1;
            pshuf_high = (uint8_t *) &dec_table_64x2[*ctrl_bits];
            len2 = pshuf_high[8 + ((*ctrl_bits >> 3) & 7)] + 1;
        } else {
            exit(1);
        }

        const __m128i comp_low = _mm_loadu_si128((const __m128i*)(*stream+1));
        const __m128i comp_high = _mm_loadu_si128((const __m128i*)(*stream+1+len1+1));
        const auto compressed = _mm256_set_m128i(comp_high, comp_low);

        const __m128i dec_shuffle_low = _mm_loadu_si128((const __m128i*)pshuf_low);
        const __m128i dec_shuffle_high = _mm_loadu_si128((const __m128i*)pshuf_high);
        const auto dec_shuffle = _mm256_set_m128i(dec_shuffle_high, dec_shuffle_low);

        __m256i data = _mm256_shuffle_epi8(compressed, dec_shuffle);
        *stream += len1 + len2 + 2;
        return data;
    } else {

        //align the bytes
        if constexpr (bytes_per_run==3){

            const __m128i data_low = _mm_loadu_si128((const __m128i*)(*stream));
            const __m128i data_high = _mm_loadu_si128((const __m128i*)(*stream+12));
            auto data = _mm256_set_m128i(data_high, data_low);

            const __m256i dec_shuff = _mm256_set_epi8(-1,11,10,9, -1,8,7,6, -1,5,4,3, -1,2,1,0,
                                                      -1,11,10,9, -1,8,7,6, -1,5,4,3, -1,2,1,0);
            data = _mm256_shuffle_epi8(data, dec_shuff);
            *stream+=24;
            return data;

        } else if constexpr (bytes_per_run==5){

            const __m128i data_low = _mm_loadu_si128((const __m128i*)(*stream));
            const __m128i data_high = _mm_loadu_si128((const __m128i*)(*stream+10));
            auto data = _mm256_set_m128i(data_high, data_low);

            const __m256i dec_shuff = _mm256_set_epi8(-1,-1,-1,9,8,7,6,5, -1,-1,-1,4,3,2,1,0,
                                                      -1,-1,-1,9,8,7,6,5, -1,-1,-1,4,3,2,1,0);
            data = _mm256_shuffle_epi8(data, dec_shuff);
            *stream+=20;
            return data;

        } else if constexpr (bytes_per_run==6){

            const __m128i data_low = _mm_loadu_si128((const __m128i*)(*stream));
            const __m128i data_high = _mm_loadu_si128((const __m128i*)(*stream+12));
            auto data = _mm256_set_m128i(data_high, data_low);

            const __m256i dec_shuff = _mm256_set_epi8(-1,-1,11,10,9,8,7,6, -1,-1,5,4,3,2,1,0,
                                                      -1,-1,11,10,9,8,7,6, -1,-1,5,4,3,2,1,0);
            data = _mm256_shuffle_epi8(data, dec_shuff);
            *stream+=24;
            return data;

        } else if constexpr (bytes_per_run==7){

            const __m128i data_low = _mm_loadu_si128((const __m128i*)(*stream));
            const __m128i data_high = _mm_loadu_si128((const __m128i*)(*stream+14));
            auto data = _mm256_set_m128i(data_high, data_low);

            const __m256i dec_shuff = _mm256_set_epi8(-1,13,12,11,10,9,8,7, -1,6,5,4,3,2,1,0,
                                                      -1,13,12,11,10,9,8,7, -1,6,5,4,3,2,1,0);
            data = _mm256_shuffle_epi8(data, dec_shuff);
            *stream+=28;
            return data;
        } else {
            const __m256i data = _mm256_loadu_si256((const __m256i*)(*stream));
            *stream+=32;
            return data;
        }
    }
}

template<bool overflow16, bool overflow32=false>
static inline std::pair<uint64_t, uint64_t> get_phi_run_avx2_8x32(const uint8_t **stream, uint64_t idx){

    //NOTE here I do not need to vbyte compress the block
    __m256i block = _mm256_loadu_si256((const __m256i*)*stream);
    *stream+=32;

    uint32_t idx_run = 0;

    //sometimes the back of the block has garbage, so I have to assume overflow
    uint32_t prev_acc = 0;
    uint32_t acc = _m256_hsum_epi8_ovf(block);

    while(acc<=idx){
        block = _mm256_loadu_si256((const __m256i*)*stream);
        *stream+=32;
        idx_run+=32;

        prev_acc = acc;
        acc += _m256_hsum_epi8_ovf(block);
    }

    idx-=prev_acc;

    uint32_t pf_sum, tmp_idx_run;
    if constexpr(overflow16 || overflow32) {
        _m256_psum_epi8_ovf(block, idx, pf_sum, tmp_idx_run);
    }else{
        uint8_t tmp[32];
        block = _mm256_add_epi8(block, _mm256_slli_si256(block, 1));
        block = _mm256_add_epi8(block, _mm256_slli_si256(block, 2));
        block = _mm256_add_epi8(block, _mm256_slli_si256(block, 4));
        block = _mm256_add_epi8(block, _mm256_slli_si256(block, 8));
        _mm256_storeu_si256((__m256i *)&tmp, block);//this is due to limitations in the shift in 256-bit lanes
        block = _mm256_add_epi16(block, _mm256_set_m128i(_mm_set1_epi8(tmp[15]), _mm_set1_epi8(0)));
        //

        const __m256i idx_mask = _mm256_cmple_epu8(block, _mm256_set1_epi8(idx));//mask for >idx
        const uint32_t less_than = _mm256_movemask_epi8(idx_mask);
        idx_run = __builtin_ctzll(~less_than);
        pf_sum = tmp[idx_run] + tmp[15]*(idx_run>15);
    }

    const uint32_t offset = static_cast<uint32_t>((*stream - 32)[tmp_idx_run]) - (pf_sum-idx);
    return std::make_pair(idx_run+tmp_idx_run, offset);
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline std::pair<uint64_t, uint64_t> get_phi_run_avx2_16x16(const uint8_t **stream, uint64_t idx){

    __m256i block = decode_block_avx2<vbyte_compressed, 1, 2>(stream);
    uint32_t prev_acc = 0;

    //sometimes the block has garbage at the end, so we have to assume overflow at the end
    uint32_t acc = _m256_hsum_epi16_ovf(block);
    uint32_t idx_run = 0;

    while(acc<=idx){
        block = decode_block_avx2<vbyte_compressed,1,2>(stream);
        idx_run+=16;
        //print16x8(bk_lengths);
        prev_acc = acc;
        acc += _m256_hsum_epi16_ovf(block);
    }

    idx-=prev_acc;
    uint32_t pf_sum, tmp_idx_run;
    uint16_t tmp[16];
    if constexpr (overflow8 || overflow16){
        _m256_psum_epi16_ovf(block, idx, pf_sum, tmp_idx_run);
    }else {
        //prefix sum without overflow
        __m256i pf_sum_vec = _mm256_add_epi16(block, _mm256_slli_si256(block, 2));
        pf_sum_vec = _mm256_add_epi16(pf_sum_vec, _mm256_slli_si256(pf_sum_vec, 4));
        pf_sum_vec = _mm256_add_epi16(pf_sum_vec, _mm256_slli_si256(pf_sum_vec, 8));
        _mm256_storeu_si256((__m256i *)&tmp, pf_sum_vec);
        pf_sum_vec = _mm256_add_epi16(pf_sum_vec, _mm256_set_m128i(_mm_set1_epi16(tmp[7]), _mm_set1_epi16(0)));
        const __m256i idx_mask = _mm256_cmple_epu16(pf_sum_vec, _mm256_set1_epi16(idx));//mask for <=idx

        uint32_t less_than = _mm256_movemask_epi8(_mm256_packs_epi16(idx_mask, _mm256_setzero_si256()));
        less_than = (less_than >> 8) | (less_than & 0xFF);//this is a hack because the way _mm256_packs_epi16 works
        idx_run = __builtin_ctzll(~less_than);
        pf_sum = tmp[idx_run] + tmp[7]*(idx_run>7);
    }

    _mm256_storeu_si256((__m256i *)&tmp, block);
    const uint16_t len = tmp[idx_run];

    uint64_t offset = len - (pf_sum-idx);
    return std::make_pair(idx_run+tmp_idx_run, offset);
}

template<bool overflow16, bool overflow32=false, bool check_head>
static inline int64_t rank_avx2_8x32(const uint8_t **stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m256i alpha_mask = _mm256_set1_epi8(alpha_m);
    const __m256i sym_vec = _mm256_set1_epi8(symbol);

    __m256i block = _mm256_loadu_si256((const __m256i*)*stream);
    *stream+=32;
    __m256i bk_lengths =  _m256_shift_right_epi8(block, sigma_bits);

    uint32_t prev_acc = 0;

    //the back of the block *might* contain garbage, so I have to assume overflow
    uint32_t acc = _m256_hsum_epi8_ovf(bk_lengths);

    uint32_t rank=0;
    while(acc<=idx){
        //compute acc rank in the previous block
        bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi8(_mm256_and_si256(block, alpha_mask), sym_vec));
        if constexpr (overflow16 || overflow32){
            rank += _m256_hsum_epi8_ovf(bk_lengths);
        } else {
            rank += _m256_hsum_epi8(bk_lengths);
        }

        block = _mm256_loadu_si256((const __m256i*)*stream);
        *stream+=32;
        bk_lengths = _m256_shift_right_epi8(block, sigma_bits);
        prev_acc = acc;
        acc += _m256_hsum_epi8_ovf(bk_lengths);
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr(overflow16 || overflow32) {
        _m256_psum_epi8_ovf(bk_lengths, idx, pf_sum, idx_run);
    } else {
        //prefix sum without overflow
        uint8_t tmp[32];
        bk_lengths = _mm256_add_epi8(bk_lengths, _mm256_slli_si256(bk_lengths, 1));
        bk_lengths = _mm256_add_epi8(bk_lengths, _mm256_slli_si256(bk_lengths, 2));
        bk_lengths = _mm256_add_epi8(bk_lengths, _mm256_slli_si256(bk_lengths, 4));
        bk_lengths = _mm256_add_epi8(bk_lengths, _mm256_slli_si256(bk_lengths, 8));
        _mm256_storeu_si256((__m256i *)&tmp, bk_lengths);//this is due to limitations in the shift in 256-bit lanes
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_set_m128i(_mm_set1_epi8(tmp[15]), _mm_set1_epi8(0)));
        //

        const __m256i idx_mask = _mm256_cmple_epu8(bk_lengths, _mm256_set1_epi8(idx));//mask for >idx
        const uint32_t less_than = _mm256_movemask_epi8(idx_mask);
        idx_run = __builtin_ctzll(~less_than);
        pf_sum = tmp[idx_run] + tmp[15]*(idx_run>15);
    }

    const uint8_t run = (*stream-32)[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    //create a mask for the index position
    //eg: idx_run = 3 yields {0:0xFF, 1:0xFF, 2:0xFF, 3:0xFF, 4:0, 5:0, ...}
    static const __m256i indices = _mm256_set_epi8(31,30,29,28,27,26,25,24,
                                                   23,22,21,20,19,18,17,16,
                                                   15,14,13,12,11,10,9,8,
                                                   7,6,5,4,3,2,1,0);
    const __m256i idx_run_mask = _mm256_cmpgt_epi8( _mm256_set1_epi8(idx_run+1), indices);
    //

    bk_lengths =  _m256_shift_right_epi8(block, sigma_bits);
    bk_lengths = _mm256_and_si256(bk_lengths, idx_run_mask);
    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi8(_mm256_and_si256(block, alpha_mask), sym_vec));

    if constexpr (overflow16 || overflow32){
        rank += _m256_hsum_epi8_ovf(bk_lengths);
    }else{
        rank += _m256_hsum_epi8(bk_lengths);
    }

    if constexpr (check_head) {
        const uint8_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx;//check if idx is the head of the run
        rank -=(pf_sum-idx) * is_same_sym;
        rank = (rank<<1) | is_head;
        rank = (rank<<1) | is_same_sym;
    } else {
        rank -=(pf_sum-idx) * (last_symbol==symbol);
    }
    return rank;
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false, bool check_head>
static inline int64_t rank_avx2_16x16(const uint8_t **stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m256i alpha_mask = _mm256_set1_epi16(alpha_m);
    const __m256i sym_vec = _mm256_set1_epi16(symbol);

    __m256i block = decode_block_avx2<vbyte_compressed, 1, 2>(stream);
    __m256i bk_lengths =  _m256_shift_right_epi16(block, sigma_bits);

    //the last block might have some garbage, so we have to assume overflow at the end
    uint32_t prev_acc = 0;
    uint32_t acc = _m256_hsum_epi16_ovf(bk_lengths);
    uint64_t rank = 0;

    while(acc<=idx){
        bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi16(_mm256_and_si256(block, alpha_mask), sym_vec));
        if constexpr (overflow8 || overflow16){
            rank += _m256_hsum_epi16_ovf(bk_lengths);
        }else{
            rank += _m256_hsum_epi16(bk_lengths);
        }

        block = decode_block_avx2<vbyte_compressed,1,2>(stream);

        bk_lengths =  _m256_shift_right_epi16(block, sigma_bits);

        prev_acc = acc;
        acc += _m256_hsum_epi16_ovf(bk_lengths);
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;
    uint16_t tmp[16];
    if constexpr (overflow8 || overflow16){
        _m256_psum_epi16_ovf(bk_lengths, idx, pf_sum, idx_run);
    }else{
        //prefix sum without overflow
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_slli_si256(bk_lengths, 2));
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_slli_si256(bk_lengths, 4));
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_slli_si256(bk_lengths, 8));
        _mm256_storeu_si256((__m256i *)&tmp, bk_lengths);
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_set_m128i(_mm_set1_epi16(tmp[7]), _mm_set1_epi16(0)));
        const __m256i idx_mask = _mm256_cmple_epu16(bk_lengths, _mm256_set1_epi16(idx));//mask for <=idx
        uint32_t less_than = _mm256_movemask_epi8(_mm256_packs_epi16(idx_mask, _mm256_setzero_si256()));
        less_than = (less_than >> 8) | (less_than & 0xFF);//this is a hack because the way _mm256_packs_epi16 works
        idx_run = __builtin_ctzll(~less_than);
        pf_sum = tmp[idx_run] + tmp[7]*(idx_run>7);
    }

    _mm256_storeu_si256((__m256i *)&tmp, block);
    const uint16_t run = tmp[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    //create a mask for the index position
    //ef: idx_run = 3 is equal of {0:0xFFFF, 1:0xFFFF, 2:0xFFFF, 3:0xFFFF, 4:0, 5:0, ...}
    static const __m256i indices = _mm256_set_epi16(15,14,13,12,
                                                    11,10,9,8,
                                                    7,6,5,4,
                                                    3,2,1,0);
    const __m256i idx_run_mask = _mm256_cmpgt_epi16( _mm256_set1_epi16(idx_run+1), indices);
    //

    bk_lengths =  _m256_shift_right_epi16(block, sigma_bits);
    bk_lengths = _mm256_and_si256(bk_lengths, idx_run_mask);
    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi16(_mm256_and_si256(block, alpha_mask), sym_vec));

    if constexpr (overflow8 || overflow16){
        rank += _m256_hsum_epi16_ovf(bk_lengths);
    }else{
        rank += _m256_hsum_epi16(bk_lengths);
    }

    if constexpr (check_head) {
        const uint16_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx;//check if idx is the head of the run
        rank -= (pf_sum-idx) * is_same_sym;
        rank = (rank<<1) | is_head;
        rank = (rank<<1) | is_same_sym;
    } else {
        rank -= (pf_sum-idx) * (last_symbol==symbol);
    }
    return (int64_t)rank;
}

template<bool vbyte_compressed, uint8_t bytes_per_run, bool check_head>
static inline int64_t rank_avx2_32x8(const uint8_t ** stream, const uint8_t sigma, uint64_t idx, uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m256i alpha_mask = _mm256_set1_epi32(alpha_m);
    const __m256i sym_vec = _mm256_set1_epi32(symbol);

    __m256i block = decode_block_avx2<vbyte_compressed, 2, bytes_per_run>(stream);
    __m256i bk_lengths =  _m256_shift_right_epi32(block, sigma_bits);

    uint32_t prev_acc=0;
    uint64_t acc= _m256_hsum_epi32(bk_lengths);
    uint64_t rank = 0;

    while(acc<=idx) {
        bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi32(_mm256_and_si256(block, alpha_mask), sym_vec));
        rank += _m256_hsum_epi32(bk_lengths);

        block = decode_block_avx2<vbyte_compressed, 2, bytes_per_run>(stream);
        bk_lengths =  _m256_shift_right_epi32(block, sigma_bits);

        prev_acc = acc;
        acc += _m256_hsum_epi32(bk_lengths);
    }

    idx-=prev_acc;

    uint32_t tmp[8];
    bk_lengths = _mm256_add_epi32(bk_lengths, _mm256_slli_si256(bk_lengths, 4));
    bk_lengths = _mm256_add_epi32(bk_lengths, _mm256_slli_si256(bk_lengths, 8));
    _mm256_storeu_si256((__m256i *)&tmp, bk_lengths);
    bk_lengths = _mm256_add_epi32(bk_lengths,_mm256_set_m128i(_mm_set1_epi32(tmp[3]), _mm_set1_epi32(0)));

    const __m256i idx_mask = _mm256_cmple_epu32(bk_lengths, _mm256_set1_epi32(idx));//mask for >idx
    const int less_than = _mm256_movemask_ps(_mm256_castsi256_ps(idx_mask));
    const uint8_t idx_run = __builtin_ctz(~less_than);
    const uint32_t pf_sum = tmp[idx_run] + tmp[3]*(idx_run>3);

    _mm256_storeu_si256((__m256i *)&tmp, block);
    const uint32_t run = tmp[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    //create a mask for the index position
    //eg: idx_run = 3 produces {0:0xFFFFFFFF, 1:0xFFFFFFFF, 2:0xFFFFFFFF, 3:0xFFFFFFFF, 4:0, 5:0, ...}
    static const __m256i indices = _mm256_set_epi32(7,6,5,4, 3,2,1,0);
    const __m256i idx_run_mask = _mm256_cmpgt_epi32( _mm256_set1_epi32(idx_run+1), indices);
    //
    bk_lengths =  _m256_shift_right_epi32(block, sigma_bits);
    bk_lengths = _mm256_and_si256(bk_lengths, idx_run_mask);
    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi32(_mm256_and_si256(block, alpha_mask), sym_vec));

    rank += _m256_hsum_epi32(bk_lengths);

    if constexpr (check_head) {
        const uint32_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx;//check if idx is the head of the run
        rank -= (pf_sum-idx) * is_same_sym;
        rank = (rank<<1) | is_head;
        rank = (rank<<1) | is_same_sym;
    } else {
        rank -= (pf_sum-idx)*(last_symbol==symbol);
    }
    return (int64_t)rank;
}

template<uint8_t bytes_per_run, bool check_head>
static inline int64_t rank_avx2_64x4(const uint8_t ** stream, uint8_t sigma, uint64_t idx, uint8_t symbol){
    return 0;
}

template<bool overflow16, bool overflow32=false, bool check_head>
static inline std::pair<int64_t, int64_t> range_rank_avx2_8x32(const uint8_t **stream, const uint8_t sigma, uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m256i alpha_mask = _mm256_set1_epi8(alpha_m);
    const __m256i sym_vec = _mm256_set1_epi8(symbol);

    __m256i block = _mm256_loadu_si256((const __m256i*)*stream);
    *stream+=32;
    __m256i bk_lengths =  _m256_shift_right_epi8(block, sigma_bits);

    uint32_t prev_acc = 0;

    //print8x32(block);
    //print8x32(_mm256_and_si256(block, alpha_mask));
    //print8x32(bk_lengths);

    //the back of the block *might* contain garbage, so I have to assume overflow
    uint32_t acc = _m256_hsum_epi8_ovf(bk_lengths);

    uint32_t rank=0;
    size_t l=0;
    while(acc<=idx_i){
        //compute acc rank in the previous block
        bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi8(_mm256_and_si256(block, alpha_mask), sym_vec));
        if constexpr (overflow16 || overflow32){
            rank += _m256_hsum_epi8_ovf(bk_lengths);
        } else {
            rank += _m256_hsum_epi8(bk_lengths);
        }

        block = _mm256_loadu_si256((const __m256i*)*stream);
        *stream+=32;
        bk_lengths = _m256_shift_right_epi8(block, sigma_bits);
        //print8x16(bk_lengths);
        prev_acc = acc;
        acc += _m256_hsum_epi8_ovf(bk_lengths);
        l++;

        //print8x32(_mm256_and_si256(block, alpha_mask));
        //print8x32(bk_lengths);
    }

    idx_i-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr(overflow16 || overflow32) {
        _m256_psum_epi8_ovf(bk_lengths, idx_i, pf_sum, idx_run);
    } else {
        //prefix sum without overflow
        uint8_t tmp[32];
        bk_lengths = _mm256_add_epi8(bk_lengths, _mm256_slli_si256(bk_lengths, 1));
        bk_lengths = _mm256_add_epi8(bk_lengths, _mm256_slli_si256(bk_lengths, 2));
        bk_lengths = _mm256_add_epi8(bk_lengths, _mm256_slli_si256(bk_lengths, 4));
        bk_lengths = _mm256_add_epi8(bk_lengths, _mm256_slli_si256(bk_lengths, 8));
        _mm256_storeu_si256((__m256i *)&tmp, bk_lengths);//this is due to limitations in the shift in 256-bit lanes
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_set_m128i(_mm_set1_epi8(tmp[15]), _mm_set1_epi8(0)));
        //

        const __m256i idx_mask = _mm256_cmple_epu8(bk_lengths, _mm256_set1_epi8(idx_i));//mask for >idx
        const uint32_t less_than = _mm256_movemask_epi8(idx_mask);
        idx_run = __builtin_ctzll(~less_than);
        pf_sum = tmp[idx_run] + tmp[15]*(idx_run>15);
    }

    *stream -= 32;
    const uint8_t run = (*stream)[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    block = _mm256_loadu_si256((const __m256i*)*stream);
    bk_lengths =  _m256_shift_right_epi8(block, sigma_bits);
    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_loadu_si256((const __m256i*)mask8x32[idx_run+1]));
    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi8(_mm256_and_si256(block, alpha_mask), sym_vec));

    if constexpr (overflow16 || overflow32){
        rank += _m256_hsum_epi8_ovf(bk_lengths);
    }else{
        rank += _m256_hsum_epi8(bk_lengths);
    }

    if constexpr (check_head) {
        const uint8_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank -=(pf_sum-idx_i) * is_same_sym;
        rank = (rank<<1) | is_head;
        rank = (rank<<1) | is_same_sym;
    } else {
        rank -=(pf_sum-idx_i) * (last_symbol==symbol);
    }

    return std::make_pair(rank, rank);
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false, bool check_head>
static inline int64_t range_rank_avx2_16x16(const uint8_t **stream, const uint8_t sigma, uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m256i alpha_mask = _mm256_set1_epi16(alpha_m);
    const __m256i sym_vec = _mm256_set1_epi16(symbol);

    const uint8_t *prev_state = *stream;

    //TODO testing
    //__m128i b = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    //print16x8(b);
    //b = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    //print16x8(b);
    //

    *stream = prev_state;
    __m256i block = decode_block_avx2<vbyte_compressed, 1, 2>(stream);
    //print16x16(block);
    __m256i bk_lengths =  _m256_shift_right_epi16(block, sigma_bits);
    //print16x16(_mm256_and_si256(block, alpha_mask));
    //print16x16(bk_lengths);

    //the last block might have some garbage, so we have to assume overflow at the end
    uint32_t prev_acc = 0;
    uint32_t acc = _m256_hsum_epi16_ovf(bk_lengths);
    uint64_t rank = 0;

    while(acc<=idx_i){
        bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi16(_mm256_and_si256(block, alpha_mask), sym_vec));
        if constexpr (overflow8 || overflow16){
            rank += _m256_hsum_epi16_ovf(bk_lengths);
        }else{
            rank += _m256_hsum_epi16(bk_lengths);
        }

        prev_state = *stream;
        block = decode_block_avx2<vbyte_compressed,1,2>(stream);
        //print16x16(block);

        bk_lengths =  _m256_shift_right_epi16(block, sigma_bits);
        //print16x16(_mm256_and_si256(block, alpha_mask));
        //print16x16(bk_lengths);

        prev_acc = acc;
        acc += _m256_hsum_epi16_ovf(bk_lengths);
    }

    idx_i-=prev_acc;
    uint32_t idx_run, pf_sum;
    uint16_t tmp[16];
    if constexpr (overflow8 || overflow16){
        _m256_psum_epi16_ovf(bk_lengths, idx_i, pf_sum, idx_run);
    }else{
        //prefix sum without overflow
        //print16x16(bk_lengths);
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_slli_si256(bk_lengths, 2));
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_slli_si256(bk_lengths, 4));
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_slli_si256(bk_lengths, 8));
        _mm256_storeu_si256((__m256i *)&tmp, bk_lengths);
        bk_lengths = _mm256_add_epi16(bk_lengths, _mm256_set_m128i(_mm_set1_epi16(tmp[7]), _mm_set1_epi16(0)));
        //print16x16(bk_lengths);
        const __m256i idx_mask = _mm256_cmple_epu16(bk_lengths, _mm256_set1_epi16(idx_i));//mask for <=idx
        uint32_t less_than = _mm256_movemask_epi8(_mm256_packs_epi16(idx_mask, _mm256_setzero_si256()));
        less_than = (less_than >> 8) | (less_than & 0xFF);//this is a hack because the way _mm256_packs_epi16 works
        idx_run = __builtin_ctzll(~less_than);
        pf_sum = tmp[idx_run] + tmp[7]*(idx_run>7);
    }

    _mm256_storeu_si256((__m256i *)&tmp, block);
    const uint16_t run = tmp[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    *stream = prev_state;
    block = decode_block_avx2<vbyte_compressed,1,2>(stream);

    //print16x16(block);
    bk_lengths =  _m256_shift_right_epi16(block, sigma_bits);
    //print16x16(bk_lengths);
    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_loadu_si256((const __m256i*)mask16x16[idx_run+1]));

    //print16x16(bk_lengths);
    //print16x16(alpha_mask);
    //print16x16(sym_vec);
    //print16x16(_mm256_and_si256(block, alpha_mask));

    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi16(_mm256_and_si256(block, alpha_mask), sym_vec));
    //print16x16(bk_lengths);

    if constexpr (overflow8 || overflow16){
        rank += _m256_hsum_epi16_ovf(bk_lengths);
    }else{
        rank += _m256_hsum_epi16(bk_lengths);
    }

    if constexpr (check_head) {
        const uint16_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank -= (pf_sum-idx_i) * is_same_sym;
        rank = (rank<<1) | is_head;
        rank = (rank<<1) | is_same_sym;
    } else {
        rank -= (pf_sum-idx_i) * (last_symbol==symbol);
    }
    return (int64_t)rank;
}

template<bool vbyte_compressed, uint8_t bytes_per_run, bool check_head>
static inline std::pair<int64_t, int64_t> range_rank_avx2_32x8(const uint8_t ** stream, const uint8_t sigma, uint64_t idx_i, uint64_t idx_j, uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m256i alpha_mask = _mm256_set1_epi32(alpha_m);
    const __m256i sym_vec = _mm256_set1_epi32(symbol);

    const uint8_t *prev_state = *stream;
    __m256i block = decode_block_avx2<vbyte_compressed, 2, bytes_per_run>(stream);
    //print32x8(block);
    __m256i bk_lengths =  _m256_shift_right_epi32(block, sigma_bits);
    //print32x8(bk_lengths);

    uint32_t prev_acc=0;
    uint64_t acc= _m256_hsum_epi32(bk_lengths);
    uint64_t rank = 0;

    while(acc<=idx_i) {
        bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi32(_mm256_and_si256(block, alpha_mask), sym_vec));
        //print32x8(bk_lengths);
        rank += _m256_hsum_epi32(bk_lengths);

        prev_state = *stream;
        block = decode_block_avx2<vbyte_compressed, 2, bytes_per_run>(stream);
        //print32x8(block);
        bk_lengths =  _m256_shift_right_epi32(block, sigma_bits);
        //print32x8(bk_lengths);

        prev_acc = acc;
        acc += _m256_hsum_epi32(bk_lengths);
    }

    idx_i-=prev_acc;

    uint32_t tmp[8];
    bk_lengths = _mm256_add_epi32(bk_lengths, _mm256_slli_si256(bk_lengths, 4));
    bk_lengths = _mm256_add_epi32(bk_lengths, _mm256_slli_si256(bk_lengths, 8));
    _mm256_storeu_si256((__m256i *)&tmp, bk_lengths);
    bk_lengths = _mm256_add_epi32(bk_lengths, _mm256_set_m128i(_mm_set1_epi32(tmp[3]), _mm_set1_epi32(0)));

    //print32x8(bk_lengths);
    const __m256i idx_mask = _mm256_cmple_epu32(bk_lengths, _mm256_set1_epi32(idx_i));//mask for >idx
    const int less_than = _mm256_movemask_ps(_mm256_castsi256_ps(idx_mask));
    const uint8_t idx_run = __builtin_ctz(~less_than);
    const uint32_t pf_sum = tmp[idx_run] + tmp[3]*(idx_run>3);

    _mm256_storeu_si256((__m256i *)&tmp, block);
    const uint32_t run = tmp[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    *stream = prev_state;
    block = decode_block_avx2<vbyte_compressed, 2, bytes_per_run>(stream);
    //print32x4(block);
    bk_lengths =  _m256_shift_right_epi32(block, sigma_bits);
    //print32x4(bk_lengths);
    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_loadu_si256((const __m256i*)mask32x8[idx_run+1]));
    //print32x4(bk_lengths);
    bk_lengths = _mm256_and_si256(bk_lengths, _mm256_cmpeq_epi32(_mm256_and_si256(block, alpha_mask), sym_vec));
    //print32x4(bk_lengths);

    rank += _m256_hsum_epi32(bk_lengths);

    if constexpr (check_head) {
        const uint32_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank -= (pf_sum-idx_i) * is_same_sym;
        rank = (rank<<1) | is_head;
        rank = (rank<<1) | is_same_sym;
    } else {
        rank -= (pf_sum-idx_i)*(last_symbol==symbol);
    }
    return std::make_pair(rank, rank);
}

template<uint8_t bytes_per_run, bool check_head>
static inline std::pair<int64_t, int64_t> range_rank_avx2_64x4(const uint8_t ** stream, uint8_t sigma, uint64_t idx_i, size_t idx_j, uint8_t symbol){
    return {0,0};
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
        //for(size_t j=0;j<16;j++){
        //    std::cout<<i<<" = "<<idx<<" -> "<<j<<" / "<<acc[j]<<std::endl;
        //}

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
        // _mm256_storeu_si256((__m256i *)&acc, pf_sum);
        //for(size_t j=0;j<16;j++){
        //    std::cout<<idx<<" -> "<<j<<" / "<<acc[j]<<std::endl;
        //}
        //_mm256_storeu_si256((__m256i *)&acc, counter);
        //for(size_t j=0;j<16;j++){
        //    std::cout<<j<<" | "<<acc[j]<<std::endl;
        //}
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
