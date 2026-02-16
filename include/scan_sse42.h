/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */

#ifndef VLBT_SCAN_SSE42_H
#define VLBT_SCAN_SSE42_H

#include <x86intrin.h>
#include <bits/stdint-uintn.h>

#include "simd_tables.h"
#include "utils.h"

#define _mm_cmpge_epu8(a, b) \
        _mm_cmpeq_epi8(_mm_max_epu8(a, b), a)
#define _mm_cmple_epu8(a, b) _mm_cmpge_epu8(b, a)

#define _mm_cmpge_epu16(a, b) \
        _mm_cmpeq_epi16(_mm_max_epu16(a, b), a)
#define _mm_cmple_epu16(a, b) _mm_cmpge_epu16(b, a)

#define _mm_cmpge_epu32(a, b) \
        _mm_cmpeq_epi32(_mm_max_epu32(a, b), a)
#define _mm_cmple_epu32(a, b) _mm_cmpge_epu32(b, a)

#define _mm_cmpge_epu64(a, b) \
        _mm_cmpeq_epi64(_mm_max_epu64(a, b), a)
#define _mm_cmple_epu64(a, b) _mm_cmpge_epu64(b, a)

static void print8x16(__m128i vec) {
    uint8_t values[16];
    _mm_storeu_si128((__m128i*)values, vec); // Unaligned store
    for (int i = 0; i < 16; ++i) {
        std::cout<<int(values[i])<<" ";
    }
    std::cout<<""<<std::endl;
}

static void print16x8(__m128i vec) {
    uint16_t values[8];
    _mm_storeu_si128((__m128i*)values, vec); // Unaligned store
    for (int i = 0; i < 8; ++i) {
        std::cout<<int(values[i])<<" ";
    }
    std::cout<<""<<std::endl;
}

static void print32x4(__m128i vec) {
    uint32_t values[4];
    _mm_storeu_si128((__m128i*)values, vec); // Unaligned store
    for (int i = 0; i < 4; ++i) {
        std::cout<<values[i]<<" ";
    }
    std::cout<<""<<std::endl;
}

static void print64x2(__m128i vec) {
    uint64_t values[2];
    _mm_storeu_si128((__m128i*)values, vec); // Unaligned store
    for (int i = 0; i < 2; ++i) {
        std::cout<<values[i]<<" ";
    }
    std::cout<<""<<std::endl;
}

static __m128i shift_right_epi8(const __m128i& input, const uint8_t shift) {
    //from https://wunkolo.github.io/post/2020/11/gf2p8affineqb-int8-shifting/
    switch (shift) {
        case 0: return input;
        case 1: return _mm_and_si128(_mm_srli_epi64(input, 1), _mm_set1_epi8(127));
        case 2: return _mm_and_si128(_mm_srli_epi64(input, 2), _mm_set1_epi8(63));
        case 3: return _mm_and_si128(_mm_srli_epi64(input, 3), _mm_set1_epi8(31));
        case 4: return _mm_and_si128(_mm_srli_epi64(input, 4), _mm_set1_epi8(15));
        case 5: return _mm_and_si128(_mm_srli_epi64(input, 5), _mm_set1_epi8(7));
        case 6: return _mm_and_si128(_mm_srli_epi64(input, 6), _mm_set1_epi8(3));
        case 7: return _mm_and_si128(_mm_srli_epi64(input, 7), _mm_set1_epi8(1));
        default: return _mm_setzero_si128(); // for shift >= 8
    }
}

static __m128i shift_right_epi16(const __m128i& input, const int shift) {
    switch (shift) {
        case 0: return input;
        case 1: return _mm_srli_epi16(input, 1);
        case 2: return _mm_srli_epi16(input, 2);
        case 3: return _mm_srli_epi16(input, 3);
        case 4: return _mm_srli_epi16(input, 4);
        case 5: return _mm_srli_epi16(input, 5);
        case 6: return _mm_srli_epi16(input, 6);
        case 7: return _mm_srli_epi16(input, 7);
        default: return _mm_setzero_si128(); // for shift >= 8
    }
}

static __m128i shift_right_epi32(const __m128i& input, const int shift) {
    switch (shift) {
        case 0: return input;
        case 1: return _mm_srli_epi32(input, 1);
        case 2: return _mm_srli_epi32(input, 2);
        case 3: return _mm_srli_epi32(input, 3);
        case 4: return _mm_srli_epi32(input, 4);
        case 5: return _mm_srli_epi32(input, 5);
        case 6: return _mm_srli_epi32(input, 6);
        case 7: return _mm_srli_epi32(input, 7);
        default: return _mm_setzero_si128(); // for shift >= 8
    }
}

static __m128i shift_right_epi64(const __m128i& input, const int shift) {
    switch (shift) {
        case 0: return input;
        case 1: return _mm_srli_epi64(input, 1);
        case 2: return _mm_srli_epi64(input, 2);
        case 3: return _mm_srli_epi64(input, 3);
        case 4: return _mm_srli_epi64(input, 4);
        case 5: return _mm_srli_epi64(input, 5);
        case 6: return _mm_srli_epi64(input, 6);
        case 7: return _mm_srli_epi64(input, 7);
        default: return _mm_setzero_si128(); // for shift >= 8
    }
}

static void psum_epi8_ovf(const __m128i& input, const uint32_t& idx, uint32_t& pf_sum, uint32_t& idx_run) {

    __m128i halves[2];

    const __m128i idx_vec = _mm_set1_epi16(idx);

    halves[0] = _mm_shuffle_epi8(input, _mm_set_epi8(-1,7, -1,6, -1,5, -1,4,
                                                       -1,3, -1,2, -1,1, -1,0));
    halves[1] = _mm_shuffle_epi8(input, _mm_set_epi8(-1,15,-1,14, -1,13, -1,12,
                                                        -1,11, -1,10, -1,9, -1,8));

    halves[0] = _mm_add_epi16(halves[0], _mm_slli_si128(halves[0], 2));
    halves[0] = _mm_add_epi16(halves[0], _mm_slli_si128(halves[0], 4));
    halves[0] = _mm_add_epi16(halves[0], _mm_slli_si128(halves[0], 8));

    //print16x8(halves[0]);

    __m128i idx_mask = _mm_cmple_epu16(halves[0], idx_vec);//mask for <=idx
    uint64_t lt_low = (uint32_t) _mm_movemask_epi8(idx_mask);
    //print16x8(idx_mask);

    halves[1] = _mm_add_epi16(halves[1], _mm_shuffle_epi8(halves[0],
                                                          _mm_set_epi8(-1,-1, -1,-1, -1,-1, -1,-1,
                                                                       -1,-1, -1,-1, -1,-1, 15,14)));
    halves[1] = _mm_add_epi16(halves[1], _mm_slli_si128(halves[1], 2));
    halves[1] = _mm_add_epi16(halves[1], _mm_slli_si128(halves[1], 4));
    halves[1] = _mm_add_epi16(halves[1], _mm_slli_si128(halves[1], 8));

    //print16x8(halves[1]);

    idx_mask = _mm_cmple_epu16(halves[1], idx_vec);//mask for <=idx
    //print16x8(idx_mask);
    const uint64_t lt_high = (uint32_t) _mm_movemask_epi8(idx_mask);

    const uint64_t lt = lt_low | (lt_high<<16);//two bits per element
    idx_run = __builtin_ctzll(~lt)>>1;

    uint16_t pf_sum_vec[8]={0};
    _mm_storeu_si128((__m128i*)pf_sum_vec, halves[idx_run>>3]);
    pf_sum = pf_sum_vec[idx_run & 7];
}

static uint32_t hsum_epi8_ovf(const __m128i& input) {
    const __m128i sum1 = _mm_add_epi16(_mm_shuffle_epi8(input, _mm_set_epi8(-1,7,  -1,6,  -1,5,  -1,4,  -1,3,  -1,2,  -1,1, -1,0)),
                                       _mm_shuffle_epi8(input, _mm_set_epi8(-1,15, -1,14, -1,13, -1,12, -1,11, -1,10, -1,9, -1,8)));
    //print16x8(sum1);
    const __m128i sum2 = _mm_add_epi16(sum1, _mm_srli_si128(sum1, 2));
    //print16x8(sum2);
    const __m128i sum3 = _mm_add_epi16(sum2, _mm_srli_si128(sum2, 4));
    //print16x8(sum3);
    const __m128i sum4 = _mm_add_epi16(sum3, _mm_srli_si128(sum3, 8));
    //print16x8(sum4);
    return static_cast<uint16_t>(_mm_cvtsi128_si32(sum4));//extract the lowest 32 bits
}

static uint32_t hsum_epi8(const __m128i& input) {
    __m128i sad = _mm_sad_epu8(input, _mm_setzero_si128());
    sad = _mm_add_epi64(sad, _mm_srli_si128(sad, 8));
    return _mm_cvtsi128_si32(sad);
}

static void psum_epi16_ovf(const __m128i& input, const uint32_t& idx, uint32_t& pf_sum, uint32_t& idx_run) {

    __m128i halves[2];

    halves[0] = _mm_shuffle_epi8(input, _mm_set_epi8(-1,-1,7,6, -1,-1,5,4, -1,-1,3,2, -1,-1,1,0));
    halves[1] = _mm_shuffle_epi8(input, _mm_set_epi8(-1,-1,15,14, -1,-1,13,12, -1,-1,11,10, -1,-1,9,8));

    halves[0] = _mm_add_epi32(halves[0], _mm_slli_si128(halves[0], 4));
    halves[0] = _mm_add_epi32(halves[0], _mm_slli_si128(halves[0], 8));

    //print32x4(halves[0]);

    const __m128i idx_vec =  _mm_set1_epi32(idx);
    __m128i idx_mask = _mm_cmple_epu32(halves[0], idx_vec);//mask for >idx
    //uint64_t less_than = (uint32_t)_mm_cvtsi128_si32(_mm_shuffle_epi8(idx_mask,
    //                                                                  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1,
    //                                                                               -1,-1,-1,-1, 12,8,4,0)));
    //print32x4(idx_mask);
    uint64_t lt_low = (uint32_t)_mm_movemask_epi8(idx_mask);//4 bits represent one element

    halves[1] = _mm_add_epi32(halves[1], _mm_shuffle_epi8(halves[0], _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1,
                                                                                  -1,-1,-1,-1, 15,14,13,12)));
    halves[1] = _mm_add_epi32(halves[1], _mm_slli_si128(halves[1], 4));
    halves[1] = _mm_add_epi32(halves[1], _mm_slli_si128(halves[1], 8));

    //print32x4(halves[1]);
    //print32x4(idx_vec);

    idx_mask = _mm_cmple_epu32(halves[1], idx_vec);//4 bits represent one element
    //uint64_t lt = (uint32_t)_mm_cvtsi128_si32(_mm_shuffle_epi8(idx_mask, _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1,
    //                                                                                  -1,-1,-1,-1, 12,8,4,0)));
    const uint64_t lt_high = (uint32_t)_mm_movemask_epi8(idx_mask);//4 bits represent one element

    lt_low |= lt_high<<16;
    //print32x4(idx_mask);
    idx_run = __builtin_ctzll(~lt_low)>>2;//

    uint32_t pf_sum_vec[8];
    _mm_storeu_si128((__m128i*)pf_sum_vec, halves[idx_run>>2]);
    pf_sum = pf_sum_vec[idx_run & 3];
}

static uint32_t hsum_epi16_ovf(const __m128i& input) {
    const __m128i sum = _mm_add_epi32(_mm_shuffle_epi8(input, _mm_set_epi8(-1,-1,7,6, -1,-1,5,4, -1,-1,3,2, -1,-1,1,0)),
                                      _mm_shuffle_epi8(input, _mm_set_epi8(-1,-1,15,14, -1,-1,13,12, -1,-1,11,10, -1,-1,9,8)));
    const __m128i sum2 = _mm_hadd_epi32(sum, sum);
    uint64_t res = _mm_cvtsi128_si64(sum2);
    res = (res & 0xFFFFFFFF) + (res>>32);
    return static_cast<uint32_t>(res);
}

static uint32_t hsum_epi16(const __m128i& input) {
    const __m128i sum1 = _mm_add_epi16(input, _mm_srli_si128(input, 2));
    const __m128i sum2 = _mm_add_epi16(sum1, _mm_srli_si128(sum1, 4));
    const __m128i sum3 = _mm_add_epi16(sum2, _mm_srli_si128(sum2, 8));
    return static_cast<uint16_t>(_mm_cvtsi128_si32(sum3));
}

static uint32_t hsum_epi32(const __m128i& input) {
    const __m128i sum1 = _mm_add_epi32(input, _mm_srli_si128(input, 4));
    const __m128i sum2 = _mm_add_epi32(sum1, _mm_srli_si128(sum1, 8));
    return _mm_cvtsi128_si32(sum2);
}

static uint64_t hsum_sum_epi64(const __m128i& input) {
    const __m128i sum = _mm_add_epi32(input, _mm_srli_si128(input, 8));
    return _mm_cvtsi128_si64(sum);
}

template<bool vbyte_compressed, uint8_t ctr_width, uint8_t bytes_per_run>
static __m128i decode_block_sse42(const uint8_t **stream) {

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
        const __m128i compressed =  _mm_loadu_si128((const __m128i*)(*stream+1));
        const __m128i dec_shuffle = _mm_loadu_si128((const __m128i*)pshuf);
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

template<bool overflow16, bool overflow32=false, bool get_run_id=false>
static auto inv_select_sse42_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx) {

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    __m128i block = _mm_loadu_si128((const __m128i*)*stream);
    *stream+=16;
    __m128i bk_lengths =  shift_right_epi8(block, sigma_bits);
    uint32_t prev_acc = 0;

    //print8x16(bk_lengths);

    //sometimes the back of the block has garbage, so I have to assume overflow
    uint32_t acc = hsum_epi8_ovf(bk_lengths);

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

    if constexpr (get_run_id){
        const auto run_id= l*16 + idx_run;
        int64_t options[2] = {-1, static_cast<int64_t>(run_id)};
        return std::make_tuple(rank, sym, options[idx==pf_sum]);
    } else {
        return std::make_pair(rank, sym);
    }
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false, bool get_run_id=false>
static auto inv_select_sse42_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx) {

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    __m128i block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    __m128i bk_lengths =  shift_right_epi16(block, sigma_bits);
    uint32_t prev_acc = 0;

    //some time the block has some garbage, so we have to assume overflow at the end
    uint32_t acc = hsum_epi16_ovf(bk_lengths);

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

    if constexpr (get_run_id){
        const auto run_id= l*8 + idx_run;
        int64_t options[2] = {-1, static_cast<int64_t>(run_id)};
        return std::make_tuple(rank, sym, options[idx==pf_sum]);
    }else{
        return std::make_pair(rank, sym);
    }

}

template<bool vbyte_compressed, uint8_t bytes_per_run, bool get_run_id=false>
static auto inv_select_sse42_32x4(const uint8_t ** stream, uint8_t sigma, uint64_t idx) {

    //TODO: assert idx fits 2 bytes
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    __m128i block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    __m128i bk_lengths =  shift_right_epi32(block, sigma_bits);
    uint32_t prev_acc=0;

    uint32_t acc = hsum_epi32(bk_lengths);
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
    rank += idx-pf_sum;

    if constexpr (get_run_id){
        const auto run_id= l*4 + idx_run;
        int64_t options[2] = {-1, static_cast<int64_t>(run_id)};
        return std::make_tuple(rank, sym, options[idx==pf_sum]);
    }else{
        return std::make_pair(rank, sym);
    }
}

//byte compressed by default
template<uint8_t bytes_per_run, bool get_run_id=false>
static auto inv_select_sse42_64x2(const uint8_t ** stream, uint8_t sigma, uint64_t idx){
    if constexpr (get_run_id){
        return std::make_tuple<int64_t, uint8_t, uint64_t>(0,0, 0);
    }else{
        return std::make_pair<int64_t, uint8_t>(0,0);
    }
}

template<bool overflow16, bool overflow32=false>
static uint8_t access_sse42_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    __m128i block = _mm_loadu_si128((const __m128i*)*stream);
    *stream+=16;
    __m128i bk_lengths =  shift_right_epi8(block, sigma_bits);

    //print8x16(bk_lengths);

    uint32_t prev_acc = 0;
    //sometimes the back of the block has garbage, so I have to assume overflow
    uint32_t acc = hsum_epi8_ovf(bk_lengths);

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
    uint32_t idx_run;
    uint8_t alpha_m = (1UL << sigma_bits)-1, run;

    if constexpr(overflow16){
        uint32_t pf_sum;
        psum_epi8_ovf(bk_lengths, idx, pf_sum, idx_run);
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
    }
    run = stream_start[(l*16)+idx_run];
    return run & alpha_m;
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static uint8_t access_sse42_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    const uint8_t sigma_bits = sym_width(sigma);

    __m128i block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    __m128i bk_lengths =  shift_right_epi16(block, sigma_bits);

    uint32_t prev_acc = 0;
    //some time the block has some garbage, so we have to assume overflow at the end
    uint32_t acc = hsum_epi16_ovf(bk_lengths);

    while(acc<=idx){
        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        bk_lengths =  shift_right_epi16(block, sigma_bits);
        //print16x8(bk_lengths);

        prev_acc = acc;
        acc += hsum_epi16_ovf(bk_lengths);
    }

    idx-=prev_acc;
    uint32_t idx_run;

    uint8_t alpha_m = (1UL << sigma_bits)-1;

    if constexpr (overflow8){
        uint32_t pf_sum;
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
    auto run = (uint16_t)_mm_cvtsi128_si32(run_vec);
    return run & alpha_m;
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static uint8_t access_sse42_32x4(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    const uint8_t sigma_bits = sym_width(sigma);

    __m128i block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    __m128i bk_lengths =  shift_right_epi32(block, sigma_bits);

    uint32_t prev_acc=0;
    uint32_t acc= hsum_epi32(bk_lengths);
    //print32x4(bk_lengths);

    while(acc<=idx){
        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        bk_lengths =  shift_right_epi32(block, sigma_bits);

        //print32x4(bk_lengths);

        prev_acc = acc;
        acc += hsum_epi32(bk_lengths);
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

    const __m128i shuff_idxs = _mm_set_epi8(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
    const __m128i shuff = _mm_add_epi8(shuff_idxs, _mm_set1_epi8(idx_run<<2));

    //print8x16(shuff);

    const __m128i run_vec = _mm_shuffle_epi8(block, shuff);
    const uint32_t run = _mm_cvtsi128_si32(run_vec);

    return run & alpha_m;
}

template<uint8_t bytes_per_run>
static uint8_t access_sse42_64x2(const uint8_t **stream, uint8_t sigma, uint64_t idx){
    return 0;
}

template<bool overflow16, bool overflow32=false>
static std::pair<uint64_t, uint64_t> get_phi_run_sse42_8x16(const uint8_t **stream, uint64_t idx){

    //NOTE here I do not need to vbyte compress the block
    __m128i block = _mm_loadu_si128((const __m128i*)*stream);
    *stream+=16;

    uint32_t idx_run = 0;

    //print8x16(bk_lengths);
    //sometimes the back of the block has garbage, so I have to assume overflow
    uint32_t prev_acc = 0;
    uint32_t acc = hsum_epi8_ovf(block);

    while(acc<=idx){
        block = _mm_loadu_si128((const __m128i*)*stream);
        *stream+=16;
        idx_run+=16;

        //print8x16(bk_lengths);
        prev_acc = acc;
        acc += hsum_epi8_ovf(block);
    }

    idx-=prev_acc;

    uint32_t pf_sum, tmp_idx_run;
    if constexpr(overflow16){
        psum_epi8_ovf(block, idx, pf_sum, tmp_idx_run);
    }else{
        //prefix sum without overflow
        block = _mm_add_epi8(block, _mm_slli_si128(block, 1));
        block = _mm_add_epi8(block, _mm_slli_si128(block, 2));
        block = _mm_add_epi8(block, _mm_slli_si128(block, 4));
        block = _mm_add_epi8(block, _mm_slli_si128(block, 8));
        //print8x16(bk_lengths);
        const __m128i idx_mask = _mm_cmple_epu8(block, _mm_set1_epi8(idx));//mask for >idx
        const uint16_t less_than = _mm_movemask_epi8(idx_mask);
        tmp_idx_run = __builtin_ctzll(~less_than);

        const __m128i shuff = _mm_set1_epi8(idx_run);
        const __m128i pf_sum_vec = _mm_shuffle_epi8(block, shuff);
        pf_sum = ((uint8_t*)&pf_sum_vec)[0];
    }

    const uint32_t offset = static_cast<uint32_t>((*stream - 16)[tmp_idx_run]) - (pf_sum-idx);
    return std::make_pair(idx_run+tmp_idx_run, offset);
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static std::pair<uint64_t, uint64_t> get_phi_run_sse42_16x8(const uint8_t **stream, uint64_t idx){

    __m128i block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    uint32_t prev_acc = 0;
    //some time the block has some garbage, so we have to assume overflow at the end
    uint32_t acc = hsum_epi16_ovf(block);
    uint32_t idx_run =0;

    while(acc<=idx){
        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        idx_run+=8;
        //print16x8(bk_lengths);
        prev_acc = acc;
        acc += hsum_epi16_ovf(block);
    }

    idx-=prev_acc;
    uint32_t pf_sum, tmp_idx_run, len;

    if constexpr (overflow8){
        psum_epi16_ovf(block, idx, pf_sum, tmp_idx_run);
        const __m128i shuff = _mm_add_epi8(_mm_set_epi64x(0x100010001000100ULL, 0x100010001000100ULL),
                                           _mm_set1_epi8(tmp_idx_run<<1));
        len = (uint16_t)_mm_cvtsi128_si32(_mm_shuffle_epi8(block, shuff));
    }else{
        //prefix sum without overflow
        __m128i pf_sum_vec = _mm_add_epi16(block, _mm_slli_si128(block, 2));
        pf_sum_vec = _mm_add_epi16(pf_sum_vec, _mm_slli_si128(pf_sum_vec, 4));
        pf_sum_vec = _mm_add_epi16(pf_sum_vec, _mm_slli_si128(pf_sum_vec, 8));
        //print16x8(bk_lengths);

        const __m128i idx_mask = _mm_cmple_epu16(pf_sum_vec, _mm_set1_epi16(idx));//mask for <=idx

        //const auto less_than = (uint8_t)_mm_movemask_epi8(_mm_packs_epi16(idx_mask, _mm_setzero_si128()));
        //tmp_idx_run = __builtin_ctzll(~less_than);

        const int less_than = _mm_movemask_epi8(idx_mask);
        tmp_idx_run = __builtin_ctz(~less_than)>>1;

        //_mm_set_epi64x(0x100010001000100ULL, 0x100010001000100ULL) is equal to set {0, 1, 0, 1, 0, 1, ...}
        const __m128i shuff = _mm_add_epi8(_mm_set_epi64x(0x100010001000100ULL, 0x100010001000100ULL),
                                           _mm_set1_epi8(tmp_idx_run<<1));
        pf_sum = (uint16_t)_mm_cvtsi128_si32(_mm_shuffle_epi8(pf_sum_vec, shuff));
        len = (uint16_t)_mm_cvtsi128_si32(_mm_shuffle_epi8(block, shuff));
    }

    //print8x16(shuff);
    uint64_t offset = len - (pf_sum-idx);
    return std::make_pair(idx_run+tmp_idx_run, offset);
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static std::pair<uint64_t, uint64_t> get_phi_run_sse42_32x4(const uint8_t **stream, uint64_t idx){

    __m128i block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    uint32_t prev_acc=0;
    uint32_t acc= hsum_epi32(block);
    uint32_t idx_run =0;
    //print32x4(bk_lengths);

    while(acc<=idx){
        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        prev_acc = acc;
        idx_run+=4;
        acc += hsum_epi32(block);
    }

    idx-=prev_acc;

    //prefix sum
    __m128i pf_sum_vec = _mm_add_epi32(block, _mm_slli_si128(block, 4));
    pf_sum_vec = _mm_add_epi32(pf_sum_vec, _mm_slli_si128(pf_sum_vec, 8));
    //

    //print32x4(bk_lengths);

    const __m128i idx_mask = _mm_cmple_epu32(pf_sum_vec, _mm_set1_epi32(idx));//mask for >idx
    const uint64_t less_than = _mm_cvtsi128_si64(_mm_packs_epi32(idx_mask, _mm_setzero_si128()));
    const uint8_t tmp_idx_run = __builtin_ctzll(~less_than)>>4;

    //_mm_set_epi64x(0x302010003020100ULL, 0x302010003020100uLL) equal to set {0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4}
    const __m128i shuff = _mm_add_epi8(_mm_set_epi64x(0x302010003020100ULL, 0x302010003020100ULL),
                                       _mm_set1_epi8(tmp_idx_run<<2));
    const uint32_t pf_sum = _mm_cvtsi128_si32(_mm_shuffle_epi8(pf_sum_vec, shuff));
    const uint32_t len = _mm_cvtsi128_si32(_mm_shuffle_epi8(block, shuff));

    const uint32_t offset = len - static_cast<uint32_t>(pf_sum - idx);
    return std::make_pair(idx_run+tmp_idx_run, offset);
}

template<uint8_t bytes_per_run>
static std::pair<uint64_t, uint64_t> get_phi_run_sse42_64x2(const uint8_t **stream, uint64_t idx){
    return std::make_pair(0, 0);
}

template<bool overflow16, bool overflow32=false, bool check_head>
static int64_t rank_sse42_8x16(const uint8_t **stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi8(alpha_m);
    const __m128i sym_vec = _mm_set1_epi8(symbol);

    __m128i block = _mm_loadu_si128((const __m128i*)*stream);
    *stream+=16;
    __m128i bk_lengths =  shift_right_epi8(block, sigma_bits);

    uint32_t prev_acc = 0;
    //the back of the block *might* contain garbage, so I have to assume overflow
    uint32_t acc = hsum_epi8_ovf(bk_lengths);

    uint32_t rank=0;
    //size_t l=0;
    while(acc<=idx){
        //compute acc rank in the previous block
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank += hsum_epi8_ovf(bk_lengths);
        } else {
            rank += hsum_epi8(bk_lengths);
        }

        block = _mm_loadu_si128((const __m128i*)*stream);
        *stream+=16;
        bk_lengths =  shift_right_epi8(block, sigma_bits);
        //print8x16(bk_lengths);
        prev_acc = acc;
        acc += hsum_epi8_ovf(bk_lengths);
        //l++;
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr(overflow16){
        psum_epi8_ovf(bk_lengths, idx, pf_sum, idx_run);
    } else{
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
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
        pf_sum = ((uint8_t*)&pf_sum_vec)[0];
    }

    const uint8_t run = (*stream-16)[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    //create a mask for the index position
    //eg: idx_run = 3 yields {0:0xFF, 1:0xFF, 2:0xFF, 3:0xFF, 4:0, 5:0, ...}
    //static const __m128i indices = _mm_set_epi8(15,14,13,12,11,10,9,8, 7,6,5,4,3,2,1,0);
    static const __m128i indices = _mm_set_epi64x(0xF0E0D0C0B0A0908,0x706050403020100ULL);
    const __m128i idx_run_mask = _mm_cmpgt_epi8( _mm_set1_epi8(idx_run+1), indices);

    bk_lengths =  shift_right_epi8(block, sigma_bits);
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));

    if constexpr (overflow16){
        rank += hsum_epi8_ovf(bk_lengths);
    }else{
        rank += hsum_epi8(bk_lengths);
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
static int64_t rank_sse42_16x8(const uint8_t **stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi16(alpha_m);
    const __m128i sym_vec = _mm_set1_epi16(symbol);

    __m128i block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    __m128i bk_lengths =  shift_right_epi16(block, sigma_bits);

    //the last block might have some garbage, so we have to assume overflow at the end
    uint32_t prev_acc = 0;
    uint32_t acc = hsum_epi16_ovf(bk_lengths);
    uint64_t rank = 0;

    while(acc<=idx){
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank += hsum_epi16_ovf(bk_lengths);
        }else{
            rank += hsum_epi16(bk_lengths);
        }

        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        bk_lengths =  shift_right_epi16(block, sigma_bits);

        prev_acc = acc;
        acc += hsum_epi16_ovf(bk_lengths);
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx, pf_sum, idx_run);
    } else {
        //prefix sum without overflow
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 2));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 4));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 8));

        const __m128i idx_mask = _mm_cmple_epu16(bk_lengths, _mm_set1_epi16(idx));//mask for <=idx
        const auto less_than = (uint8_t)_mm_movemask_epi8(_mm_packs_epi16(idx_mask, _mm_setzero_si128()));
        idx_run = __builtin_ctzll(~less_than);
    }

    //_mm_set_epi64x(0x100010001000100ULL, 0x100010001000100ULL) is equal to set {0, 1, 0, 1, 0, 1, ...}
    const __m128i shuff = _mm_add_epi8(_mm_set_epi64x(0x100010001000100ULL, 0x100010001000100ULL),
                                       _mm_set1_epi8(idx_run<<1));

    const __m128i run_vec = _mm_shuffle_epi8(block, shuff);

    const auto run = (uint16_t)_mm_cvtsi128_si32(run_vec);
    const uint8_t last_symbol = run & alpha_m;

    if constexpr (!overflow8){
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
        pf_sum = (uint16_t)_mm_cvtsi128_si32(pf_sum_vec);
    }

    //create a mask for the index position
    //eg: idx_run = 3 yields {0:0xFFFF, 1:0xFFFF, 2:0xFFFF, 3:0xFFFF, 4:0, 5:0, ...}
    //static const __m128i indices = _mm_set_epi16(7,6,5,4,3,2,1,0);
    static const __m128i indices = _mm_set_epi64x(0x7000600050004ULL, 0x3000200010000ULL);
    const __m128i idx_run_mask = _mm_cmpgt_epi16( _mm_set1_epi16(idx_run+1), indices);

    bk_lengths =  shift_right_epi16(block, sigma_bits);
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));

    if constexpr (overflow8){
        rank += hsum_epi16_ovf(bk_lengths);
    }else{
        rank += hsum_epi16(bk_lengths);
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
static int64_t rank_sse42_32x4(const uint8_t ** stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi32(alpha_m);
    const __m128i sym_vec = _mm_set1_epi32(symbol);

    __m128i block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    __m128i bk_lengths =  shift_right_epi32(block, sigma_bits);

    uint32_t prev_acc=0;
    uint64_t acc=hsum_epi32(bk_lengths);
    uint64_t rank=0;

    while(acc<=idx){
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));
        rank += hsum_epi32(bk_lengths);

        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        bk_lengths =  shift_right_epi32(block, sigma_bits);

        prev_acc = acc;
        acc += hsum_epi32(bk_lengths);
    }

    idx-=prev_acc;

    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 4));
    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 8));

    const __m128i idx_mask = _mm_cmple_epu32(bk_lengths, _mm_set1_epi32(idx));//mask for >idx

    const uint64_t less_than = _mm_cvtsi128_si64(_mm_packs_epi32(idx_mask, _mm_setzero_si128()));
    const uint8_t idx_run = __builtin_ctzll(~less_than)>>4;

    //_mm_set_epi64x(0x302010003020100ULL, 0x302010003020100uLL) equal to set {0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4}
    const __m128i shuff = _mm_add_epi8(_mm_set_epi64x(0x302010003020100ULL, 0x302010003020100ULL),
                                       _mm_set1_epi8(idx_run<<2));

    const __m128i run_vec = _mm_shuffle_epi8(block, shuff);
    const uint32_t run = _mm_cvtsi128_si32(run_vec);
    const uint8_t last_symbol = run & alpha_m;

    const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
    const uint32_t pf_sum = _mm_cvtsi128_si32(pf_sum_vec);

    //create a mask for the index position
    //eg: idx_run = 3 yields {0:0xFFFFFFFF, 1:0xFFFFFFFF, 2:0xFFFFFFFF, 3:0xFFFFFFFF, 4:0, 5:0, ...}
    static const __m128i indices = _mm_set_epi64x(0x300000002ULL, 0x100000000ULL);
    //static const __m128i indices = _mm_set_epi32(3,2,1,0);
    const __m128i idx_run_mask = _mm_cmpgt_epi32( _mm_set1_epi32(idx_run+1), indices);

    bk_lengths =  shift_right_epi32(block, sigma_bits);
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));

    rank += hsum_epi32(bk_lengths);

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
static int64_t rank_sse42_64x2(const uint8_t ** stream, uint8_t sigma, uint64_t idx, uint8_t symbol){
    return 0;
}

template<bool overflow16, bool overflow32=false, bool check_head>
static std::pair<uint64_t, uint64_t> range_rank_sse42_8x16(const uint8_t **stream, const uint8_t sigma,
                                                           uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi8(alpha_m);
    const __m128i sym_vec = _mm_set1_epi8(symbol);

    __m128i block = _mm_loadu_si128((const __m128i*)*stream);
    *stream+=16;
    __m128i bk_lengths =  shift_right_epi8(block, sigma_bits);

    uint32_t prev_acc = 0;
    //print8x16(bk_lengths);
    //the back of the block *might* contain garbage, so I have to assume overflow
    uint32_t acc = hsum_epi8_ovf(bk_lengths);

    uint32_t rank_i=0;
    //size_t l=0;
    while(acc<=idx_i){
        //compute acc rank in the previous block
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank_i += hsum_epi8_ovf(bk_lengths);
        } else {
            rank_i += hsum_epi8(bk_lengths);
        }

        block = _mm_loadu_si128((const __m128i*)*stream);
        *stream+=16;
        bk_lengths =  shift_right_epi8(block, sigma_bits);
        //print8x16(bk_lengths);
        prev_acc = acc;
        acc += hsum_epi8_ovf(bk_lengths);
        //l++;
    }

    uint64_t rank_j= rank_i;

    idx_i-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr(overflow16){
        psum_epi8_ovf(bk_lengths, idx_i, pf_sum, idx_run);
    } else{
        //prefix sum without overflow
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 1));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 2));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 4));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 8));
        //print8x16(bk_lengths);
        const __m128i idx_mask = _mm_cmple_epu8(bk_lengths, _mm_set1_epi8(idx_i));//mask for >idx
        const uint16_t less_than = _mm_movemask_epi8(idx_mask);
        idx_run = __builtin_ctzll(~less_than);

        //print8x16(idx_mask);
        const __m128i shuff = _mm_set1_epi8(idx_run);
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
        pf_sum = ((uint8_t*)&pf_sum_vec)[0];
    }

    uint8_t run = (*stream-16)[idx_run];
    uint8_t last_symbol = run & alpha_m;

    //create a mask for the index position
    //eg: idx_run = 3 yields {0:0xFF, 1:0xFF, 2:0xFF, 3:0xFF, 4:0, 5:0, ...}
    //static const __m128i indices = _mm_set_epi8(15,14,13,12,11,10,9,8, 7,6,5,4,3,2,1,0);
    static const __m128i indices = _mm_set_epi64x(0xF0E0D0C0B0A0908,0x706050403020100ULL);
    __m128i idx_run_mask = _mm_cmpgt_epi8( _mm_set1_epi8(idx_run+1), indices);

    bk_lengths =  shift_right_epi8(block, sigma_bits);
    //bk_lengths = _mm_and_si128(bk_lengths, _mm_loadu_si128((const __m128i*)mask8x16[idx_run+1]));
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));

    if constexpr (overflow16){
        rank_i += hsum_epi8_ovf(bk_lengths);
    }else{
        rank_i += hsum_epi8(bk_lengths);
    }

    if constexpr (check_head) {
        const uint8_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank_i -=(pf_sum-idx_i) * is_same_sym;
        rank_i = (rank_i<<1) | is_head;
        rank_i = (rank_i<<1) | is_same_sym;
    } else {
        rank_i -=(pf_sum-idx_i) * (last_symbol==symbol);
    }

    //========
    //deal with j
    bk_lengths =  shift_right_epi8(block, sigma_bits);
    //acc and the block are still valid, no need to recompute them
    while(acc<=idx_j){
        //compute acc rank in the previous block
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank_j += hsum_epi8_ovf(bk_lengths);
        } else {
            rank_j += hsum_epi8(bk_lengths);
        }

        block = _mm_loadu_si128((const __m128i*)*stream);
        *stream+=16;
        bk_lengths =  shift_right_epi8(block, sigma_bits);
        //print8x16(bk_lengths);
        prev_acc = acc;
        acc += hsum_epi8_ovf(bk_lengths);
        //l++;
    }

    idx_j-=prev_acc;

    if constexpr(overflow16){
        psum_epi8_ovf(bk_lengths, idx_j, pf_sum, idx_run);
    } else{
        //prefix sum without overflow
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 1));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 2));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 4));
        bk_lengths = _mm_add_epi8(bk_lengths, _mm_slli_si128(bk_lengths, 8));
        //print8x16(bk_lengths);
        const __m128i idx_mask = _mm_cmple_epu8(bk_lengths, _mm_set1_epi8(idx_j));//mask for >idx
        const uint16_t less_than = _mm_movemask_epi8(idx_mask);
        idx_run = __builtin_ctzll(~less_than);

        //print8x16(idx_mask);
        const __m128i shuff = _mm_set1_epi8(idx_run);
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
        pf_sum = ((uint8_t*)&pf_sum_vec)[0];
    }

    run = (*stream-16)[idx_run];
    last_symbol = run & alpha_m;

    idx_run_mask = _mm_cmpgt_epi8( _mm_set1_epi8(idx_run+1), indices);

    bk_lengths =  shift_right_epi8(block, sigma_bits);
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));

    if constexpr (overflow16){
        rank_j += hsum_epi8_ovf(bk_lengths);
    }else{
        rank_j += hsum_epi8(bk_lengths);
    }
    rank_j -=(pf_sum-idx_j) * (last_symbol==symbol);

    return std::make_pair(rank_i, rank_j);
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false, bool check_head>
static std::pair<uint64_t, uint64_t> range_rank_sse42_16x8(const uint8_t **stream, const uint8_t sigma,
                                                           uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi16(alpha_m);
    const __m128i sym_vec = _mm_set1_epi16(symbol);

    __m128i block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    __m128i bk_lengths =  shift_right_epi16(block, sigma_bits);

    //the last block might have some garbage, so we have to assume overflow at the end
    uint32_t prev_acc = 0;
    uint32_t acc = hsum_epi16_ovf(bk_lengths);
    uint64_t rank_i = 0;

    while(acc<=idx_i){
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank_i += hsum_epi16_ovf(bk_lengths);
        }else{
            rank_i += hsum_epi16(bk_lengths);
        }

        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        //print16x8(block);
        bk_lengths =  shift_right_epi16(block, sigma_bits);
        //print16x8(bk_lengths);

        prev_acc = acc;
        acc += hsum_epi16_ovf(bk_lengths);
    }

    uint64_t rank_j = rank_i;

    idx_i-=prev_acc;

    uint32_t idx_run, pf_sum;
    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx_i, pf_sum, idx_run);
    }else{
        //prefix sum without overflow
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 2));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 4));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 8));
        //print16x8(bk_lengths);

        const __m128i idx_mask = _mm_cmple_epu16(bk_lengths, _mm_set1_epi16(idx_i));//mask for <=idx
        const auto less_than = (uint8_t)_mm_movemask_epi8(_mm_packs_epi16(idx_mask, _mm_setzero_si128()));
        idx_run = __builtin_ctzll(~less_than);
        //print16x8(idx_mask);
    }

    //_mm_set_epi64x(0x100010001000100ULL, 0x100010001000100ULL) is equal to set {0, 1, 0, 1, 0, 1, ...}
    const __m128i shuff_i = _mm_add_epi8(_mm_set_epi64x(0x100010001000100ULL, 0x100010001000100ULL),
                                       _mm_set1_epi8(idx_run<<1));

    __m128i run_vec = _mm_shuffle_epi8(block, shuff_i);
    auto run = (uint16_t)_mm_cvtsi128_si32(run_vec);
    uint8_t last_symbol = run & alpha_m;

    if constexpr (!overflow8){
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff_i);
        pf_sum = (uint16_t)_mm_cvtsi128_si32(pf_sum_vec);
    }

    //create a mask for the index position
    //eg: idx_run = 3 yields {0:0xFFFF, 1:0xFFFF, 2:0xFFFF, 3:0xFFFF, 4:0, 5:0, ...}
    //static const __m128i indices = _mm_set_epi16(7,6,5,4,3,2,1,0);
    static const __m128i indices = _mm_set_epi64x(0x7000600050004ULL, 0x3000200010000ULL);
    __m128i idx_run_mask = _mm_cmpgt_epi16( _mm_set1_epi16(idx_run+1), indices);

    bk_lengths =  shift_right_epi16(block, sigma_bits);
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));

    if constexpr (overflow8){
        rank_i += hsum_epi16_ovf(bk_lengths);
    }else{
        rank_i += hsum_epi16(bk_lengths);
    }

    if constexpr (check_head) {
        const uint16_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank_i -= (pf_sum-idx_i) * is_same_sym;
        rank_i = (rank_i<<1) | is_head;
        rank_i = (rank_i<<1) | is_same_sym;
    } else {
        rank_i -= (pf_sum-idx_i) * (last_symbol==symbol);
    }

    //======
    //deal with j
    bk_lengths =  shift_right_epi16(block, sigma_bits);
    //acc and the block are still valid, no need to recompute them
    while(acc<=idx_j){
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank_j += hsum_epi16_ovf(bk_lengths);
        }else{
            rank_j += hsum_epi16(bk_lengths);
        }

        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        bk_lengths =  shift_right_epi16(block, sigma_bits);

        prev_acc = acc;
        acc += hsum_epi16_ovf(bk_lengths);
    }

    idx_j-=prev_acc;
    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx_j, pf_sum, idx_run);
    }else{
        //prefix sum without overflow
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 2));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 4));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 8));

        const __m128i idx_mask = _mm_cmple_epu16(bk_lengths, _mm_set1_epi16(idx_j));//mask for <=idx
        const auto less_than = (uint8_t)_mm_movemask_epi8(_mm_packs_epi16(idx_mask, _mm_setzero_si128()));
        idx_run = __builtin_ctzll(~less_than);
    }

    const __m128i shuff_j = _mm_add_epi8(_mm_set_epi64x(0x100010001000100ULL, 0x100010001000100ULL),
                                   _mm_set1_epi8(idx_run<<1));
    run_vec = _mm_shuffle_epi8(block, shuff_j);
    run = (uint16_t)_mm_cvtsi128_si32(run_vec);
    last_symbol = run & alpha_m;

    if constexpr (!overflow8){
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff_j);
        pf_sum = (uint16_t)_mm_cvtsi128_si32(pf_sum_vec);
    }

    idx_run_mask = _mm_cmpgt_epi16( _mm_set1_epi16(idx_run+1), indices);
    bk_lengths =  shift_right_epi16(block, sigma_bits);
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));

    if constexpr (overflow8){
        rank_j += hsum_epi16_ovf(bk_lengths);
    }else{
        rank_j += hsum_epi16(bk_lengths);
    }
    rank_j -= (pf_sum-idx_j) * (last_symbol==symbol);

    return std::make_pair(rank_i, rank_j);
}

template<bool vbyte_compressed, uint8_t bytes_per_run, bool check_head>
static std::pair<uint64_t, uint64_t> range_rank_sse42_32x4(const uint8_t ** stream, const uint8_t sigma,
                                                           uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi32(alpha_m);
    const __m128i sym_vec = _mm_set1_epi32(symbol);

    __m128i block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    __m128i bk_lengths =  shift_right_epi32(block, sigma_bits);

    uint32_t prev_acc=0;
    uint64_t acc=hsum_epi32(bk_lengths);
    uint64_t rank_i=0;

    while(acc<=idx_i){
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));
        rank_i += hsum_epi32(bk_lengths);

        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        bk_lengths =  shift_right_epi32(block, sigma_bits);

        prev_acc = acc;
        acc += hsum_epi32(bk_lengths);
    }

    uint64_t rank_j=rank_i;

    idx_i-=prev_acc;
    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 4));
    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 8));

    __m128i idx_mask = _mm_cmple_epu32(bk_lengths, _mm_set1_epi32(idx_i));//mask for >idx

    int64_t less_than = _mm_cvtsi128_si64(_mm_packs_epi32(idx_mask, _mm_setzero_si128()));
    uint8_t idx_run = __builtin_ctzll(~less_than)>>4;

    //_mm_set_epi64x(0x302010003020100ULL, 0x302010003020100uLL) equal to set {0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4}
    const __m128i shuff_i = _mm_add_epi8(_mm_set_epi64x(0x302010003020100ULL, 0x302010003020100ULL),
                                       _mm_set1_epi8(idx_run<<2));

    __m128i run_vec = _mm_shuffle_epi8(block, shuff_i);
    int32_t run = _mm_cvtsi128_si32(run_vec);
    uint8_t last_symbol = run & alpha_m;

    __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff_i);
    uint32_t pf_sum = _mm_cvtsi128_si32(pf_sum_vec);

    //create a mask for the index position
    //eg: idx_run = 3 yields {0:0xFFFFFFFF, 1:0xFFFFFFFF, 2:0xFFFFFFFF, 3:0xFFFFFFFF, 4:0, 5:0, ...}
    static const __m128i indices = _mm_set_epi64x(0x300000002ULL, 0x100000000ULL);
    __m128i idx_run_mask = _mm_cmpgt_epi32( _mm_set1_epi32(idx_run+1), indices);

    bk_lengths =  shift_right_epi32(block, sigma_bits);
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));

    rank_i += hsum_epi32(bk_lengths);

    if constexpr (check_head) {
        const uint32_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank_i -= (pf_sum-idx_i) * is_same_sym;
        rank_i = (rank_i<<1) | is_head;
        rank_i = (rank_i<<1) | is_same_sym;
    } else {
        rank_i -= (pf_sum-idx_i)*(last_symbol==symbol);
    }

    //=== deal with j
    bk_lengths =  shift_right_epi32(block, sigma_bits);
    //acc and the block are still valid, no need to recompute them
    while(acc<=idx_j){
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));
        rank_j += hsum_epi32(bk_lengths);

        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        bk_lengths =  shift_right_epi32(block, sigma_bits);

        prev_acc = acc;
        acc += hsum_epi32(bk_lengths);
    }

    idx_j-=prev_acc;
    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 4));
    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 8));

    idx_mask = _mm_cmple_epu32(bk_lengths, _mm_set1_epi32(idx_j));//mask for >idx
    less_than = _mm_cvtsi128_si64(_mm_packs_epi32(idx_mask, _mm_setzero_si128()));
    idx_run = __builtin_ctzll(~less_than)>>4;

    const __m128i shuff_j = _mm_add_epi8(_mm_set_epi64x(0x302010003020100ULL, 0x302010003020100ULL),
                                       _mm_set1_epi8(idx_run<<2));
    run_vec = _mm_shuffle_epi8(block, shuff_j);
    run = _mm_cvtsi128_si32(run_vec);
    last_symbol = run & alpha_m;

    pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff_j);
    pf_sum = _mm_cvtsi128_si32(pf_sum_vec);

    idx_run_mask = _mm_cmpgt_epi32( _mm_set1_epi32(idx_run+1), indices);

    bk_lengths =  shift_right_epi32(block, sigma_bits);
    bk_lengths = _mm_and_si128(bk_lengths, idx_run_mask);
    bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));

    rank_j += hsum_epi32(bk_lengths);
    rank_j -= (pf_sum-idx_j)*(last_symbol==symbol);

    return std::make_pair(rank_i, rank_j);
}

template<uint8_t bytes_per_run, bool check_head>
static std::pair<uint64_t, uint64_t> range_rank_sse42_64x2(const uint8_t ** stream, const uint8_t sigma, uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){
    return std::make_pair(0,0);
}

static size_t first_run_sse42_8x16(const uint8_t **stream, const uint8_t sigma, uint8_t symbol) {

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi8(alpha_m);
    const __m128i sym_vec = _mm_set1_epi8(symbol);

    __m128i block = _mm_loadu_si128((const __m128i*)*stream);
    *stream+=16;
    __m128i sym_mask = _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec);
    bool has_sym = !_mm_testz_si128(sym_mask, sym_mask);;

    size_t run_idx = 0;
    while (!has_sym){
        run_idx+=16;
        block = _mm_loadu_si128((const __m128i*)*stream);
        *stream+=16;
        sym_mask = _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec);
        has_sym = !_mm_testz_si128(sym_mask, sym_mask);
    }

    const int matches = _mm_movemask_epi8(sym_mask);
    const uint8_t first = __builtin_ctzl(matches);
    return run_idx + first;
}

template<bool vbyte_compressed>
static inline size_t first_run_sse42_16x8(const uint8_t **stream, const uint8_t sigma, uint8_t symbol) {

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi16(alpha_m);
    const __m128i sym_vec = _mm_set1_epi16(symbol);

    __m128i block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    __m128i sym_mask = _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec);
    bool has_sym = !_mm_testz_si128(sym_mask, sym_mask);;

    size_t run_idx = 0;
    while(!has_sym) {
        run_idx+=8;
        block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
        sym_mask = _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec);
        has_sym = !_mm_testz_si128(sym_mask, sym_mask);;
    }

    const auto matches = _mm_movemask_epi8(_mm_packs_epi16(sym_mask, _mm_setzero_si128()));
    const uint8_t first = __builtin_ctzl(matches);
    return run_idx + first;
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline size_t first_run_sse42_32x4(const uint8_t **stream, const uint8_t sigma, uint8_t symbol) {
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi32(alpha_m);
    const __m128i sym_vec = _mm_set1_epi32(symbol);

    __m128i block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    __m128i sym_mask = _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec);
    bool has_sym = !_mm_testz_si128(sym_mask, sym_mask);;

    size_t run_idx = 0;
    while(!has_sym) {
        run_idx+=4;
        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        sym_mask = _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec);
        has_sym = !_mm_testz_si128(sym_mask, sym_mask);;
    }

    const uint64_t matches = _mm_cvtsi128_si64(_mm_packs_epi32(sym_mask, _mm_setzero_si128()));
    const uint8_t first = __builtin_ctzll(matches)>>4;
    return run_idx + first;
}

template<uint8_t bytes_per_run>
static inline size_t first_run_sse42_64x2(const uint8_t **stream, const uint8_t sigma, uint8_t symbol){
    return 0;
}

//the function succ_neon_*** returns the run id for the leftmost run labeled "symbol" after position "idx".
//the output is a pair (run_id (int64_t), rank (uint64_t)), where "rank" is the number of times "symbol" occurs before "i"
//NOTE: run_id = -1 if there is no run labeled "symbol" from *i* onwards
//NOTE: if "idx" is the head of its run and that run is labeled "symbol", then run_id is the id for that run
template<bool overflow16, bool overflow32>
static inline std::pair<int64_t, uint64_t> succ_sse42_8x16(const uint8_t **stream, const size_t n_runs,
                                                           const uint8_t sigma, uint64_t idx, uint8_t symbol) {
    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi8(alpha_m);
    const __m128i sym_vec = _mm_set1_epi8(symbol);

    __m128i block = _mm_loadu_si128((const __m128i*)*stream);
    *stream+=16;
    __m128i bk_lengths =  shift_right_epi8(block, sigma_bits);

    uint32_t prev_acc = 0;
    //print8x16(bk_lengths);
    //the back of the block *might* contain garbage, so I have to assume overflow
    uint32_t acc = hsum_epi8_ovf(bk_lengths);

    uint32_t rank=0;
    uint64_t run_id=0;
    while(acc<=idx){
        //compute acc rank in the previous block
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank += hsum_epi8_ovf(bk_lengths);
        } else {
            rank += hsum_epi8(bk_lengths);
        }

        block = _mm_loadu_si128((const __m128i*)*stream);
        *stream+=16;
        bk_lengths =  shift_right_epi8(block, sigma_bits);
        //print8x16(bk_lengths);
        prev_acc = acc;
        acc += hsum_epi8_ovf(bk_lengths);
        run_id+=16;
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr(overflow16){
        psum_epi8_ovf(bk_lengths, idx, pf_sum, idx_run);
    } else{
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
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
        pf_sum = ((uint8_t*)&pf_sum_vec)[0];
    }

    *stream -= 16;
    const uint8_t run = (*stream)[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    block = _mm_loadu_si128((const __m128i*)*stream);
    bk_lengths =  shift_right_epi8(block, sigma_bits);
    const __m128i mask_pref = _mm_loadu_si128((const __m128i*)mask8x16[idx_run+1]);
    bk_lengths = _mm_and_si128(bk_lengths, mask_pref);
    __m128i sym_mask = _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec);
    bk_lengths = _mm_and_si128(bk_lengths, sym_mask);

    if constexpr (overflow16){
        rank += hsum_epi8_ovf(bk_lengths);
    }else{
        rank += hsum_epi8(bk_lengths);
    }

    const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
    rank -=(pf_sum-idx) * is_same_sym;

    //return if "symbol" matches the symbol where "i" lies
    int64_t opts[2] ={-1, static_cast<int64_t>(run_id)+idx_run};
    if(is_same_sym){
        bool is_head = (pf_sum-(run>>sigma_bits))==idx;//"i" is also the head position within the run
        return {opts[is_head], rank};
    }

    const __m128i mask_suff =  _mm_xor_si128(mask_pref, _mm_set1_epi32(-1));
    sym_mask = _mm_and_si128(sym_mask, mask_suff);
    bool has_sym = !_mm_testz_si128(sym_mask, sym_mask);
    run_id+=16;
    while (run_id<n_runs && !has_sym) {
        *stream+=16;
        block = _mm_loadu_si128((const __m128i*)*stream);
        sym_mask = _mm_cmpeq_epi8(_mm_and_si128(block, alpha_mask), sym_vec);
        has_sym = !_mm_testz_si128(sym_mask, sym_mask);
        run_id+=16;
    }

    const int matches = _mm_movemask_epi8(sym_mask);
    const uint8_t first = __builtin_ctzl(matches);
    run_id = run_id - 16  + first;
    opts[1] = static_cast<int64_t>(run_id);

    return {opts[has_sym && run_id<n_runs], rank};
}

template<bool vbyte_compressed, bool overflow8, bool overflow16>
static inline std::pair<uint64_t, uint64_t> succ_sse42_16x8(const uint8_t **stream, const size_t n_runs,
                                                            const uint8_t sigma, uint64_t idx, uint8_t symbol) {
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi16(alpha_m);
    const __m128i sym_vec = _mm_set1_epi16(symbol);

    const uint8_t *prev_state = *stream;
    __m128i block = decode_block_sse42<vbyte_compressed, 1, 2>(stream);
    //print16x8(block);
    __m128i bk_lengths =  shift_right_epi16(block, sigma_bits);
    //print16x8(bk_lengths);

    //the last block might have some garbage, so we have to assume overflow at the end
    uint32_t prev_acc = 0;
    uint32_t acc = hsum_epi16_ovf(bk_lengths);
    uint64_t rank = 0;
    uint64_t run_id = 0;

    while(acc<=idx){
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec));
        if constexpr (overflow16){
            rank += hsum_epi16_ovf(bk_lengths);
        }else{
            rank += hsum_epi16(bk_lengths);
        }
        run_id+=8;

        prev_state = *stream;
        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        //print16x8(block);
        bk_lengths =  shift_right_epi16(block, sigma_bits);
        //print16x8(bk_lengths);

        prev_acc = acc;
        acc += hsum_epi16_ovf(bk_lengths);
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx, pf_sum, idx_run);
    }else{
        //prefix sum without overflow
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 2));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 4));
        bk_lengths = _mm_add_epi16(bk_lengths, _mm_slli_si128(bk_lengths, 8));
        //print16x8(bk_lengths);

        const __m128i idx_mask = _mm_cmple_epu16(bk_lengths, _mm_set1_epi16(idx));//mask for <=idx
        auto less_than = (uint8_t)_mm_movemask_epi8(_mm_packs_epi16(idx_mask, _mm_setzero_si128()));
        idx_run = __builtin_ctzll(~less_than);
        //print16x8(idx_mask);
    }

    const __m128i shuff = _mm_add_epi8(_mm_set_epi8(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0),
                                       _mm_set1_epi8(idx_run<<1));
    //print8x16(shuff);

    const __m128i run_vec = _mm_shuffle_epi8(block, shuff);
    const auto run = (uint16_t)_mm_cvtsi128_si32(run_vec);
    const uint8_t last_symbol = run & alpha_m;

    if constexpr (!overflow8){
        const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
        pf_sum = (uint16_t)_mm_cvtsi128_si32(pf_sum_vec);
    }

    *stream = prev_state;
    block = decode_block_sse42<vbyte_compressed,1,2>(stream);
    bk_lengths =  shift_right_epi16(block, sigma_bits);
    const __m128i mask_pref = _mm_loadu_si128((const __m128i*)mask16x8[idx_run+1]);
    bk_lengths = _mm_and_si128(bk_lengths, mask_pref);
    __m128i sym_mask = _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec);
    bk_lengths = _mm_and_si128(bk_lengths, sym_mask);

    if constexpr (overflow8){
        rank += hsum_epi16_ovf(bk_lengths);
    }else{
        rank += hsum_epi16(bk_lengths);
    }
    const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
    rank -= (pf_sum-idx) * is_same_sym;

    int64_t opts[2] ={-1, (int64_t)run_id+idx_run};
    if(is_same_sym){
        const uint16_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_head = (pf_sum-len)==idx;//check if idx is the head of the run
        return {opts[is_head], rank};
    }

    const __m128i mask_suff =  _mm_xor_si128(mask_pref, _mm_set1_epi32(-1));
    sym_mask = _mm_and_si128(sym_mask, mask_suff);
    bool has_sym = !_mm_testz_si128(sym_mask, sym_mask);
    run_id+=8;

    while(run_id<n_runs && !has_sym) {
        block = decode_block_sse42<vbyte_compressed,1,2>(stream);
        sym_mask = _mm_cmpeq_epi16(_mm_and_si128(block, alpha_mask), sym_vec);
        has_sym = !_mm_testz_si128(sym_mask, sym_mask);
        run_id+=8;
    }

    const auto matches = _mm_movemask_epi8(_mm_packs_epi16(sym_mask, _mm_setzero_si128()));
    const uint8_t first = __builtin_ctzl(matches);

    run_id = run_id - 8  + first;
    opts[1] = static_cast<int64_t>(run_id);

    return {opts[has_sym && run_id<n_runs], rank};
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint64_t> succ_sse42_32x4(const uint8_t ** stream, const size_t n_runs,
                                                            const uint8_t sigma, uint64_t idx, uint8_t symbol) {

    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const __m128i alpha_mask = _mm_set1_epi32(alpha_m);
    const __m128i sym_vec = _mm_set1_epi32(symbol);

    const uint8_t *prev_state = *stream;
    __m128i block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    //print32x4(block);
    __m128i bk_lengths =  shift_right_epi32(block, sigma_bits);
    //print32x4(bk_lengths);

    uint32_t prev_acc=0;
    uint64_t acc=hsum_epi32(bk_lengths);
    uint64_t rank=0;
    uint64_t run_id=0;

    while(acc<=idx){
        bk_lengths = _mm_and_si128(bk_lengths, _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec));
        //print32x4(bk_lengths);
        rank += hsum_epi32(bk_lengths);
        run_id+=4;

        prev_state = *stream;
        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        //print32x4(block);
        bk_lengths =  shift_right_epi32(block, sigma_bits);
        //print32x4(bk_lengths);

        prev_acc = acc;
        acc += hsum_epi32(bk_lengths);
    }

    idx-=prev_acc;

    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 4));
    bk_lengths = _mm_add_epi32(bk_lengths, _mm_slli_si128(bk_lengths, 8));

    //print32x4(bk_lengths);

    __m128i idx_mask = _mm_cmple_epu32(bk_lengths, _mm_set1_epi32(idx));//mask for >idx

    //print32x4(idx_mask);

    uint64_t less_than = _mm_cvtsi128_si64(_mm_packs_epi32(idx_mask, _mm_setzero_si128()));
    uint8_t idx_run = __builtin_ctzll(~less_than)>>4;

    const __m128i shuff_idxs = _mm_set_epi8(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
    const __m128i shuff = _mm_add_epi8(shuff_idxs, _mm_set1_epi8(idx_run<<2));

    //print8x16(shuff);

    const __m128i run_vec = _mm_shuffle_epi8(block, shuff);
    uint32_t run = _mm_cvtsi128_si32(run_vec);
    uint8_t last_symbol = run & alpha_m;

    const __m128i pf_sum_vec = _mm_shuffle_epi8(bk_lengths, shuff);
    uint32_t pf_sum = _mm_cvtsi128_si32(pf_sum_vec);

    *stream = prev_state;
    block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
    bk_lengths =  shift_right_epi32(block, sigma_bits);
    const __m128i mask_pref = _mm_loadu_si128((const __m128i*)mask32x4[idx_run+1]);
    bk_lengths = _mm_and_si128(bk_lengths, mask_pref);
    __m128i sym_mask = _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec);
    bk_lengths = _mm_and_si128(bk_lengths, sym_mask);

    rank += hsum_epi32(bk_lengths);
    const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
    rank -= (pf_sum-idx) * is_same_sym;

    int64_t opts[2] = {-1, (int64_t)run_id+idx_run};
    if(is_same_sym){
        const uint32_t len = run>>sigma_bits;
        return {opts[(pf_sum-len)==idx], rank};
    }

    const __m128i mask_suff =  _mm_xor_si128(mask_pref, _mm_set1_epi32(-1));
    sym_mask = _mm_and_si128(sym_mask, mask_suff);
    bool has_sym = !_mm_testz_si128(sym_mask, sym_mask);
    run_id+=4;

    while(run_id<n_runs && !has_sym) {
        block = decode_block_sse42<vbyte_compressed, 2, bytes_per_run>(stream);
        sym_mask = _mm_cmpeq_epi32(_mm_and_si128(block, alpha_mask), sym_vec);
        has_sym = !_mm_testz_si128(sym_mask, sym_mask);
        run_id+=4;
    }

    const uint64_t matches = _mm_cvtsi128_si64(_mm_packs_epi32(sym_mask, _mm_setzero_si128()));
    const uint8_t first = __builtin_ctzll(matches)>>4;

    run_id = run_id - 4  + first;
    opts[1] = static_cast<int64_t>(run_id);

    return {opts[has_sym && run_id<n_runs], rank};
}

template<uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint64_t> succ_sse42_64x2(const uint8_t ** stream, const size_t stream_bytes, const uint8_t sigma, uint64_t idx, uint8_t symbol){
    return {0, 0};
}

#endif //VLBT_SCAN_SSE42_H
