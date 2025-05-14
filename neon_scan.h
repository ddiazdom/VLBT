//
// Created by Diaz, Diego on 25.4.2025.
//

#ifndef BWT_DTS_BENCHMARKS_NEON_SCAN_H
#define BWT_DTS_BENCHMARKS_NEON_SCAN_H

#include <arm_neon.h>
#include "vbyte_simd_dec_tables.h"
#include "utils.h"

template<bool vbyte_compressed, uint8_t ctr_width, uint8_t bytes_per_run>
static inline uint8x16_t decode_block_neon(const uint8_t **stream){

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
        uint8x16_t compressed = vld1q_u8(*stream + 1);
        uint8x16_t dec_shuffle = vld1q_u8(pshuf);
        uint8x16_t data = vqtbl1q_u8(compressed, dec_shuffle);
        *stream += len + 1;

        return data;

    } else {
        uint8x16_t data = vld1q_u8(*stream);
        //align the bytes
        if constexpr (bytes_per_run==3){
            const int8x16_t dec_shuff = {0,1,2,-1, 3,4,5,-1, 6,7,8,-1, 9,10,11,-1};
            data = vqtbl1q_u8(data, dec_shuff);
            *stream+=12;
        } else if constexpr (bytes_per_run==5){
            const int8x16_t dec_shuff = {0,1,2,3,4,-1,-1,-1, 5,6,7,8,9,-1,-1,-1};
            data = vqtbl1q_u8(data, dec_shuff);
            *stream+=10;
        } else if constexpr (bytes_per_run==6){
            const int8x16_t dec_shuff = {0,1,2,3,4,5,-1,-1, 6,7,8,9,10,11,-1,-1};
            data = vqtbl1q_u8(data, dec_shuff);
            *stream+=12;
        } else if constexpr (bytes_per_run==7){
            const int8x16_t dec_shuff = {0,1,2,3,4,5,6,-1, 7,8,9,10,11,12,13,-1};
            data = vqtbl1q_u8(data, dec_shuff);
            *stream+=14;
        } else{
            *stream+=16;
        }
        return data;
    }
}

static inline void psum_epi8_ovf(uint8x16_t input, uint32_t idx, uint32_t& pf_sum, uint32_t& idx_run) {
    uint16x8_t halves[2];
    uint16x8_t idx_vec = vdupq_n_u16(idx);

    const int8x16_t l_shuff = {0,-1, 1,-1, 2,-1, 3,-1,
                               4,-1, 5,-1, 6,-1, 7,-1};
    const int8x16_t h_shuff ={8,-1, 9,-1, 10,-1, 11,-1,
                              12,-1, 13,-1, 14,-1, 15,-1};

    halves[0] = vqtbl1q_u8(input, l_shuff);
    halves[1] = vqtbl1q_u8(input, h_shuff);

    halves[0] = vaddq_u16(vextq_u16(vdupq_n_u16(0), halves[0], 7), halves[0]);
    halves[0] = vaddq_u16(vextq_u16(vdupq_n_u16(0), halves[0], 6), halves[0]);
    halves[0] = vaddq_u16(vextq_u16(vdupq_n_u16(0), halves[0], 4), halves[0]);

    uint16x8_t idx_mask = vcgtq_u16(halves[0], idx_vec);//mask for >idx
    uint8x8_t res = vshrn_n_u16(idx_mask, 4);
    uint64_t lt_low = vget_lane_u64(vreinterpret_u64_u8(res), 0);

    const int8x16_t pf_shuff = {14,15, -1,-1, -1,-1, -1,-1,
                                -1,-1, -1,-1, -1,-1, -1,-1};

    halves[1] = vaddq_u16(halves[1], vqtbl1q_u8(halves[0], pf_shuff));

    halves[1] = vaddq_u16(vextq_u16(vdupq_n_u16(0), halves[1], 7), halves[1]);
    halves[1] = vaddq_u16(vextq_u16(vdupq_n_u16(0), halves[1], 6), halves[1]);
    halves[1] = vaddq_u16(vextq_u16(vdupq_n_u16(0), halves[1], 4), halves[1]);

    idx_mask = vcgtq_u16(halves[1], idx_vec);//mask for >idx
    res = vshrn_n_u16(idx_mask, 4);
    uint64_t lt_high = vget_lane_u64(vreinterpret_u64_u8(res), 0);

    uint64_t lt = __builtin_ctzll(lt_low);
    lt += __builtin_ctzll(lt_high)*(lt==64);
    idx_run = lt>>3;

    const uint8x16_t shuff_idxs = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    const uint8x16_t ext_shuff = vaddq_u16(shuff_idxs, vdupq_n_u8((idx_run & 7)<<1));
    pf_sum = vgetq_lane_u16(vqtbl1q_u8(halves[idx_run>>3], ext_shuff), 0);
}

static inline void psum_epi16_ovf(uint16x8_t input, uint32_t idx, uint32_t& pf_sum, uint32_t& idx_run) {

    uint32x4_t halves[2];
    uint32x4_t idx_vec = vdupq_n_u32(idx);

    const int8x16_t l_shuff ={0,1,-1,-1, 2,3,-1,-1,
                              4,5,-1,-1, 6,7,-1,-1};

    const int8x16_t h_shuff ={8,9,-1,-1, 10,11,-1,-1,
                              12,13,-1,-1, 14,15,-1,-1};

    halves[0] = vqtbl1q_u8(input, l_shuff);
    halves[1] = vqtbl1q_u8(input, h_shuff);

    halves[0] = vaddq_u32(vextq_u32(vdupq_n_u32(0), halves[0], 3), halves[0]);
    halves[0] = vaddq_u32(vextq_u32(vdupq_n_u32(0), halves[0], 2), halves[0]);

    const int8x16_t idx_shuff = {0,4,8,12, -1,-1,-1,-1,
                                 -1,-1, -1,-1, -1,-1,-1,-1};
    uint32x4_t idx_mask = vcgtq_u32(halves[0], idx_vec);//mask for >idx
    auto lt_low = (uint64_t)vgetq_lane_u32(vqtbl1q_u8(idx_mask, idx_shuff), 0);

    const int8x16_t pf_shuff = {12,13,14,15, -1,-1,-1,-1,
                                -1,-1, -1,-1, -1,-1,-1,-1};
    halves[1] = vaddq_u32(halves[1], vqtbl1q_u8(halves[0], pf_shuff));

    halves[1] = vaddq_u32(vextq_u32(vdupq_n_u32(0), halves[1], 3), halves[1]);
    halves[1] = vaddq_u32(vextq_u32(vdupq_n_u32(0), halves[1], 2), halves[1]);

    idx_mask = vcgtq_u32(halves[1], idx_vec);//mask for >idx
    auto lt_high = (uint64_t)vgetq_lane_u32(vqtbl1q_u8(idx_mask, idx_shuff), 0);

    uint64_t lt = lt_low | lt_high<<32;
    idx_run = __builtin_ctzll(lt)>>3;

    const uint8x16_t shuff_idxs = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    const uint8x16_t ext_shuff = vaddq_u32(shuff_idxs, vdupq_n_u8((idx_run & 3)<<2));
    pf_sum = vgetq_lane_u32(vqtbl1q_u8(halves[idx_run>>2], ext_shuff), 0);
}

template<bool overflow16, bool overflow32=false>
static inline std::pair<uint64_t, uint8_t> inv_select_neon_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;
    const int8x16_t alpha_shift = vdupq_n_u8(-sigma_bits);

    size_t l=0;
    uint8x16_t block =  vld1q_u8(*stream);
    *stream+=16;
    uint8x16_t bk_lengths = vshlq_u8(block, alpha_shift);
    uint32_t prev_acc=0;

    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    uint32_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        block =  vld1q_u8(*stream);
        *stream+=16;
        bk_lengths = vshlq_u8(block, alpha_shift);

        prev_acc = acc;

        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        l++;
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    uint8_t alpha_m = (1UL << sigma_bits)-1, run, sym;
    const uint8x16_t alpha_mask = vdupq_n_u8(alpha_m);
    uint8x16_t sym_vec;

    if constexpr (overflow16){
        psum_epi8_ovf(bk_lengths, idx, pf_sum, idx_run);
        run = stream_start[(l*16)+idx_run];
        sym = run & alpha_m;
        sym_vec = vdupq_n_u8(sym);
    }else{
        //prefix sum
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 15), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 14), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 12), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 8), bk_lengths);
        //

        const uint8x16_t idx_mask = vcgtq_u8(bk_lengths, vdupq_n_u8(idx));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(vreinterpretq_u16_u8(idx_mask), 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>2;

        const uint8x16_t shuff = vdupq_n_u8(idx_run);
        const uint8x16_t run_vec = vqtbl1q_u8(block, shuff);
        run = vgetq_lane_u8(run_vec, 0);

        sym = run & alpha_m;
        sym_vec = vandq_u8(run_vec, alpha_mask);
        pf_sum = vgetq_lane_u8(vqtbl1q_u8(bk_lengths, shuff), 0);
    }

    pf_sum-= run>>sigma_bits;

    *stream = stream_start;
    //compute the rank answer
    uint64_t rank = 0;
    for(size_t p=0;p<l;p++){
        block = vld1q_u8(*stream);
        *stream+=16;
        bk_lengths = vshlq_u8(block, alpha_shift);

        bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));

        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    block = vld1q_u8(*stream);
    bk_lengths = vshlq_u8(block, alpha_shift);

    bk_lengths = vandq_u8(bk_lengths, vld1q_u8(mask8x16[idx_run]));
    bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    rank +=idx-pf_sum;
    return {rank, sym};
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline std::pair<uint64_t, uint8_t> inv_select_neon_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    const uint8_t sigma_bits = sym_width(sigma);
    const int16x8_t alpha_shift = vdupq_n_u16(-sigma_bits);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    uint16x8_t block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    uint16x8_t bk_lengths = vshlq_u16(block, alpha_shift);
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    uint64_t prev_acc=0, acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        bk_lengths = vshlq_u16(block, alpha_shift);

        prev_acc = acc;

        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        l++;
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    uint16_t alpha_m = (1UL << sigma_bits)-1;
    const uint16x8_t alpha_mask = vdupq_n_u16(alpha_m);

    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx, pf_sum, idx_run);
    } else {
        //prefix sum
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 7), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 6), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 4), bk_lengths);
        //

        const uint16x8_t idx_mask = vcgtq_u16(bk_lengths, vdupq_n_u16(idx));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(idx_mask, 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>3;
    }

    const uint8x16_t shuff_idxs = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    const uint8x16_t shuff = vaddq_u16(shuff_idxs, vdupq_n_u8(idx_run<<1));

    const uint16x8_t run_vec = vreinterpretq_u16_u8(vqtbl1q_u8(block, shuff));
    uint16_t run = vgetq_lane_u16(run_vec, 0);
    uint8_t sym = run & alpha_m;
    const uint16x8_t sym_vec = vandq_u16(run_vec, alpha_mask);

    if constexpr (!overflow8){
        pf_sum = vgetq_lane_u16(vreinterpretq_u16_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);
    }

    pf_sum-= run>>sigma_bits;

    *stream = stream_start;
    //compute the rank answer
    uint64_t rank = 0;
    for(size_t p=0;p<l;p++){

        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        bk_lengths = vshlq_u16(block, alpha_shift);

        bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));

        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    bk_lengths = vshlq_u16(block, alpha_shift);

    bk_lengths = vandq_u16(bk_lengths, vld1q_u16(mask16x8[idx_run]));
    bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    rank +=idx-pf_sum;
    return {rank, sym};
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint8_t> inv_select_neon_32x4(const uint8_t ** stream, uint8_t sigma, uint64_t idx){

    const uint8_t sigma_bits = sym_width(sigma);
    const int32x4_t alpha_shift = vdupq_n_u32(-sigma_bits);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    uint32x4_t block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    uint32x4_t bk_lengths = vshlq_u32(block, alpha_shift);
    uint64x2_t tmp = vpaddlq_u32(bk_lengths);
    uint64_t prev_acc=0, acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        bk_lengths = vshlq_u32(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(bk_lengths);
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
        l++;
    }

    idx-=prev_acc;

    //prefix sum
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 3), bk_lengths);
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 2), bk_lengths);
    //

    const uint32x4_t idx_mask = vcgtq_u32(bk_lengths, vdupq_n_u32(idx));// mask for >idx
    const uint16x4_t res = vshrn_n_u32(idx_mask, 16);
    const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u16(res), 0);
    uint8_t idx_run = __builtin_ctzll(less_than)>>4;
    uint32_t alpha_m = (1UL << sigma_bits)-1;
    const uint32x4_t alpha_mask = vdupq_n_u32(alpha_m);

    const uint8x16_t shuff_idxs = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    const uint8x16_t shuff = vaddq_u32(shuff_idxs, vdupq_n_u8(idx_run<<2));

    const uint32x4_t run_vec = vreinterpretq_u32_u8(vqtbl1q_u8(block, shuff));
    const uint32x4_t sym_vec = vandq_u32(run_vec, alpha_mask);

    uint32_t run = vgetq_lane_u32(run_vec, 0);
    uint32_t pf_sum = vgetq_lane_u32(vreinterpretq_u32_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);
    pf_sum-= run>>sigma_bits;
    uint8_t sym = run & alpha_m;

    *stream = stream_start;

    //compute the rank answer
    uint64_t rank = 0;
    for(size_t p=0;p<l;p++){
        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        bk_lengths = vshlq_u32(block, alpha_shift);

        bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));

        tmp = vpaddlq_u32(bk_lengths);
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    bk_lengths = vshlq_u32(block, alpha_shift);

    bk_lengths = vandq_u32(bk_lengths, vld1q_u32(mask32x4[idx_run]));
    bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(bk_lengths);
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    rank +=idx-pf_sum;

    return {rank, sym};
}

template<uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint8_t> inv_select_neon_64x2(const uint8_t ** stream, uint8_t sigma, uint64_t idx){
    return {0,0};
}


template<bool overflow16, bool overflow32=false>
static inline uint8_t access_neon_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t *stream_start = *stream;
    const int8x16_t alpha_shift = vdupq_n_u8(-sigma_bits);

    size_t l=0;
    uint8x16_t block =  vld1q_u8(*stream);
    *stream+=16;
    uint8x16_t bk_lengths = vshlq_u8(block, alpha_shift);
    uint64_t prev_acc=0;

    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        block =  vld1q_u8(*stream);
        *stream+=16;
        bk_lengths = vshlq_u8(block, alpha_shift);

        prev_acc = acc;

        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        l++;
    }

    idx-=prev_acc;
    uint32_t idx_run;
    uint8_t alpha_m = (1UL << sigma_bits)-1, run;

    if constexpr (overflow16){
        uint32_t pf_sum;
        psum_epi8_ovf(bk_lengths, idx, pf_sum, idx_run);
    }else {
        //prefix sum
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 15), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 14), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 12), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 8), bk_lengths);
        //
        const uint8x16_t idx_mask = vcgtq_u8(bk_lengths, vdupq_n_u8(idx));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(vreinterpretq_u16_u8(idx_mask), 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>2;
    }
    run = stream_start[(l*16)+idx_run];
    return run & alpha_m;
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline uint8_t access_neon_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    const uint8_t sigma_bits = sym_width(sigma);
    const int16x8_t alpha_shift = vdupq_n_u16(-sigma_bits);

    size_t l=0;
    uint16x8_t block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    uint16x8_t bk_lengths = vshlq_u16(block, alpha_shift);
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    uint64_t prev_acc=0, acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        bk_lengths = vshlq_u16(block, alpha_shift);

        prev_acc = acc;

        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        l++;
    }

    idx-=prev_acc;
    uint32_t idx_run;

    uint16_t alpha_m = (1UL << sigma_bits)-1;

    if constexpr (overflow8){
        uint32_t pf_sum;
        psum_epi16_ovf(bk_lengths, idx, pf_sum, idx_run);
    }else {
        //prefix sum
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 7), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 6), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 4), bk_lengths);
        //

        const uint16x8_t idx_mask = vcgtq_u16(bk_lengths, vdupq_n_u16(idx));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(idx_mask, 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than) >> 3;
    }

    const uint8x16_t shuff_idxs = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    const uint8x16_t shuff = vaddq_u16(shuff_idxs, vdupq_n_u8(idx_run<<1));

    const uint16x8_t run_vec = vreinterpretq_u16_u8(vqtbl1q_u8(block, shuff));
    return vgetq_lane_u16(run_vec, 0) & alpha_m;
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline uint8_t access_neon_32x4(const uint8_t ** stream, uint8_t sigma, uint64_t idx){

    //TODO: assert idx fits 4 bytes
    const uint8_t sigma_bits = sym_width(sigma);
    const int32x4_t alpha_shift = vdupq_n_u32(-sigma_bits);

    size_t l=0;
    uint32x4_t block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    uint32x4_t bk_lengths = vshlq_u32(block, alpha_shift);
    uint64x2_t tmp = vpaddlq_u32(bk_lengths);
    uint64_t prev_acc=0, acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        bk_lengths = vshlq_u32(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(bk_lengths);
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
        l++;
    }

    idx-=prev_acc;

    //prefix sum
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 3), bk_lengths);
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 2), bk_lengths);
    //

    const uint32x4_t idx_mask = vcgtq_u32(bk_lengths, vdupq_n_u32(idx));// mask for >idx
    const uint16x4_t res = vshrn_n_u32(idx_mask, 16);
    const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u16(res), 0);
    uint8_t idx_run = __builtin_ctzll(less_than)>>4;

    uint32_t alpha_m = (1UL << sigma_bits)-1;
    const uint8x16_t shuff_idxs = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    const uint8x16_t shuff = vaddq_u32(shuff_idxs, vdupq_n_u8(idx_run<<2));
    const uint32x4_t run_vec = vreinterpretq_u32_u8(vqtbl1q_u8(block, shuff));

    return vgetq_lane_u32(run_vec, 0) & alpha_m;
}

template<uint8_t bytes_per_run>
static inline uint8_t access_neon_64x2(const uint8_t ** stream, uint8_t sigma, uint64_t idx){
    return 0;
}

template<bool overflow16, bool overflow32=false>
static inline int64_t rank_neon_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx, uint8_t symbol){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    uint8_t alpha_m = (1UL << sigma_bits)-1;
    const int8x16_t alpha_shift = vdupq_n_u8(-sigma_bits);
    const uint8x16_t alpha_mask = vdupq_n_u8(alpha_m);
    uint8x16_t sym_vec = vdupq_n_u8(symbol);

    uint8x16_t block = vld1q_u8(*stream);
    *stream+=16;
    uint8x16_t bk_lengths = vshlq_u8(block, alpha_shift);

    uint32_t prev_acc=0;
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    uint32_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    uint32_t prev_rank=0;
    uint8x16_t sym_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));
    tmp  = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(sym_lengths)));
    uint32_t rank = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    size_t l=0;
    while(acc<=idx){
        block =  vld1q_u8(*stream);
        *stream+=16;
        bk_lengths = vshlq_u8(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        prev_rank = rank;
        sym_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(sym_lengths)));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
        l++;
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr (overflow16){
        psum_epi8_ovf(bk_lengths, idx, pf_sum, idx_run);
    }else{
        //prefix sum
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 15), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 14), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 12), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 8), bk_lengths);
        //

        const uint8x16_t idx_mask = vcgtq_u8(bk_lengths, vdupq_n_u8(idx));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(vreinterpretq_u16_u8(idx_mask), 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>2;

        const uint8x16_t shuff = vdupq_n_u8(idx_run);
        pf_sum = vgetq_lane_u8(vqtbl1q_u8(bk_lengths, shuff), 0);
    }

    *stream -=16;
    uint8_t run = (*stream)[idx_run];
    uint8_t last_symbol = run & alpha_m;

    block = vld1q_u8(*stream);
    sym_lengths = vshlq_u8(block, alpha_shift);
    sym_lengths = vandq_u8(sym_lengths, vld1q_u8(mask8x16[idx_run+1]));
    sym_lengths = vandq_u8(sym_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(sym_lengths)));
    rank = prev_rank + vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    rank -=(pf_sum-idx) * (last_symbol==symbol);
    return rank;
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline int64_t rank_neon_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx, uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const int16x8_t alpha_shift = vdupq_n_u16(-sigma_bits);
    uint16_t alpha_m = (1UL << sigma_bits)-1;
    const uint16x8_t alpha_mask = vdupq_n_u16(alpha_m);
    const uint16x8_t sym_vec = vdupq_n_u16(symbol);

    const uint8_t *prev_state = *stream;

    uint16x8_t block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    uint16x8_t bk_lengths = vshlq_u16(block, alpha_shift);
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(bk_lengths));

    uint64_t prev_acc=0;
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    uint64_t prev_rank=0;
    uint16x8_t sym_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));
    tmp = vpaddlq_u32(vpaddlq_u16(sym_lengths));
    uint64_t rank = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        prev_state = *stream;
        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        bk_lengths = vshlq_u16(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        sym_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(sym_lengths));
        prev_rank = rank;
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx-=prev_acc;
    uint32_t idx_run, pf_sum;

    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx, pf_sum, idx_run);
    } else {
        //prefix sum
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 7), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 6), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 4), bk_lengths);
        //

        const uint16x8_t idx_mask = vcgtq_u16(bk_lengths, vdupq_n_u16(idx));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(idx_mask, 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>3;
    }

    const uint8x16_t shuff_idxs = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    const uint8x16_t shuff = vaddq_u16(shuff_idxs, vdupq_n_u8(idx_run<<1));

    const uint16x8_t run_vec = vreinterpretq_u16_u8(vqtbl1q_u8(block, shuff));
    uint16_t run = vgetq_lane_u16(run_vec, 0);
    uint8_t last_symbol = run & alpha_m;

    if constexpr (!overflow8){
        pf_sum = vgetq_lane_u16(vreinterpretq_u16_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);
    }

    *stream = prev_state;
    block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    sym_lengths = vshlq_u16(block, alpha_shift);
    sym_lengths = vandq_u16(sym_lengths, vld1q_u16(mask16x8[idx_run+1]));
    sym_lengths = vandq_u16(sym_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(sym_lengths));
    rank = prev_rank + vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    rank -=(pf_sum-idx)*(last_symbol==symbol);

    return (int64_t)rank;
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline int64_t rank_neon_32x4(const uint8_t ** stream, uint8_t sigma, uint64_t idx, uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const int32x4_t alpha_shift = vdupq_n_u32(-sigma_bits);
    uint32_t alpha_m = (1UL << sigma_bits)-1;
    const uint32x4_t alpha_mask = vdupq_n_u32(alpha_m);
    const uint32x4_t sym_vec = vdupq_n_u32(symbol);

    const uint8_t *prev_state = *stream;

    uint32x4_t block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    uint32x4_t bk_lengths = vshlq_u32(block, alpha_shift);
    uint64x2_t tmp = vpaddlq_u32(bk_lengths);

    uint64_t prev_acc=0;
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    uint64_t prev_rank=0;
    uint32x4_t sym_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));
    tmp = vpaddlq_u32(sym_lengths);
    uint64_t rank = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        prev_state = *stream;
        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        bk_lengths = vshlq_u32(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(bk_lengths);
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        prev_rank = rank;
        sym_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(sym_lengths);
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx-=prev_acc;

    //prefix sum
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 3), bk_lengths);
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 2), bk_lengths);
    //

    const uint32x4_t idx_mask = vcgtq_u32(bk_lengths, vdupq_n_u32(idx));// mask for >idx
    const uint16x4_t res = vshrn_n_u32(idx_mask, 16);
    const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u16(res), 0);
    uint8_t idx_run = __builtin_ctzll(less_than)>>4;

    const uint8x16_t shuff_idxs = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    const uint8x16_t shuff = vaddq_u32(shuff_idxs, vdupq_n_u8(idx_run<<2));

    const uint32x4_t run_vec = vreinterpretq_u32_u8(vqtbl1q_u8(block, shuff));
    uint32_t run = vgetq_lane_u32(run_vec, 0);
    uint8_t last_symbol = run & alpha_m;

    uint32_t pf_sum = vgetq_lane_u32(vreinterpretq_u32_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);

    *stream = prev_state;
    block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    sym_lengths = vshlq_u32(block, alpha_shift);
    sym_lengths = vandq_u32(sym_lengths, vld1q_u32(mask32x4[idx_run+1]));
    sym_lengths = vandq_u32(sym_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(sym_lengths);
    rank = prev_rank + vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    rank -=(pf_sum-idx)*(last_symbol==symbol);

    return (int64_t)rank;
}

template<uint8_t bytes_per_run>
static inline int64_t rank_neon_64x2(const uint8_t ** stream, uint8_t sigma, uint64_t idx, uint8_t symbol){
    return 0;
}

#endif //BWT_DTS_BENCHMARKS_NEON_SCAN_H
