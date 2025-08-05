//
// Created by Diaz, Diego on 25.4.2025.
//

#ifndef VLBT_SCAN_NEON_H
#define VLBT_SCAN_NEON_H

#include <arm_neon.h>
#include "simd_tables.h"
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

static inline void psum_epi8_ovf(uint8x16_t input, const uint32_t idx, uint32_t& pf_sum, uint32_t& idx_run) {
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

template<bool overflow16, bool overflow32, bool get_run_id=false>
static inline auto inv_select_neon_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    //std::pair<uint64_t, uint8_t>//without run_id
    //std::tuple<uint64_t, uint8_t, int64_t>//with run_id

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

    if constexpr (get_run_id){
        const auto run_id= l*16 + idx_run;
        int64_t options[2] = {-1, static_cast<int64_t>(run_id)};
        return std::make_tuple(rank, sym, options[idx==pf_sum]);
    }else{
        return std::make_pair(rank, sym);
    }
}

template<bool vbyte_compressed, bool overflow8, bool overflow16, bool get_run_id=false>
static inline auto inv_select_neon_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx){

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
    const uint16_t run = vgetq_lane_u16(run_vec, 0);
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

    if constexpr (get_run_id){
        const auto run_id= l*8 + idx_run;
        int64_t options[2] = {-1, static_cast<int64_t>(run_id)};
        return std::make_tuple(rank, sym, options[idx==pf_sum]);
    }else{
        return std::make_pair(rank, sym);
    }
}

template<bool vbyte_compressed, uint8_t bytes_per_run, bool get_run_id=false>
static inline auto inv_select_neon_32x4(const uint8_t ** stream, const uint8_t sigma, uint64_t idx){

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
    const uint8_t idx_run = __builtin_ctzll(less_than)>>4;
    const uint32_t alpha_m = (1UL << sigma_bits)-1;
    const uint32x4_t alpha_mask = vdupq_n_u32(alpha_m);

    const uint8x16_t shuff_idxs = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    const uint8x16_t shuff = vaddq_u32(shuff_idxs, vdupq_n_u8(idx_run<<2));

    const uint32x4_t run_vec = vreinterpretq_u32_u8(vqtbl1q_u8(block, shuff));
    const uint32x4_t sym_vec = vandq_u32(run_vec, alpha_mask);

    const uint32_t run = vgetq_lane_u32(run_vec, 0);
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

    if constexpr (get_run_id){
        const auto run_id= l*4 + idx_run;
        int64_t options[2] = {-1, static_cast<int64_t>(run_id)};
        return std::make_tuple(rank, sym, options[idx==pf_sum]);
    }else{
        return std::make_pair(rank, sym);
    }
}

template<uint8_t bytes_per_run, bool get_run_id=false>
static inline auto inv_select_neon_64x2(const uint8_t ** stream, uint8_t sigma, uint64_t idx){
    if constexpr (get_run_id){
        return std::make_tuple<int64_t, uint8_t, uint64_t>(0,0, 0);
    }else{
        return std::make_pair<int64_t, uint8_t>(0,0);
    }
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

    const uint8_t sigma_bits = sym_width(sigma);
    const int32x4_t alpha_shift = vdupq_n_u32(-sigma_bits);

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
static inline std::pair<uint64_t, uint64_t> get_phi_run_neon_8x16(const uint8_t **stream, uint64_t idx){

    //NOTE here I do not need to vbyte compress the block
    uint8x16_t block =  vld1q_u8(*stream);
    *stream+=16;
    uint64_t prev_acc=0;
    uint32_t idx_run=0;

    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(block)));
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    while(acc<=idx){
        block =  vld1q_u8(*stream);
        *stream+=16;
        idx_run+=16;
        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(block)));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx-=prev_acc;

    uint32_t pf_sum, tmp_idx_run;
    if constexpr (overflow16){
        psum_epi8_ovf(block, idx, pf_sum, tmp_idx_run);
    }else {
        //prefix sum
        block = vaddq_u8(vextq_u8(vdupq_n_u8(0), block, 15), block);
        block = vaddq_u8(vextq_u8(vdupq_n_u8(0), block, 14), block);
        block = vaddq_u8(vextq_u8(vdupq_n_u8(0), block, 12), block);
        block = vaddq_u8(vextq_u8(vdupq_n_u8(0), block, 8), block);
        //
        const uint8x16_t idx_mask = vcgtq_u8(block, vdupq_n_u8(idx));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(vreinterpretq_u16_u8(idx_mask), 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        tmp_idx_run = __builtin_ctzll(less_than)>>2;

        const uint8x16_t shuff = vdupq_n_u8(tmp_idx_run);
        pf_sum = vgetq_lane_u8(vqtbl1q_u8(block, shuff), 0);
    }
    const uint32_t offset = static_cast<uint32_t>((*stream - 16)[tmp_idx_run]) - (pf_sum-idx);
    return std::make_pair(idx_run+tmp_idx_run, offset);
}

template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline std::pair<uint64_t, uint64_t> get_phi_run_neon_16x8(const uint8_t **stream, uint64_t idx){

    uint16x8_t block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(block));
    uint64_t prev_acc=0, acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    uint32_t idx_run=0;

    while(acc<=idx){
        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        idx_run+=8;
        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(block));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx-=prev_acc;
    uint32_t pf_sum, tmp_idx_run, len;

    static const uint8x16_t shuff_idxs = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    if constexpr (overflow8){
        psum_epi16_ovf(block, idx, pf_sum, tmp_idx_run);
        const uint8x16_t shuff = vaddq_u16(shuff_idxs, vdupq_n_u8(tmp_idx_run<<1));
        const uint16x8_t run_vec = vreinterpretq_u16_u8(vqtbl1q_u8(block, shuff));
        len = vgetq_lane_u16(run_vec, 0);
    } else {
        //prefix sum
        uint16x8_t pf_sum_vec = vaddq_u16(vextq_u16(vdupq_n_u16(0), block, 7), block);
        pf_sum_vec = vaddq_u16(vextq_u16(vdupq_n_u16(0), pf_sum_vec, 6), pf_sum_vec);
        pf_sum_vec = vaddq_u16(vextq_u16(vdupq_n_u16(0), pf_sum_vec, 4), pf_sum_vec);
        //

        const uint16x8_t idx_mask = vcgtq_u16(pf_sum_vec, vdupq_n_u16(idx));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(idx_mask, 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        tmp_idx_run = __builtin_ctzll(less_than) >> 3;

        const uint8x16_t shuff = vaddq_u16(shuff_idxs, vdupq_n_u8(tmp_idx_run<<1));

        pf_sum = vgetq_lane_u16(vreinterpretq_u16_u8(vqtbl1q_u8(pf_sum_vec, shuff)), 0);

        const uint16x8_t run_vec = vreinterpretq_u16_u8(vqtbl1q_u8(block, shuff));
        len = vgetq_lane_u16(run_vec, 0);
    }

    uint64_t offset = len - (pf_sum-idx);
    return std::make_pair(idx_run+tmp_idx_run, offset);
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint64_t> get_phi_run_neon_32x4(const uint8_t ** stream, uint64_t idx){

    uint32x4_t block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    uint64x2_t tmp = vpaddlq_u32(block);
    uint64_t prev_acc=0, acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    uint32_t idx_run=0;

    while(acc<=idx){
        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        prev_acc = acc;
        idx_run+=4;
        tmp = vpaddlq_u32(block);
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx-=prev_acc;

    //prefix sum
    uint32x4_t pf_sum_vec = vaddq_u32(vextq_u32(vdupq_n_u32(0), block, 3), block);
    pf_sum_vec = vaddq_u32(vextq_u32(vdupq_n_u32(0), pf_sum_vec, 2), pf_sum_vec);
    //

    const uint32x4_t idx_mask = vcgtq_u32(pf_sum_vec, vdupq_n_u32(idx));// mask for >idx
    const uint16x4_t res = vshrn_n_u32(idx_mask, 16);
    const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u16(res), 0);
    const uint8_t tmp_idx_run = __builtin_ctzll(less_than)>>4;

    const uint8x16_t shuff_idxs = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    const uint8x16_t shuff = vaddq_u32(shuff_idxs, vdupq_n_u8(tmp_idx_run<<2));

    const uint32_t pf_sum = vgetq_lane_u32(vreinterpretq_u32_u8(vqtbl1q_u8(pf_sum_vec, shuff)), 0);

    const uint32x4_t run_vec = vreinterpretq_u32_u8(vqtbl1q_u8(block, shuff));

    const uint32_t offset = vgetq_lane_u32(run_vec, 0) - static_cast<uint32_t>(pf_sum - idx);
    return std::make_pair(idx_run+tmp_idx_run, offset);
}

template<uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint64_t> get_phi_run_neon_64x2(const uint8_t ** stream, uint64_t idx){
    return std::make_pair(0, 0);
}

template<bool overflow16, bool overflow32, bool check_head>
static inline int64_t rank_neon_8x16(const uint8_t **stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const int8x16_t alpha_shift = vdupq_n_u8(-sigma_bits);
    const uint8x16_t alpha_mask = vdupq_n_u8(alpha_m);
    const uint8x16_t sym_vec = vdupq_n_u8(symbol);

    uint8x16_t block = vld1q_u8(*stream);
    *stream+=16;
    uint8x16_t bk_lengths = vshlq_u8(block, alpha_shift);

    uint32_t prev_acc=0;
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    uint32_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    uint32_t rank=0;
    while(acc<=idx){
        bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block =  vld1q_u8(*stream);
        *stream+=16;
        bk_lengths = vshlq_u8(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
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

    const uint8_t run = (*stream-16)[idx_run];
    const uint8_t last_symbol = run & alpha_m;

    //block = vld1q_u8(*stream);
    bk_lengths = vshlq_u8(block, alpha_shift);
    bk_lengths = vandq_u8(bk_lengths, vld1q_u8(mask8x16[idx_run+1]));
    bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    if constexpr (check_head) {
        const uint8_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx;//check if idx is the head of the run
        rank -= (pf_sum-idx) * is_same_sym;
        rank = (rank<<1) | is_head;
        rank = (rank<<1) | is_same_sym;
    } else {
        rank -=(pf_sum-idx) * (last_symbol==symbol);
    }
    return rank;
}

template<bool vbyte_compressed, bool overflow8, bool overflow16, bool check_head>
static inline int64_t rank_neon_16x8(const uint8_t **stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const int16x8_t alpha_shift = vdupq_n_u16(-sigma_bits);
    const uint16_t alpha_m = (1UL << sigma_bits)-1;
    const uint16x8_t alpha_mask = vdupq_n_u16(alpha_m);
    const uint16x8_t sym_vec = vdupq_n_u16(symbol);

    uint16x8_t block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    uint16x8_t bk_lengths = vshlq_u16(block, alpha_shift);

    uint64_t prev_acc=0;
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    uint64_t rank = 0;

    while(acc<=idx){
        bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        bk_lengths = vshlq_u16(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
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
    const uint16_t run = vgetq_lane_u16(run_vec, 0);
    const uint8_t last_symbol = run & alpha_m;

    if constexpr (!overflow8){
        pf_sum = vgetq_lane_u16(vreinterpretq_u16_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);
    }

    bk_lengths = vshlq_u16(block, alpha_shift);
    bk_lengths = vandq_u16(bk_lengths, vld1q_u16(mask16x8[idx_run+1]));
    bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    if constexpr (check_head) {
        const uint16_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx;//check if idx is the head of the run
        rank -=(pf_sum-idx) * is_same_sym;
        rank = (rank<<1) | is_head;
        rank = (rank<<1) | is_same_sym;
    } else {
        rank -=(pf_sum-idx) * (last_symbol==symbol);
    }
    return (int64_t)rank;
}

template<bool vbyte_compressed, uint8_t bytes_per_run, bool check_head>
static inline int64_t rank_neon_32x4(const uint8_t ** stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const int32x4_t alpha_shift = vdupq_n_u32(-sigma_bits);
    const uint32_t alpha_m = (1UL << sigma_bits)-1;
    const uint32x4_t alpha_mask = vdupq_n_u32(alpha_m);
    const uint32x4_t sym_vec = vdupq_n_u32(symbol);

    uint32x4_t block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    uint32x4_t bk_lengths = vshlq_u32(block, alpha_shift);
    uint64x2_t tmp = vpaddlq_u32(bk_lengths);

    uint64_t prev_acc=0;
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    uint64_t rank = 0;

    while(acc<=idx){
        bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(bk_lengths);
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        bk_lengths = vshlq_u32(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(bk_lengths);
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx-=prev_acc;

    //prefix sum
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 3), bk_lengths);
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 2), bk_lengths);
    //

    const uint32x4_t idx_mask = vcgtq_u32(bk_lengths, vdupq_n_u32(idx));// mask for >idx
    const uint16x4_t res = vshrn_n_u32(idx_mask, 16);
    const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u16(res), 0);
    const uint8_t idx_run = __builtin_ctzll(less_than)>>4;

    const uint8x16_t shuff_idxs = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
    const uint8x16_t shuff = vaddq_u32(shuff_idxs, vdupq_n_u8(idx_run<<2));

    const uint32x4_t run_vec = vreinterpretq_u32_u8(vqtbl1q_u8(block, shuff));
    const uint32_t run = vgetq_lane_u32(run_vec, 0);
    const uint8_t last_symbol = run & alpha_m;

    const uint32_t pf_sum = vgetq_lane_u32(vreinterpretq_u32_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);

    bk_lengths = vshlq_u32(block, alpha_shift);
    bk_lengths = vandq_u32(bk_lengths, vld1q_u32(mask32x4[idx_run+1]));
    bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(bk_lengths);
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    if constexpr (check_head) {
        const uint32_t len = run >> sigma_bits;//get the length of the run where idx falls
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

template<uint8_t bytes_per_run, bool check_head>
static inline int64_t rank_neon_64x2(const uint8_t ** stream, const uint8_t sigma, uint64_t idx, const uint8_t symbol){
    return 0;
}

template<bool overflow16, bool overflow32, bool check_head>
static inline std::pair<uint64_t, uint64_t> range_rank_neon_8x16(const uint8_t **stream, const uint8_t sigma,
                                                                 uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const uint8_t alpha_m = (1UL << sigma_bits)-1;
    const int8x16_t alpha_shift = vdupq_n_u8(-sigma_bits);
    const uint8x16_t alpha_mask = vdupq_n_u8(alpha_m);
    const uint8x16_t sym_vec = vdupq_n_u8(symbol);

    uint8x16_t block = vld1q_u8(*stream);
    *stream+=16;
    uint8x16_t bk_lengths = vshlq_u8(block, alpha_shift);

    uint32_t prev_acc = 0;
    uint64x2_t tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    uint32_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    uint32_t rank_i=0;
    while(acc<=idx_i){
        bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        rank_i += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block =  vld1q_u8(*stream);
        *stream+=16;
        bk_lengths = vshlq_u8(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    uint32_t rank_j = rank_i;

    idx_i-=prev_acc;

    uint32_t idx_run, pf_sum;
    if constexpr (overflow16){
        psum_epi8_ovf(bk_lengths, idx_i, pf_sum, idx_run);
    }else{
        //prefix sum
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 15), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 14), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 12), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 8), bk_lengths);
        //

        const uint8x16_t idx_mask = vcgtq_u8(bk_lengths, vdupq_n_u8(idx_i));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(vreinterpretq_u16_u8(idx_mask), 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>2;

        const uint8x16_t shuff = vdupq_n_u8(idx_run);
        pf_sum = vgetq_lane_u8(vqtbl1q_u8(bk_lengths, shuff), 0);
    }

    uint8_t run = (*stream-16)[idx_run];
    uint8_t last_symbol = run & alpha_m;

    bk_lengths = vshlq_u8(block, alpha_shift);
    bk_lengths = vandq_u8(bk_lengths, vld1q_u8(mask8x16[idx_run+1]));
    bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    rank_i += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    if constexpr (check_head) {
        const uint8_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank_i -= (pf_sum-idx_i) * is_same_sym;
        rank_i = (rank_i<<1) | is_head;
        rank_i = (rank_i<<1) | is_same_sym;
    } else {
        rank_i -=(pf_sum-idx_i) * (last_symbol==symbol);
    }

    //now we deal with j
    //==================
    bk_lengths = vshlq_u8(block, alpha_shift);
    //acc and block are still valid, no need to recompute them
    while(acc<=idx_j){
        bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        rank_j += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block =  vld1q_u8(*stream);
        *stream+=16;
        bk_lengths = vshlq_u8(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx_j-=prev_acc;

    if constexpr (overflow16){
        psum_epi8_ovf(bk_lengths, idx_j, pf_sum, idx_run);
    }else{
        //prefix sum
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 15), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 14), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 12), bk_lengths);
        bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 8), bk_lengths);
        //

        const uint8x16_t idx_mask = vcgtq_u8(bk_lengths, vdupq_n_u8(idx_j));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(vreinterpretq_u16_u8(idx_mask), 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>2;

        const uint8x16_t shuff = vdupq_n_u8(idx_run);
        pf_sum = vgetq_lane_u8(vqtbl1q_u8(bk_lengths, shuff), 0);
    }

    run = (*stream-16)[idx_run];
    last_symbol = run & alpha_m;

    bk_lengths = vshlq_u8(block, alpha_shift);
    bk_lengths = vandq_u8(bk_lengths, vld1q_u8(mask8x16[idx_run+1]));
    bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    rank_j += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    rank_j -=(pf_sum-idx_j) * (last_symbol==symbol);

    return std::make_pair(rank_i, rank_j);
}

template<bool vbyte_compressed, bool overflow8, bool overflow16, bool check_head>
static inline std::pair<uint64_t, uint64_t> range_rank_neon_16x8(const uint8_t **stream, const uint8_t sigma,
                                                                 uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const int16x8_t alpha_shift = vdupq_n_u16(-sigma_bits);
    const uint16_t alpha_m = (1UL << sigma_bits)-1;
    const uint16x8_t alpha_mask = vdupq_n_u16(alpha_m);
    const uint16x8_t sym_vec = vdupq_n_u16(symbol);

    uint16x8_t block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    uint16x8_t bk_lengths = vshlq_u16(block, alpha_shift);

    uint64_t prev_acc=0;
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    uint64_t rank_i = 0;

    while(acc<=idx_i){
        bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        rank_i += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        bk_lengths = vshlq_u16(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    uint64_t rank_j = rank_i;

    idx_i-=prev_acc;
    uint32_t idx_run, pf_sum;
    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx_i, pf_sum, idx_run);
    } else {
        //prefix sum
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 7), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 6), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 4), bk_lengths);
        //

        const uint16x8_t idx_mask = vcgtq_u16(bk_lengths, vdupq_n_u16(idx_i));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(idx_mask, 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>3;
    }

    const uint8x16_t shuff_idxs = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};

    uint8x16_t shuff = vaddq_u16(shuff_idxs, vdupq_n_u8(idx_run<<1));
    uint16x8_t run_vec = vreinterpretq_u16_u8(vqtbl1q_u8(block, shuff));
    uint16_t run = vgetq_lane_u16(run_vec, 0);
    uint8_t last_symbol = run & alpha_m;

    if constexpr (!overflow8) {
        pf_sum = vgetq_lane_u16(vreinterpretq_u16_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);
    }

    bk_lengths = vshlq_u16(block, alpha_shift);
    bk_lengths = vandq_u16(bk_lengths, vld1q_u16(mask16x8[idx_run+1]));
    bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    rank_i += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    if constexpr (check_head) {
        uint16_t len = run >> sigma_bits;//get the length of the run where idx falls
        bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank_i -=(pf_sum-idx_i) * is_same_sym;
        rank_i = (rank_i<<1) | is_head;
        rank_i = (rank_i<<1) | is_same_sym;
    } else {
        rank_i -=(pf_sum-idx_i) * (last_symbol==symbol);
    }

    //=== deal with j
    bk_lengths = vshlq_u16(block, alpha_shift);
    //acc and block are still valid, no need to recompute them
    while(acc<=idx_j){
        bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        rank_j += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        bk_lengths = vshlq_u16(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx_j-=prev_acc;
    if constexpr (overflow8){
        psum_epi16_ovf(bk_lengths, idx_j, pf_sum, idx_run);
    } else {
        //prefix sum
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 7), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 6), bk_lengths);
        bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 4), bk_lengths);
        //

        const uint16x8_t idx_mask = vcgtq_u16(bk_lengths, vdupq_n_u16(idx_j));//mask for >idx
        const uint8x8_t res = vshrn_n_u16(idx_mask, 4);
        const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        idx_run = __builtin_ctzll(less_than)>>3;
    }

    shuff = vaddq_u16(shuff_idxs, vdupq_n_u8(idx_run<<1));
    run_vec = vreinterpretq_u16_u8(vqtbl1q_u8(block, shuff));
    run = vgetq_lane_u16(run_vec, 0);
    last_symbol = run & alpha_m;

    if constexpr (!overflow8) {
        pf_sum = vgetq_lane_u16(vreinterpretq_u16_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);
    }

    bk_lengths = vshlq_u16(block, alpha_shift);
    bk_lengths = vandq_u16(bk_lengths, vld1q_u16(mask16x8[idx_run+1]));
    bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    rank_j += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    rank_j -= (pf_sum-idx_j) * (last_symbol==symbol);

    return std::make_pair(rank_i, rank_j);
}

template<bool vbyte_compressed, uint8_t bytes_per_run, bool check_head>
static inline std::pair<uint64_t, uint64_t> range_rank_neon_32x4(const uint8_t ** stream, const uint8_t sigma,
                                                                uint64_t idx_i, uint64_t idx_j, const uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const int32x4_t alpha_shift = vdupq_n_u32(-sigma_bits);
    const uint32_t alpha_m = (1UL << sigma_bits)-1;
    const uint32x4_t alpha_mask = vdupq_n_u32(alpha_m);
    const uint32x4_t sym_vec = vdupq_n_u32(symbol);

    uint32x4_t block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    uint32x4_t bk_lengths = vshlq_u32(block, alpha_shift);
    uint64x2_t tmp = vpaddlq_u32(bk_lengths);

    uint64_t prev_acc=0;
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    uint64_t rank_i = 0;

    while(acc<=idx_i){
        bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(bk_lengths);
        rank_i += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        bk_lengths = vshlq_u32(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(bk_lengths);
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    uint64_t rank_j = rank_i;
    idx_i-=prev_acc;

    //prefix sum
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 3), bk_lengths);
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 2), bk_lengths);
    //

    uint32x4_t idx_mask = vcgtq_u32(bk_lengths, vdupq_n_u32(idx_i));// mask for >idx
    uint16x4_t res = vshrn_n_u32(idx_mask, 16);
    uint64_t less_than = vget_lane_u64(vreinterpret_u64_u16(res), 0);
    uint8_t idx_run = __builtin_ctzll(less_than)>>4;

    const uint8x16_t shuff_idxs = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};

    uint8x16_t shuff = vaddq_u32(shuff_idxs, vdupq_n_u8(idx_run<<2));
    uint32x4_t run_vec = vreinterpretq_u32_u8(vqtbl1q_u8(block, shuff));
    uint32_t run = vgetq_lane_u32(run_vec, 0);
    uint8_t last_symbol = run & alpha_m;

    uint32_t pf_sum = vgetq_lane_u32(vreinterpretq_u32_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);

    bk_lengths = vshlq_u32(block, alpha_shift);
    bk_lengths = vandq_u32(bk_lengths, vld1q_u32(mask32x4[idx_run+1]));
    bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(bk_lengths);
    rank_i += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    if constexpr (check_head) {
        const uint32_t len = run >> sigma_bits;//get the length of the run where idx falls
        const bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
        const bool is_head = is_same_sym && (pf_sum-len)==idx_i;//check if idx is the head of the run
        rank_i -= (pf_sum-idx_i) * is_same_sym;
        rank_i = (rank_i<<1) | is_head;
        rank_i = (rank_i<<1) | is_same_sym;
    } else {
        rank_i -= (pf_sum-idx_i) * (last_symbol==symbol);
    }

    //=== deal with j
    bk_lengths = vshlq_u32(block, alpha_shift);
    //acc and block are still valid, no need to recompute them
    while(acc<=idx_j){
        bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(bk_lengths);
        rank_j += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        bk_lengths = vshlq_u32(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(bk_lengths);
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    idx_j-=prev_acc;

    //prefix sum
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 3), bk_lengths);
    bk_lengths = vaddq_u32(vextq_u32(vdupq_n_u32(0), bk_lengths, 2), bk_lengths);
    //

    idx_mask = vcgtq_u32(bk_lengths, vdupq_n_u32(idx_j));// mask for >idx
    res = vshrn_n_u32(idx_mask, 16);
    less_than = vget_lane_u64(vreinterpret_u64_u16(res), 0);
    idx_run = __builtin_ctzll(less_than)>>4;

    shuff = vaddq_u32(shuff_idxs, vdupq_n_u8(idx_run<<2));
    run_vec = vreinterpretq_u32_u8(vqtbl1q_u8(block, shuff));
    run = vgetq_lane_u32(run_vec, 0);
    last_symbol = run & alpha_m;

    pf_sum = vgetq_lane_u32(vreinterpretq_u32_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);

    bk_lengths = vshlq_u32(block, alpha_shift);
    bk_lengths = vandq_u32(bk_lengths, vld1q_u32(mask32x4[idx_run+1]));
    bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));

    tmp = vpaddlq_u32(bk_lengths);
    rank_j += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    rank_j -= (pf_sum-idx_j) * (last_symbol==symbol);

    return std::make_pair(rank_i, rank_j);
}

template<uint8_t bytes_per_run, bool check_head>
static inline std::pair<int64_t, int64_t> range_rank_neon_64x2(const uint8_t ** stream, const uint8_t sigma, uint64_t idx_i, uint64_t idx_j, uint8_t symbol){
    return {0,0};
}

static inline size_t first_run_neon_8x16(const uint8_t **stream, const uint8_t sigma, uint8_t symbol){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    uint8_t alpha_m = (1UL << sigma_bits)-1;
    const uint8x16_t alpha_mask = vdupq_n_u8(alpha_m);
    uint8x16_t sym_vec = vdupq_n_u8(symbol);

    uint8x16_t block = vld1q_u8(*stream);
    *stream+=16;
    uint8x16_t sym_mask = vceqq_u8(vandq_u8(block, alpha_mask), sym_vec);
    bool has_sym = vmaxvq_u8(sym_mask)==0xFF;
    size_t run_idx = 0;

    while(!has_sym){
        run_idx += 16;
        block = vld1q_u8(*stream);
        *stream += 16;
        sym_mask = vceqq_u8(vandq_u8(block, alpha_mask), sym_vec);
        has_sym = vmaxvq_u8(sym_mask)==0xFF;
    }

    const uint16x8_t sym_mask_16 = vreinterpretq_u16_u8(sym_mask);
    const uint8x8_t res = vshrn_n_u16(sym_mask_16, 4);
    const uint64_t matches = vget_lane_u64(vreinterpret_u64_u8(res), 0);
    uint8_t first = __builtin_ctzll(matches)>>2;

    return run_idx+first;
}

template<bool vbyte_compressed>
static inline size_t first_run_neon_16x8(const uint8_t **stream, const uint8_t sigma, uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    uint16_t alpha_m = (1UL << sigma_bits)-1;
    const uint16x8_t alpha_mask = vdupq_n_u16(alpha_m);
    const uint16x8_t sym_vec = vdupq_n_u16(symbol);

    uint16x8_t block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    uint16x8_t sym_mask = vceqq_u16(vandq_u16(block, alpha_mask), sym_vec);
    bool has_sym = vmaxvq_u16(sym_mask)==0xFFFF;
    size_t run_idx = 0;

    while(!has_sym){
        run_idx+=8;
        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        sym_mask = vceqq_u16(vandq_u16(block, alpha_mask), sym_vec);
        has_sym = vmaxvq_u16(sym_mask)==0xFFFF;
    }

    const uint8x8_t res = vshrn_n_u16(sym_mask, 4);
    const uint64_t matches = vget_lane_u64(vreinterpret_u64_u8(res), 0);
    uint8_t first = __builtin_ctzll(matches)>>3;
    return run_idx+first;
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline size_t first_run_neon_32x4(const uint8_t **stream, const uint8_t sigma, uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    uint32_t alpha_m = (1UL << sigma_bits)-1;
    const uint32x4_t alpha_mask = vdupq_n_u32(alpha_m);
    const uint32x4_t sym_vec = vdupq_n_u32(symbol);

    //const uint8_t *prev_state = *stream;

    uint32x4_t block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
    uint32x4_t sym_mask = vceqq_u32(vandq_u32(block, alpha_mask), sym_vec);
    bool has_sym = vmaxvq_u32(sym_mask)==0xFFFFFFFF;
    size_t run_idx = 0;

    while(!has_sym){
        run_idx+=4;
        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        sym_mask = vceqq_u32(vandq_u32(block, alpha_mask), sym_vec);
        has_sym = vmaxvq_u32(sym_mask)==0xFFFFFFFF;
    }

    const uint16x4_t res = vshrn_n_u32(sym_mask, 16);
    const uint64_t matches = vget_lane_u64(vreinterpret_u64_u16(res), 0);
    const uint8_t first= __builtin_ctzll(matches)>>4;
    return run_idx+first;
}

template<uint8_t bytes_per_run>
static inline size_t first_run_neon_64x2(const uint8_t **stream, const uint8_t sigma, uint8_t symbol){
    return 0;
}

//the function succ_neon_*** returns the run id for the leftmost run labeled "symbol" after position "idx".
//the output is a pair (run_id (int64_t), rank (uint64_t)), where "rank" is the number of times "symbol" occurs before "i"
//NOTE: run_id = -1 if there is no run labeled "symbol" from *i* onwards
//NOTE: if "idx" is the head of its run and that run is labeled "symbol", then run_id is the id for that run
template<bool overflow16, bool overflow32>
static inline std::pair<int64_t, uint64_t> succ_neon_8x16(const uint8_t **stream, const size_t n_runs,
                                                          const uint8_t sigma, uint64_t idx, uint8_t symbol){

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

    uint32_t rank=0;
    uint64_t run_id=0;
    while(acc<=idx){
        bk_lengths = vandq_u8(bk_lengths, vceqq_u8(vandq_u8(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        block =  vld1q_u8(*stream);
        *stream+=16;
        bk_lengths = vshlq_u8(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
        run_id+=16;
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
    bk_lengths = vshlq_u8(block, alpha_shift);
    const uint8x16_t mask_pref = vld1q_u8(mask8x16[idx_run+1]);
    bk_lengths = vandq_u8(bk_lengths, mask_pref);
    uint8x16_t sym_mask = vceqq_u8(vandq_u8(block, alpha_mask), sym_vec);
    bk_lengths = vandq_u8(bk_lengths, sym_mask);

    tmp = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    //check if the symbol of the run where idx falls matches the query symbol
    bool is_same_sym = last_symbol==symbol;
    rank -= (pf_sum-idx) * is_same_sym;

    //return if "symbol" matches the symbol where "i" lies
    int64_t opts[2] ={-1, (int64_t)run_id+idx_run};
    if(is_same_sym){
        bool is_head = (pf_sum-(run>>sigma_bits))==idx;//"i" is also the head position within the run
        return {opts[is_head], rank};
    }

    //the run where "idx" lies does not contain the successor information, we have scan forward to find it
    uint16x8_t mask_suff = vmvnq_u8(mask_pref);//only keep the runs after idx_run
    sym_mask = vandq_u8(sym_mask, mask_suff);
    bool has_sym = vmaxvq_u8(sym_mask)==0xFF;
    run_id+=16;

    while(run_id<n_runs && !has_sym){
        *stream+=16;
        block = vld1q_u8(*stream);
        sym_mask = vceqq_u8(vandq_u8(block, alpha_mask), sym_vec);
        has_sym = vmaxvq_u8(sym_mask)==0xFF;
        run_id+=16;
    }

    //get the leftmost occurrence of the symbol
    const uint16x8_t sym_mask_16 = vreinterpretq_u16_u8(sym_mask);
    const uint8x8_t res = vshrn_n_u16(sym_mask_16, 4);
    const uint64_t matches = vget_lane_u64(vreinterpret_u64_u8(res), 0);
    uint8_t first = __builtin_ctzll(matches)>>2;
    run_id = run_id - 16 + first;

    opts[1] = (int64_t)run_id;
    return {opts[has_sym && run_id<n_runs], rank};
}

template<bool vbyte_compressed, bool overflow8, bool overflow16>
static inline std::pair<uint64_t, uint64_t> succ_neon_16x8(const uint8_t **stream, const size_t n_runs,
                                                           const uint8_t sigma, uint64_t idx, uint8_t symbol){

    const uint8_t sigma_bits = sym_width(sigma);
    const int16x8_t alpha_shift = vdupq_n_u16(-sigma_bits);
    uint16_t alpha_m = (1UL << sigma_bits)-1;
    const uint16x8_t alpha_mask = vdupq_n_u16(alpha_m);
    const uint16x8_t sym_vec = vdupq_n_u16(symbol);

    const uint8_t *prev_state = *stream;

    uint16x8_t block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
    uint16x8_t bk_lengths = vshlq_u16(block, alpha_shift);

    uint64_t prev_acc=0;
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    uint64_t rank = 0;
    uint64_t run_id = 0;

    while(acc<=idx){
        bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
        run_id+=8;

        prev_state = *stream;
        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        bk_lengths = vshlq_u16(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
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
    bk_lengths = vshlq_u16(block, alpha_shift);
    const uint16x8_t mask_pref = vld1q_u16(mask16x8[idx_run+1]);
    bk_lengths = vandq_u16(bk_lengths, mask_pref);
    uint16x8_t sym_mask = vceqq_u16(vandq_u16(block, alpha_mask), sym_vec);
    bk_lengths = vandq_u16(bk_lengths, sym_mask);

    tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

    //check if the symbol of the run where idx falls matches the query symbol
    bool is_same_sym = last_symbol==symbol;
    rank -=(pf_sum-idx) * is_same_sym;

    int64_t opts[2] ={-1, (int64_t)run_id+idx_run};
    if(is_same_sym){
        uint16_t len = run >> sigma_bits;//get the length of the run where idx falls
        bool is_head = is_same_sym && (pf_sum-len)==idx;//check if idx is the head of the run
        return {opts[is_head], rank};
    }

    uint8x16_t mask_suff = vmvnq_u16(mask_pref);//only keep the runs after idx_run
    sym_mask = vandq_u16(sym_mask, mask_suff);
    bool has_sym = vmaxvq_u16(sym_mask)==0xFFFF;
    run_id+=8;

    while(run_id<n_runs && !has_sym){
        block = vreinterpretq_u16_u8(decode_block_neon<vbyte_compressed, 1, 2>(stream));
        sym_mask = vceqq_u16(vandq_u16(block, alpha_mask), sym_vec);
        has_sym = vmaxvq_u16(sym_mask)==0xFFFF;
        run_id+=8;
    }

    const uint8x8_t res = vshrn_n_u16(sym_mask, 4);
    const uint64_t matches = vget_lane_u64(vreinterpret_u64_u8(res), 0);
    uint8_t first = __builtin_ctzll(matches)>>3;
    run_id = run_id - 8 + first;
    opts[1] = (int64_t) run_id;

    return {opts[has_sym && run_id<n_runs], rank};
}

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint64_t> succ_neon_32x4(const uint8_t ** stream, const size_t n_runs,
                                                           const uint8_t sigma, uint64_t idx, uint8_t symbol){

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
    uint64_t rank = 0;
    uint64_t run_id=0;

    while(acc<=idx){
        bk_lengths = vandq_u32(bk_lengths, vceqq_u32(vandq_u32(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(bk_lengths);
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
        run_id+=4;

        prev_state = *stream;
        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        bk_lengths = vshlq_u32(block, alpha_shift);

        prev_acc = acc;
        tmp = vpaddlq_u32(bk_lengths);
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
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
    bk_lengths = vshlq_u32(block, alpha_shift);
    const uint32x4_t mask_pref = vld1q_u32(mask32x4[idx_run+1]);
    bk_lengths = vandq_u32(bk_lengths, mask_pref);
    uint32x4_t sym_mask = vceqq_u32(vandq_u32(block, alpha_mask), sym_vec);
    bk_lengths = vandq_u32(bk_lengths, sym_mask);

    tmp = vpaddlq_u32(bk_lengths);
    rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    bool is_same_sym = last_symbol==symbol;//check if the symbol of the run where idx falls matches the query symbol
    rank -= (pf_sum-idx) * is_same_sym;

    int64_t opts[2] = {-1, (int64_t)run_id+idx_run};
    if(is_same_sym){
        uint32_t len = run>>sigma_bits;
        bool is_head = is_same_sym && (pf_sum-len)==idx;//check if idx is the head of the run
        return {opts[is_head], rank};
    }

    uint32x4_t mask_suff =  vmvnq_u32(mask_pref);
    sym_mask = vandq_u32(sym_mask, mask_suff);
    bool has_sym = vmaxvq_u32(sym_mask)==0xFFFFFFFF;
    run_id+=4;

    while(run_id<n_runs && !has_sym){
        block = vreinterpretq_u32_u8(decode_block_neon<vbyte_compressed, 2, bytes_per_run>(stream));
        sym_mask = vceqq_u32(vandq_u32(block, alpha_mask), sym_vec);
        has_sym = vmaxvq_u32(sym_mask)==0xFFFFFFFF;
        run_id+=4;
    }

    const uint16x4_t res2 = vshrn_n_u32(sym_mask, 16);
    const uint64_t matches = vget_lane_u64(vreinterpret_u64_u16(res2), 0);
    uint8_t first = __builtin_ctzll(matches)>>4;
    run_id = run_id - 4 + first;

    opts[1] = (int64_t)run_id;
    return {opts[has_sym && run_id<n_runs], rank};
}

template<uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint64_t> succ_neon_64x2(const uint8_t ** stream, const size_t stream_bytes, const uint8_t sigma, uint64_t idx, uint8_t symbol){
    return {0, 0};
}
#endif //VLBT_SCAN_NEON_H
