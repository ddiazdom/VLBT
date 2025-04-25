//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef RLBWT_DYBL_H
#define RLBWT_DYBL_H

#include <cmath>
#include "bitstream.h"

template<bool vbyte_compressed>
static inline uint64_t inv_select_scalar(const uint16_t* stream, uint8_t sym, uint64_t idx){
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
#include "vbyte_simd_dec_tables.h"

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

static inline std::pair<uint64_t, uint8_t> inv_select_neon_8x16(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    //NOTE here I do not need to vbyte compress the block
    const uint8_t sigma_bits = sym_width(sigma);
    const int8x16_t alpha_shift = vdupq_n_u8(-sigma_bits);
    const uint8_t *stream_start = *stream;

    size_t l=0;
    uint8x16_t block =  vld1q_u8(*stream);
    *stream+=16;
    uint8x16_t bk_lengths = vshlq_u8(block, alpha_shift);
    uint64x2_t tmp  = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(bk_lengths)));
    uint64_t prev_acc=0, acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

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

    //prefix sum
    bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 15), bk_lengths);
    bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 14), bk_lengths);
    bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 12), bk_lengths);
    bk_lengths = vaddq_u8(vextq_u8(vdupq_n_u8(0), bk_lengths, 8), bk_lengths);
    //

    const uint8x16_t idx_mask = vcgtq_u8(bk_lengths, vdupq_n_u8(idx));//mask for >idx
    const uint8x8_t res = vshrn_n_u16(vreinterpretq_u16_u8(idx_mask), 4);
    const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
    uint8_t idx_run = __builtin_ctzll(less_than)>>2;

    uint8_t alpha_m = (1UL << sigma_bits)-1;
    const uint8x16_t alpha_mask = vdupq_n_u8(alpha_m);

    const uint8x16_t shuff = vdupq_n_u8(idx_run);
    const uint8x16_t run_vec = vqtbl1q_u8(block, shuff);
    const uint8x16_t sym_vec = vandq_u8(run_vec, alpha_mask);

    uint8_t run = vgetq_lane_u8(run_vec, 0);
    uint8_t pf_sum = vgetq_lane_u8(vqtbl1q_u8(bk_lengths, shuff), 0);
    pf_sum-= run>>sigma_bits;
    uint8_t sym = run & alpha_m;

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

template<bool vbyte_compressed>
static inline std::pair<uint64_t, uint8_t> inv_select_neon_16x8(const uint8_t **stream, uint8_t sigma, uint64_t idx){

    //TODO: assert idx fits 2 bytes
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

    //prefix sum
    bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 7), bk_lengths);
    bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 6), bk_lengths);
    bk_lengths = vaddq_u16(vextq_u16(vdupq_n_u16(0), bk_lengths, 4), bk_lengths);
    //

    const uint16x8_t idx_mask = vcgtq_u16(bk_lengths, vdupq_n_u16(idx));//mask for >idx
    const uint8x8_t res = vshrn_n_u16(idx_mask, 4);
    const uint64_t less_than = vget_lane_u64(vreinterpret_u64_u8(res), 0);
    uint8_t idx_run = __builtin_ctzll(less_than)>>3;
    uint16_t alpha_m = (1UL << sigma_bits)-1;
    const uint16x8_t alpha_mask = vdupq_n_u16(alpha_m);

    const uint8x16_t shuff_idxs = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    const uint8x16_t shuff = vaddq_u16(shuff_idxs, vdupq_n_u8(idx_run<<1));

    const uint16x8_t run_vec = vreinterpretq_u16_u8(vqtbl1q_u8(block, shuff));
    const uint16x8_t sym_vec = vandq_u16(run_vec, alpha_mask);

    uint16_t run = vgetq_lane_u16(run_vec, 0);
    uint16_t pf_sum = vgetq_lane_u16(vreinterpretq_u16_u8(vqtbl1q_u8(bk_lengths, shuff)), 0);
    pf_sum-= run>>sigma_bits;
    uint8_t sym = run & alpha_m;

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

    //TODO: assert idx fits 4 bytes
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

template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline std::pair<uint64_t, uint8_t> inv_select_neon_64x2(const uint8_t ** stream, uint8_t sigma, uint64_t idx){
    return {0,0};
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

typedef bitstream<size_t> stream_type;

template<size_t b_size, size_t b_runs, size_t s_factor>
struct rlbwt_vlb {

    static constexpr size_t block_size = b_size;
    static constexpr size_t scale_factor = s_factor;
    static constexpr size_t max_block_runs = b_runs;
    static constexpr uint8_t int_pt_width=7;//number of bits we use to encode the number of bits we use to encode pointers
    static constexpr uint8_t run_width = (sizeof(unsigned long)*8) - __builtin_clzl(b_runs-1);
    static constexpr uint64_t run_mask = (1UL<<run_width)-1UL;
    static constexpr uint8_t leaf_enc_width=4;

    uint64_t tot_syms=0;//total symbols in the text
    uint64_t orig_runs=0;//original number of runs in the BWT
    uint64_t eff_runs=0;//number of runs in the representation
    uint64_t max_freq=0;//frequency of the most frequent symbol
    uint64_t header_bytes=0;//bytes used for the header
    uint64_t lfs_bits=0;//number of bits at the beginning of the stream used by the ext succ/pred info of the low freq symbols.
    uint8_t sigma=0;//size of the effective text alphabet
    uint16_t ext_pt_width=0;//number of bits we use store the pointers to the trees
    uint8_t mtd_bits=0;//number of bits to encode the maximum tree distance
    std::vector<uint8_t> packed_alpha;//map the symbols from byte to eff alphabet
    std::vector<uint8_t> unpacked_alpha;//map eff alphabet to the original alphabet

    uint8_t levels=0;//maximum number of levels
    stream_type stream;//stream with the data

    rlbwt_vlb():levels(size_t(ceil(log(b_size)/log(s_factor)) - ceil(log(b_runs)/log(s_factor)))+1){
        //TODO static asserts in block_size, scale_factor, and b_runs
        // logarithm function to calculate value
        float lg = log(b_size) / log(s_factor);
        assert(lg==floor(lg));
        float lg2 = log(b_runs) / log(s_factor);
        assert(lg2==floor(lg2));
    }

    size_t serialize(std::ostream & ofs) const {
        size_t written_bytes = 0;
        written_bytes += serialize_elm(ofs, tot_syms);
        written_bytes += serialize_elm(ofs, orig_runs);
        written_bytes += serialize_elm(ofs, eff_runs);
        written_bytes += serialize_elm(ofs, max_freq);
        written_bytes += serialize_elm(ofs, header_bytes);
        written_bytes += serialize_elm(ofs, lfs_bits);
        written_bytes += serialize_elm(ofs, sigma);
        written_bytes += serialize_elm(ofs, ext_pt_width);
        written_bytes += serialize_elm(ofs, mtd_bits);

        written_bytes += serialize_plain_vector(ofs, packed_alpha);
        written_bytes += serialize_plain_vector(ofs, unpacked_alpha);

        written_bytes+=stream.serialize(ofs);
        return written_bytes;
    }

    void load(std::istream & ifs){
        load_elm(ifs, tot_syms);
        load_elm(ifs, orig_runs);
        load_elm(ifs, eff_runs);
        load_elm(ifs, max_freq);
        load_elm(ifs, header_bytes);
        load_elm(ifs, lfs_bits);
        load_elm(ifs, sigma);
        load_elm(ifs, ext_pt_width);
        load_elm(ifs, mtd_bits);
        load_plain_vector(ifs, packed_alpha);
        load_plain_vector(ifs, unpacked_alpha);
        stream.load(ifs);
    }

    [[nodiscard]] inline std::pair<uint64_t, uint8_t> inverse_select(size_t i) const {

        //initialize the block size
        size_t bk_sz = block_size;

        //get the block where index i lies
        uint64_t child = i/bk_sz;

        //skip the bits with the ext. succ/pred info of low-freq symbols
        size_t bit_pos = lfs_bits;

        //get the effective block where i lies and its byte offset within the stream
        //[p..p+ext_pt_width-1] is the area where the pointer information of bk lies in the stream
        size_t p = bit_pos+(ext_pt_width*child);
        p = stream.read(p, p+ext_pt_width-1);

        //std::cout<<p<<" "<<int(run_mask)<<" "<<int(run_width)<<" byte_pos"<<int(p>>run_width)<<" offset:"<<int(p & run_mask)<<std::endl;

        child = child-(p & run_mask);//eff child in the representation where i lies
        bit_pos =  (header_bytes + (p >> run_width))*8;//bit position where child begins in the stream

        uint8_t lvl=0;
        uint64_t sigma_pos[10]={0};
        uint64_t rank_pos[10]={0};
        uint8_t node_sigma[10]={0};
        uint8_t rank_width[10]={0};

        //skip ext succ/pred information
        size_t n_samps = stream.pop_count(bit_pos, bit_pos+2*sigma-1);
        bit_pos+=2*sigma;
        uint8_t w = stream.read(bit_pos, bit_pos+mtd_bits-1);//number bits we use to encode the tree distances for child
        bit_pos+=mtd_bits;
        bit_pos+=n_samps*w;//skip the n_samp tree distances
        bit_pos= INT_CEIL(bit_pos, 8)*8;//next byte-aligned position
        //
        node_sigma[lvl] = sigma;

        //read the node header
        lvl++;
        bool is_leaf = stream.read_bit(bit_pos++);
        node_sigma[lvl] = stream.pop_count(bit_pos, bit_pos+node_sigma[lvl-1]-1);
        sigma_pos[lvl]=bit_pos;
        bit_pos+=node_sigma[lvl-1];
        rank_pos[lvl] = bit_pos;
        rank_width[lvl] = sym_width(max_freq);
        bit_pos+=rank_width[lvl]*node_sigma[lvl];//skip rank information

        i-= child*bk_sz;//relative position of i within the child block

        while(!is_leaf){
            bk_sz/=scale_factor;
            child = i/bk_sz;
            assert(child<scale_factor);

            size_t child_info = stream.read(bit_pos, bit_pos+scale_factor-1);
            bit_pos+=scale_factor;

            size_t n_children = __builtin_popcount(child_info);//number of eff children
            assert(n_children>0);

            child_info &= (1<<(child+1))-1;//clean the bits marking the right siblings
            child = __builtin_popcount(child_info)-1;//eff child (zero-based)
            size_t n_real_lsib = 63-__builtin_clzll(child_info);//= select_1(child_info, (eff child)+1)-1
            i-=n_real_lsib*bk_sz;//number of symbols before child within the node

            bit_pos+=node_sigma[lvl]*n_children;//skip int succ/pred info

            //read how many bits we use to encode the pointers to the children
            size_t p_width = stream.read(bit_pos, bit_pos+int_pt_width-1);
            bit_pos+=int_pt_width;

            p = bit_pos+(child*p_width);
            p = stream.read(p, p+p_width-1);

            //skip the pointer to the children and position the bit in the next byte-aligned position
            bit_pos = INT_CEIL((bit_pos+(n_children*p_width)), 8)*8;
            //add the bit offset. Now bit_pos points to child
            bit_pos+= p*8;

            //start reading the header of child (there is no ext succ/pred info)
            lvl++;
            is_leaf = stream.read_bit(bit_pos++);
            node_sigma[lvl] = stream.pop_count(bit_pos, bit_pos+node_sigma[lvl-1]-1);
            sigma_pos[lvl]=bit_pos;
            bit_pos+=node_sigma[lvl-1];
            rank_pos[lvl] = bit_pos;
            rank_width[lvl] = sym_width(bk_sz*scale_factor);
            bit_pos+=rank_width[lvl]*node_sigma[lvl];//skip rank information
        }

        uint8_t leaf_enc = stream.read(bit_pos, bit_pos+leaf_enc_width-1);
        bit_pos+= leaf_enc_width;
        auto const *runs_stream = ((uint8_t *)stream.stream)+(INT_CEIL(bit_pos, 8));

        std::pair<uint64_t, uint8_t> rank_answer;

        //scan the runs in the leaf according to the leaf encoding
        switch (leaf_enc) {
            case 1:
                rank_answer = inv_select_neon_8x16(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 2:
                rank_answer = inv_select_neon_16x8<false>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 2 bytes (no vbyte)
                break;
            case 3:
                rank_answer = inv_select_neon_32x4<false, 3>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 3 bytes (no vbyte)
                break;
            case 4:
                rank_answer = inv_select_neon_32x4<false, 4>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 4 bytes (no vbyte)
                break;
            case 5:
                rank_answer = inv_select_neon_64x2<false, 5>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 5 bytes (no vbyte)
                break;
            case 7:
                rank_answer = inv_select_neon_16x8<true>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 2 bytes (vbyte)
                break;
            case 8:
                rank_answer = inv_select_neon_32x4<true, 3>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 3 bytes (vbyte)
                break;
            case 9:
                rank_answer = inv_select_neon_32x4<true, 4>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 10:
                rank_answer = inv_select_neon_64x2<true, 5>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 5 bytes (vbyte)
                break;
            default:
                std::cout<<"Undefined encoding"<<std::endl;
                exit(1);
        }

        bit_pos = rank_pos[lvl]+(rank_answer.second*rank_width[lvl]);

        rank_answer.first += stream.read(bit_pos, bit_pos+rank_width[lvl]-1);//add rank information
        rank_answer.second = stream.select(sigma_pos[lvl], sigma_pos[lvl]+node_sigma[lvl-1]-1, rank_answer.second+1);//update the symbol
        lvl--;

        //std::cout<<symbol<<std::endl;
        while(lvl>0){
            bit_pos = rank_pos[lvl]+(rank_answer.second*rank_width[lvl]);
            rank_answer.first += stream.read(bit_pos, bit_pos+rank_width[lvl]-1);
            rank_answer.second = stream.select(sigma_pos[lvl], sigma_pos[lvl]+node_sigma[lvl-1]-1, rank_answer.second+1);
            //std::cout<<stream.pop_count(sigma_pos[lvl], sigma_pos[lvl]+node_sigma[lvl-1]-1)<<" "<<symbol<<std::endl;
            lvl--;
        }

        return rank_answer;
    }

    [[nodiscard]] inline uint8_t eff2byte(uint8_t eff_sym) const {
        assert(eff_sym<sigma);
        return unpacked_alpha[eff_sym];
    }

    [[nodiscard]] inline uint64_t size() const {
        return tot_syms;
    }
};

#endif //RLBWT_DYBL_H