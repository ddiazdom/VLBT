//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef RLBWT_DYBL_H
#define RLBWT_DYBL_H

#include <cmath>
#include "bitstream.h"

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
#include "vbyte_simd_dec_tables.h"

static inline uint64_t rank8x16_neon(const uint8_t * stream, uint8_t sigma, uint8_t sym, uint64_t idx){
   return 0;
}

#ifdef __aarch64__
typedef uint8x16_t decode_t;
#else
typedef uint8x8x2_t decode_t;
#endif

static inline decode_t  _decode_vbyte_neon(const uint8_t key, const uint8_t *dataPtrPtr) {

    uint8_t len;
    uint8_t *pshuf = (uint8_t *)&shuffleTable[key];
    uint8x16_t decodingShuffle = vld1q_u8(pshuf);

    uint8x16_t compressed = vld1q_u8(dataPtrPtr);
#ifdef AVOIDLENGTHLOOKUP
// this avoids the dependency on lengthTable, see https://github.com/lemire/streamvbyte/issues/12
    len = pshuf[12 + (key >> 6)] + 1;
#else
    len = lengthTable[key];
#endif

#ifdef __aarch64__
    uint8x16_t data = vqtbl1q_u8(compressed, decodingShuffle);
#else
    uint8x8x2_t codehalves = {{vget_low_u8(compressed),
                               vget_high_u8(compressed)}};
    uint8x8x2_t data = {{vtbl2_u8(codehalves, vget_low_u8(decodingShuffle)),
                         vtbl2_u8(codehalves, vget_high_u8(decodingShuffle))}};
#endif
    dataPtrPtr += len;
    return data;
}


static inline uint64_t rank16x8_neon(const uint16_t * stream, uint8_t sigma, uint8_t sym, uint64_t idx){

    uint8_t sigma_bits = sym_width(sigma);
    uint8_t sigma_mask = (1UL << sigma_bits)-1;

    const uint16x8_t sym_vec = vdupq_n_u16(sym);
    const uint16x8_t alpha_mask = vdupq_n_u16(sigma_mask);
    const uint16x8_t alpha_shift = vdupq_n_u16(-sigma_bits);

    size_t i=0;

    uint16x8_t block = vld1q_u16(&stream[0]);
    uint16x8_t bk_lengths = vshlq_u16(block, alpha_shift);
    uint64x2_t tmp   = vpaddlq_u32(vpaddlq_u16(bk_lengths));
    uint64_t acc = vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0], rank=0;

    while(acc<idx){
        bk_lengths = vandq_u16(bk_lengths, vceqq_u16(vandq_u16(block, alpha_mask), sym_vec));
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        rank += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];

        i+=8;
        block = vld1q_u16(&stream[i]);
        bk_lengths = vshlq_u16(block, alpha_shift);
        tmp = vpaddlq_u32(vpaddlq_u16(bk_lengths));
        acc += vadd_u64(vget_high_u64(tmp), vget_low_u64(tmp))[0];
    }

    return rank;
}

static inline uint64_t rank4x32_neon(const uint32_t * stream, uint8_t sigma, uint8_t sym, uint64_t idx){
    return 0;
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
        stream.load(ifs);
    }

    inline std::pair<size_t, uint8_t> inverse_select(size_t i) const {

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

        //skip ext succ/pred information
        size_t n_samps = stream.pop_count(bit_pos, bit_pos+2*sigma-1);
        bit_pos+=2*sigma;
        uint8_t w = stream.read(bit_pos, bit_pos+mtd_bits-1);//number bits we use to encode the tree distances for child
        bit_pos+=mtd_bits;
        bit_pos+=n_samps*w;//skip the n_samp tree distances
        bit_pos= INT_CEIL(bit_pos, 8)*8;//next byte-aligned position
        //

        //read the node header
        bool is_leaf = stream.read_bit(bit_pos++);
        size_t parent_sigma = sigma;
        size_t node_sigma = stream.pop_count(bit_pos, bit_pos+parent_sigma-1);
        bit_pos+=parent_sigma;
        bit_pos+=sym_width(max_freq)*node_sigma;//skip rank information
        i-= child*bk_sz;//relative position of i within the child block

        while(!is_leaf){
            bk_sz/=scale_factor;
            child = i/bk_sz;
            assert(child<scale_factor);

            size_t child_info = stream.read(bit_pos, bit_pos+scale_factor-1);
            bit_pos+=scale_factor;

            size_t n_children = __builtin_popcount(child_info);//number of eff children

            child_info &= (1<<(child+1))-1;//clean the bits marking the right siblings
            child = __builtin_popcount(child_info)-1;//eff child (zero-based)
            size_t n_real_lsib = 63-__builtin_clzll(child_info);//= select_1(child_info, (eff child)+1)-1
            i-=n_real_lsib*bk_sz;//number of symbols before child within the node

            bit_pos+=node_sigma*n_children;//skip int succ/pred info

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
            is_leaf = stream.read_bit(bit_pos++);
            parent_sigma = node_sigma;
            node_sigma = stream.pop_count(bit_pos, bit_pos+parent_sigma-1);
            bit_pos+=parent_sigma;
            bit_pos+=sym_width(bk_sz)*node_sigma;//skip rank information
        }

        uint8_t leaf_enc = stream.read(bit_pos, bit_pos+leaf_enc_width-1);

        bit_pos+= leaf_enc_width;
        bit_pos = INT_CEIL(bit_pos, 8);

        return {0,0};
    }
};

#endif //RLBWT_DYBL_H