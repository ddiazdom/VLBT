//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef RLBWT_DYBL_H
#define RLBWT_DYBL_H

#include <cmath>
#include "bitstream.h"

template<bool vbyte_compressed>
static inline uint64_t inv_select_scl_8(const uint16_t* stream, uint8_t sym, uint64_t idx){
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

template<bool vbyte_compressed>
static inline uint64_t inv_select_scl_16(const uint16_t* stream, uint8_t sym, uint64_t idx){
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

template<bool vbyte_compressed>
static inline uint64_t inv_select_scl_32(const uint16_t* stream, uint8_t sym, uint64_t idx){
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

template<bool vbyte_compressed>
static inline uint64_t inv_select_scl_64(const uint16_t* stream, uint8_t sym, uint64_t idx){
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

#if defined(__ARM_NEON__)
#include "neon_scan.h"
#define INV_SELECT_8 inv_select_neon_8x16
#define INV_SELECT_16 inv_select_neon_16x8
#define INV_SELECT_32 inv_select_neon_32x4
#define INV_SELECT_64 inv_select_neon_64x2
#elif defined(__AVX2__)
#include "avx2_scan.h"
#define INV_SELECT_8 inv_select_avx2_8x16
#define INV_SELECT_16 inv_select_avx2_16x8
#define INV_SELECT_32 inv_select_avx2_32x4
#define INV_SELECT_64 inv_select_avx2_64x2
#elif defined(__SSE4_2__)
#include "sse42_scan.h
#define INV_SELECT_8 inv_select_sse42_8x16
#define INV_SELECT_16 inv_select_sse42_16x8
#define INV_SELECT_32 inv_select_sse42_32x4
#define INV_SELECT_64 inv_select_sse42_64x2
#else
#define INV_SELECT_8 inv_select_scl_8
#define INV_SELECT_16 inv_select_scl_16
#define INV_SELECT_32 inv_select_scl_32
#define INV_SELECT_64 inv_select_scl_64
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

    static inline void find_leaf(size_t i){

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
                rank_answer = INV_SELECT_8(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 2:
                rank_answer = INV_SELECT_16<false>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 2 bytes (no vbyte)
                break;
            case 3:
                rank_answer = INV_SELECT_32<false, 3>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 3 bytes (no vbyte)
                break;
            case 4:
                rank_answer = INV_SELECT_32<false, 4>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 4 bytes (no vbyte)
                break;
            case 5:
                rank_answer = INV_SELECT_64<false, 5>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 5 bytes (no vbyte)
                break;
            case 7:
                rank_answer = INV_SELECT_16<true>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 2 bytes (vbyte)
                break;
            case 8:
                rank_answer = INV_SELECT_32<true, 3>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 3 bytes (vbyte)
                break;
            case 9:
                rank_answer = INV_SELECT_32<true, 4>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 10:
                rank_answer = INV_SELECT_64<true, 5>(reinterpret_cast<const uint8_t **>(&runs_stream), node_sigma[lvl], i);//runs use 5 bytes (vbyte)
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

    [[nodiscard]] inline uint8_t operator[](size_t i) const {
        assert(i<tot_syms);
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