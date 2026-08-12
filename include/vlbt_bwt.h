/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */

#ifndef VLBT_BWT_TH_H
#define VLBT_BWT_TH_H

#include <cmath>
#include <vector>
#include "bit_stream.h"
#include "def_scan.h"
#include "vlbt_common.h"
#include "logger.h"

template<VLBT_TYPE var, size_t b_size, size_t b_runs=64, size_t s_factor=4>
class vlbt_bwt  {

public:
    static_assert(var==RLBWT || var==RLBWT_WITH_TOEHOLDS);

    static constexpr size_t block_size = b_size;
    static constexpr size_t scale_factor = s_factor;
    static constexpr size_t max_block_runs = b_runs;

    //the number of bits we use to encode the number of bits we use to encode pointers
    // 6 bits indicate that we use 2^6-1 = 63 bits to encode pointers, which is enough for any possible known application
    static constexpr uint8_t int_pt_width = 6;
    static constexpr uint8_t run_width = (sizeof(unsigned long)*8) - __builtin_clzl(b_runs-1);
    static constexpr uint8_t leaf_enc_width = 4;
    static constexpr VLBT_TYPE tag = var;

    //the number of bits to encode the number of bytes the run sequence uses:
    // number of bits we use to encode the number of bytes that the runs use in a leaf
    static constexpr uint8_t run_byte_w = ((sizeof(unsigned long)*8) - __builtin_clzl((b_runs*8) + (b_runs/8)))+1;
    typedef bit_stream<size_t> stream_type;

    struct tree_path_type{
        uint8_t lvl=0;
        uint64_t bit_pos=0;
        uint64_t sigma_pos[10]={0};
        uint64_t rank_pos[10]={0};
        uint8_t node_sigma[10]={0};
        uint8_t rank_width[10]={0};
        uint8_t leaf_enc=0;
    };

    struct inv_sel_sa_ans{
        uint8_t sym;
        int64_t rank;
        int64_t sa_samp;
    };

    struct succ_info{
        size_t bit_pos=0xffffffffffffffff;
        uint64_t rank=0;
        uint8_t symbol=0;
        uint8_t r_width=0;
        uint8_t node_sigma=0;
        size_t bk_sz=0;
    };

    uint64_t tot_syms=0;//total symbols in the text
    uint64_t orig_runs=0;//original number of runs in the BWT
    uint64_t eff_runs=0;//number of runs in the representation
    uint64_t max_freq=0;//frequency of the most frequent symbol
    uint64_t header_bytes=0;//bytes used for the header
    uint64_t lfs_bits=0;//number of bits at the beginning of the stream used by the ext succ/pred info of the low-freq symbols.
    uint8_t sigma=0;//size of the effective text alphabet
    uint16_t ext_pt_width=0;//number of bits we use store the pointers to the trees
    uint8_t mtd_bits=0;//number of bits to encode the maximum tree distance
    std::vector<uint64_t> C;//F column in the FM index
    std::vector<uint8_t> packed_alpha;//map the symbols from byte to eff alphabet
    std::vector<uint8_t> unpacked_alpha;//map eff alphabet to the original alphabet
    size_t subsamp_step=0;//subsampling parameter

    uint8_t levels=0;//maximum number of levels
    stream_type stream;//stream with the data

    vlbt_bwt():
        levels(static_cast<size_t>(ceil(LOG_BASE(s_factor, b_size)) - ceil(LOG_BASE(s_factor, b_runs))) + 1){
        // static asserts in block_size, scale_factor, and b_runs
        assert(LOG_BASE(s_factor, b_size)==floor(LOG_BASE(s_factor, b_size)));
        assert(LOG_BASE(s_factor, b_runs)==floor(LOG_BASE(s_factor, b_runs)));
    }

    const std::vector<uint8_t>& get_packed_alpha(){
        return packed_alpha;
    }

    const std::vector<uint8_t>& get_unpacked_alpha(){
        return unpacked_alpha;
    }

    size_t serialize(std::ostream & ofs) const {
        size_t written_bytes = 0;

        //this is to check the template arguments once the data structure is loaded
        written_bytes += serialize_elm(ofs, tag);
        written_bytes += serialize_elm(ofs, block_size);
        //

        written_bytes += serialize_elm(ofs, tot_syms);
        written_bytes += serialize_elm(ofs, orig_runs);
        written_bytes += serialize_elm(ofs, eff_runs);
        written_bytes += serialize_elm(ofs, max_freq);
        written_bytes += serialize_elm(ofs, header_bytes);
        written_bytes += serialize_elm(ofs, lfs_bits);
        written_bytes += serialize_elm(ofs, sigma);
        written_bytes += serialize_elm(ofs, ext_pt_width);
        written_bytes += serialize_elm(ofs, mtd_bits);
        if constexpr (tag==RLBWT_WITH_TOEHOLDS) {
            written_bytes += serialize_elm(ofs, subsamp_step);
        }

        written_bytes += serialize_plain_vector(ofs, C);
        written_bytes += serialize_plain_vector(ofs, packed_alpha);
        written_bytes += serialize_plain_vector(ofs, unpacked_alpha);

        written_bytes+=stream.serialize(ofs);
        return written_bytes;
    }

    void load(std::istream & ifs){

        VLBT_TYPE tag_tmp;
        load_elm(ifs, tag_tmp);
        if (tag_tmp!=tag) {
            LOG_ERROR("Index type mismatch: the file holds type "+std::to_string(tag_tmp)+
                      ", but type "+std::to_string(tag)+" was requested");
            exit(1);
        }
        size_t b_size_tmp;
        load_elm(ifs, b_size_tmp);
        if (b_size_tmp!=block_size) {
            LOG_ERROR("Block size mismatch: the index was built with "+std::to_string(b_size_tmp)+
                      ", but "+std::to_string(block_size)+" was requested");
            exit(1);
        }

        load_elm(ifs, tot_syms);
        load_elm(ifs, orig_runs);
        load_elm(ifs, eff_runs);
        load_elm(ifs, max_freq);
        load_elm(ifs, header_bytes);
        load_elm(ifs, lfs_bits);
        load_elm(ifs, sigma);
        load_elm(ifs, ext_pt_width);
        load_elm(ifs, mtd_bits);
        if constexpr (tag==RLBWT_WITH_TOEHOLDS) {
            load_elm(ifs, subsamp_step);
        }

        if(!ifs){
            LOG_ERROR("The index file is truncated or corrupt");
            exit(1);
        }

        load_plain_vector(ifs, C);
        load_plain_vector(ifs, packed_alpha);
        load_plain_vector(ifs, unpacked_alpha);

        stream.load(ifs);

        if(!ifs){
            LOG_ERROR("The index file is truncated or corrupt");
            exit(1);
        }
    }

    void find_prev_bit_pos(uint64_t& child_i, uint64_t& child_j, size_t& bit_pos_i, size_t& bit_pos_j) const {
        //get the effective block where "i" lies and its byte offset within the stream
        //[p..p+ext_pt_width-1] is the area where the pointer information of bk lies in the stream

        const size_t loc_lfs_bits = lfs_bits;
        const size_t loc_ext_pt_width = ext_pt_width;

        size_t p = loc_lfs_bits + loc_ext_pt_width*child_i;
        p = stream.read(p, p+loc_ext_pt_width-1);
        child_i -= p >> (run_width+1) & -(p & 1);
        p = loc_lfs_bits + loc_ext_pt_width*child_i;
        bit_pos_i = (header_bytes + (stream.read(p, p+loc_ext_pt_width-1)>>1)) * 8;

        p = loc_lfs_bits + loc_ext_pt_width*child_j;
        p = stream.read(p, p+loc_ext_pt_width-1);
        child_j -= p >> (run_width+1) & -(p & 1);
        p = loc_lfs_bits + loc_ext_pt_width*child_j;
        bit_pos_j = (header_bytes + (stream.read(p, p+loc_ext_pt_width-1)>>1)) * 8;
    }

    size_t find_prev(uint64_t& child) const {
        //get the effective block where "i" lies and its byte offset within the stream
        //[p..p+ext_pt_width-1] is the area where the pointer information of bk lies in the stream

        const size_t loc_lfs_bits = lfs_bits;
        const size_t loc_ext_pt_width = ext_pt_width;

        size_t p = loc_lfs_bits + (loc_ext_pt_width*child);
        p = stream.read(p, p+loc_ext_pt_width-1);
        //size_t offset = (p >> (run_width+1)) * ((p & 1)>0);
        //child -= offset;//eff child in the representation where "i" lies
        child -= (p >> (run_width+1)) & -(p & 1);
        p = loc_lfs_bits + loc_ext_pt_width*child;
        return (header_bytes + (stream.read(p, p+loc_ext_pt_width-1)>>1)) * 8;//bit-position where "child" begins in the stream
    }

    size_t find_next(uint64_t& child) const {

        const size_t loc_lfs_bits = lfs_bits;
        const size_t loc_ext_pt_width = ext_pt_width;

        size_t p = loc_lfs_bits + (loc_ext_pt_width*child);
        p = stream.read(p, p+loc_ext_pt_width-1);
        //size_t offset = ((p >> 1) & ((1<<run_width)-1)) * ((p & 1)>0);
        //child += offset;//eff child in the representation where "i" lies
        child += ((p >> 1) & ((1<<run_width)-1)) & -(p & 1);
        p = loc_lfs_bits + (loc_ext_pt_width * child);
        return (header_bytes + (stream.read(p, p+loc_ext_pt_width - 1) >> 1)) * 8;
    }

    void skip_ext_succ_info(size_t& bit_pos) const {
        //skip ext succ/pred information
        const size_t n_samps = stream.pop_count(bit_pos, bit_pos+sigma-1);
        bit_pos+=sigma;
        const uint8_t w = stream.read(bit_pos, bit_pos+mtd_bits-1);//number bits we use to encode the tree distances for "child"
        bit_pos+=mtd_bits;
        bit_pos+=n_samps*w;//skip the n_samp tree distances
        bit_pos= INT_CEIL(bit_pos, 8)*8;//next byte-aligned position (trees are byte-aligned)
        //
    }

    void decode_ext_succ_info(size_t& bit_pos, size_t child, uint8_t symbol) const {
        //the position where the offset for the successor is located
        const size_t succ_pos = stream.pop_count(bit_pos, bit_pos+symbol)-1;//works only because bit_stream[bit_pos+symbol] is true
        bit_pos+=sigma;
        const uint8_t w = stream.read(bit_pos, bit_pos+mtd_bits-1);//number bits we use to encode the tree distances for "child"
        bit_pos+=mtd_bits + (succ_pos*w);
        const size_t succ_child = child + stream.read(bit_pos, bit_pos+w-1);//read the offset of the successor
        size_t p = lfs_bits + (ext_pt_width*succ_child);
        p = stream.read(p, p+ext_pt_width-1)>>1;
        bit_pos = (header_bytes+p)*8;
    }

    [[nodiscard]] size_t find_low_freq_succ(size_t i, uint8_t symbol) const {

        symbol = stream.pop_count(0, symbol)-1;//this works because stream[symbol] is true
        const size_t c_bits = sigma;//c_bits + (n_symbol+1)*40 contains pointers to the areas where the info lies
        const uint8_t w = sym_width(INT_CEIL(tot_syms, block_size)*block_size);

        //read the area of the stream where the info of "symbol" lies
        size_t ptr = c_bits + symbol*40;
        int64_t first = stream.read(ptr, ptr+39);
        ptr+=40;
        int64_t last = stream.read(ptr, ptr+39)-w;
        size_t n = ((last-first)/w)+1;

        size_t idx;
        while(first<=last){
            const int64_t mid = first + static_cast<int64_t>((n / 2) * w);
            idx = stream.read(mid, mid+w-1);
            if(idx<i){
                first = mid+w;
            }else{
                last = mid-w;
            }
            n = ((last-first)/w)+1;
        }

        last+=w;
        idx = stream.read(last, last+w-1);
        uint64_t child = idx/block_size;
        return find_next(child);
    }

    template<bool check_head>
    [[nodiscard]] auto scan_leaf(uint64_t bit_pos, size_t i, size_t j, uint8_t symbol,
                                 const size_t node_sigma, const size_t rank_width) const {

        const uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
        symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//only works because bit_stream[bit_pos+symbol] is true
        bit_pos+=node_sigma;
        const size_t r_pos = bit_pos + symbol*rank_width;
        const uint64_t rank=stream.read(r_pos, r_pos+rank_width-1);
        bit_pos+=new_sigma*rank_width;

        const uint8_t leaf_enc = stream.read(bit_pos, bit_pos+leaf_enc_width-1);
        bit_pos+= leaf_enc_width;

        bool false_break;
        if constexpr (var==RLBWT_WITH_TOEHOLDS){
            if constexpr (check_head){
                bit_pos+= run_byte_w + run_width;//skip the number of bytes we use to encode the runs
                false_break = stream.read_bit(bit_pos++);//read if the leftmost run in this head is artificial
            }else{
                bit_pos+= run_byte_w + run_width + 1;//skip the metadata of the runs
            }
        }

        const uint8_t *leaf_addr = reinterpret_cast<uint8_t *>(stream.stream)+(INT_CEIL(bit_pos, 8));
        std::pair<uint64_t, uint64_t> ans;

        //LEAF encoding (bpr=bytes per run):
        //0: 1 bpr, no overflow
        //1: 1 bpr, overflow of 32 elements but not 16
        //2: 1 brp, overflow of 16 and 32 elements

        //3: 2 bpr, no_vbyte, no overflow
        //4: 2 bpr, no_vbyte, overflow 16 elements but not 8
        //5: 2 bpr, no_vbyte, overflow 8 and 16 elements
        //6: 2 bpr, vbyte, no overflow
        //7: 2 bpr, vbyte, overflow of 16 elements but not 8
        //8: 2 bpr, vbyte, overflow of 8 and 16 elements

        //9: 3 bpr, no_vbyte
        //10: 3 bpr, vbyte
        //11: 4 bpr, no_vbyte
        //12: 4 bpr, vbyte

        //13: 5 bpr, vbyte
        //14: 6 bpr, vbyte
        //15: 7 bpr, vbyte

        //scan the runs in the leaf according to the leaf encoding
        switch(leaf_enc) {
            case 0:
                ans = RANGE_RANK_8<false, false, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 1 byte (no vbyte)
                break;
            case 1:
                ans = RANGE_RANK_8<false, true, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 1 byte (no vbyte)
                break;
            case 2:
                ans = RANGE_RANK_8<true, true, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 2 bytes (no vbyte)
                break;

            case 3://template param: vbyte?, overflow8?, overflow16?
                ans = RANGE_RANK_16<false, false, false, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 3 bytes (no vbyte)
                break;
            case 4:
                ans = RANGE_RANK_16<false, false, true, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 4 bytes (no vbyte)
                break;
            case 5:
                ans = RANGE_RANK_16<false, true, true, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 5 bytes (no vbyte)
                break;
            case 6:
                ans = RANGE_RANK_16<true, false, false, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 2 bytes (vbyte)
                break;
            case 7:
                ans = RANGE_RANK_16<true, false, true, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 3 bytes (vbyte)
                break;
            case 8:
                ans = RANGE_RANK_16<true, true, true, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 4 bytes (vbyte)
                break;

            case 9://template param: vbyte?, bpr
                ans = RANGE_RANK_32<false,3, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 4 bytes (vbyte)
                break;
            case 10:
                ans = RANGE_RANK_32<true,3, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 4 bytes (vbyte)
                break;
            case 11:
                ans = RANGE_RANK_32<false,4, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 4 bytes (vbyte)
                break;
            case 12:
                ans = RANGE_RANK_32<true,4, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 4 bytes (vbyte)
                break;

            case 13://template param: bpr
                ans = RANGE_RANK_64<5, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 5 bytes (vbyte)
                break;
            case 14:
                ans = RANGE_RANK_64<6, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 5 bytes (vbyte)
                break;
            case 15:
                ans = RANGE_RANK_64<7, check_head>(&leaf_addr, new_sigma, i, j, symbol);//runs use 5 bytes (vbyte)
                break;
            default:
                LOG_ERROR("Undefined leaf encoding");
                exit(1);
        }

        if constexpr (check_head){
            //complete information
            uint64_t i_rank = ans.first;
            const bool is_same_sym = i_rank & 1UL;//the query symbol is the same as the symbol of the run where "i" falls
            i_rank>>=1;
            const bool is_head = i_rank & 1UL;//the query symbol is the same as the symbol of the run where "i" falls and "i" is the head of that run
            i_rank>>=1;
            //
            return std::make_tuple(rank + i_rank, rank + ans.second,  (!is_same_sym || (is_head && i>0) || (is_head && !false_break)));
        } else {
            return std::make_pair<uint64_t, uint64_t>(rank+ans.first, rank+ans.second);
        }
    }

    [[nodiscard]] std::pair<uint64_t, uint64_t> range_rank_no_sa_head(size_t i, size_t j, uint8_t symbol) const {

        symbol = packed_alpha[symbol];

        size_t bit_pos_i, bit_pos_j;
        uint64_t child_i = i/block_size, child_j = j/block_size;
        find_prev_bit_pos(child_i, child_j, bit_pos_i, bit_pos_j);

        uint64_t rank_i = 0, rank_j = 0;
        uint8_t rank_width = sym_width(max_freq);

        //i and j are on different trees, meaning there is no advantage, and we proceed as usual
        if(child_i != child_j) {
            const size_t prev_bit_pos_i = bit_pos_i;
            skip_ext_succ_info(bit_pos_i);
            bool is_leaf = stream.read_bit(bit_pos_i++);
            bool has_symbol = stream.read_bit(bit_pos_i+symbol);

            if(!has_symbol){
                rank_i = read_rank_from_succ(prev_bit_pos_i, i, child_i, symbol);
            }else {
                rank_i = subtree_rank(bit_pos_i, i-child_i*block_size, symbol, sigma, rank_width, block_size, is_leaf, false);
            }

            const size_t prev_bit_pos_j = bit_pos_j;
            skip_ext_succ_info(bit_pos_j);
            is_leaf = stream.read_bit(bit_pos_j++);
            has_symbol = stream.read_bit(bit_pos_j+symbol);

            if (!has_symbol) {
                rank_j = read_rank_from_succ(prev_bit_pos_j, j, child_j, symbol);
            }else {
                rank_j = subtree_rank(bit_pos_j, j-child_j*block_size, symbol, sigma, rank_width, block_size, is_leaf, false);
            }
#ifdef VLBT_TRACE_RANK
            fprintf(stderr, "[trace] different blocks: block_i=%llu (offset %llu) block_j=%llu (offset %llu) "
                            "j_has_symbol=%d -> rank_i=%llu rank_j=%llu\n",
                    (unsigned long long)child_i, (unsigned long long)(i-child_i*block_size),
                    (unsigned long long)child_j, (unsigned long long)(j-child_j*block_size),
                    (int)has_symbol, (unsigned long long)rank_i, (unsigned long long)rank_j);
#endif
            return std::make_pair(rank_i, rank_j);
        }

        size_t bit_pos = bit_pos_i;
        skip_ext_succ_info(bit_pos);

        //positions i and j fall in the same block, but the block does not contain the symbol
        if(!stream.read_bit(bit_pos+1+symbol)){
            return std::make_pair(0, 0);
        }

        //positions i and j fall in the same block and block has the symbol, but the block is a leaf
        size_t offset = child_i*block_size;
        if(stream.read_bit(bit_pos++)){
            return scan_leaf<false>(bit_pos, i-offset, j-offset, symbol, sigma, rank_width);
        }
        //

        //positions i and j fall within the same internal node
        int64_t rank = 0;
        uint8_t node_sigma = sigma;
        size_t bk_sz = block_size;

        bool traverse_common_path = false;
        size_t child_info, succ_pred_info, pos, p_width, parent_ptr_area;
        i-=offset;
        j-=offset;

#ifdef VLBT_TRACE_RANK
        //compile with -DVLBT_TRACE_RANK to dump this descent. The block disappears
        //at preprocessing time otherwise, so normal builds are unaffected
        fprintf(stderr, "[trace] descent starts: packed_sym=%u block=%llu i=%llu j=%llu (block-relative)\n",
                (unsigned)symbol, (unsigned long long)(offset/block_size),
                (unsigned long long)i, (unsigned long long)j);
#endif

        do {
            const uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//this works only because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            const size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            node_sigma=new_sigma;

            child_info = stream.read(bit_pos, bit_pos+scale_factor-1);//children info
            bit_pos += scale_factor;

            const size_t n_children = __builtin_popcount(child_info);//number of eff children
            assert(n_children>0);
            bk_sz/=scale_factor;
            rank_width = sym_width(bk_sz*scale_factor);

            child_i = i/bk_sz;
            child_j = j/bk_sz;
            assert(child_i<scale_factor && child_j<scale_factor);

            //read the effective child for i
            size_t eff_c_info = child_info & ((1<<(child_i+1))-1);//clean the bits marking the right siblings
            child_i = __builtin_popcount(eff_c_info)-1;//eff child (zero-based)
            size_t n_real_lsib = 63-__builtin_clzll(eff_c_info);//= select_1(child_info, (eff child)+1)-1
            i-=n_real_lsib*bk_sz;//number of symbols before child i within the node

            //read the effective child for j
            eff_c_info = child_info & ((1<<(child_j+1))-1);//clean the bits marking the right siblings
            child_j = __builtin_popcount(eff_c_info)-1;//eff child (zero-based)
            n_real_lsib = 63-__builtin_clzll(eff_c_info);//= select_1(child_info, (eff child)+1)-1
            j-=n_real_lsib*bk_sz;//number of symbols before child j within the node

            //read which child_i has the symbol
            succ_pred_info = bit_pos + symbol*n_children;
            succ_pred_info = stream.read(succ_pred_info, succ_pred_info+n_children-1);
            bit_pos+=new_sigma*n_children;
            //

            //read how many bits we use to encode the pointers to the children
            p_width = stream.read(bit_pos, bit_pos+int_pt_width-1);
            bit_pos+=int_pt_width;
            //

            parent_ptr_area = bit_pos;
            pos = INT_CEIL((bit_pos+(n_children*p_width)), 8)*8;

            //read the pointer to the child for i
            const size_t p = parent_ptr_area+child_i*p_width;
            bit_pos = pos + stream.read(p, p+p_width-1)*8;//add the bit offset. now bit_pos points to child
            //

            //conditions to continue descending
            //child_i==child_j means we descend over the same node
            //(succ_pred_inf>>child_i) & 1 means child_i has "symbol"
            //!stream.read_bit(bit_pos) means child_i is an internal node
            traverse_common_path = child_i==child_j && ((succ_pred_info >> child_i) & 1) && !stream.read_bit(bit_pos++);
            //
#ifdef VLBT_TRACE_RANK
            fprintf(stderr, "[trace]   level bk_sz=%zu: child_i=%llu child_j=%llu  child_info=0x%llx  "
                            "succ_pred_info=0x%llx  n_children=%zu  rank_so_far=%lld  i=%llu j=%llu  continue=%d\n",
                    bk_sz, (unsigned long long)child_i, (unsigned long long)child_j,
                    (unsigned long long)child_info, (unsigned long long)succ_pred_info,
                    n_children, (long long)rank, (unsigned long long)i, (unsigned long long)j,
                    (int)traverse_common_path);
#endif
        } while(traverse_common_path);

        //entering this if means the range of siblings i,i+1,...,j does not contain the symbol
        if((succ_pred_info>>child_i & ((1<<(child_j-child_i+1))-1))==0) {
#ifdef VLBT_TRACE_RANK
            fprintf(stderr, "[trace]   no sibling in [child_i, child_j] holds the symbol -> (%lld, %lld)\n",
                    (long long)rank_i, (long long)rank_j);
#endif
            return std::make_pair(rank_i, rank_j);
        }
        //

        rank_i = rank;
        rank_j = rank;

        //read pred info
        child_info |= 1<<scale_factor; //avoid corner cases for bitwise operations
        const size_t sp_info_i = succ_pred_info & ((1<<(child_i+1))-1);//remove right siblings of child_i
        if(sp_info_i!=0) {//rank for child_i still incomplete (i.e., at least one of the leftmost i children of the parent contains "symbol")

            bool same_child = child_i==child_j;

            //get the rightmost predecessor<=child_i containing "symbol"
            size_t pred = 63-__builtin_clzll(sp_info_i);
            size_t start = stream_type::select64(child_info, pred+1);
            size_t n_real_lsib = __builtin_ctzll(child_info>>(start+1))+1;
            bool following_pred = pred<child_i;
            i = !following_pred*i + following_pred*(bk_sz*n_real_lsib-1);
            child_i = pred;
            //

            //go to the next element
            const size_t p = parent_ptr_area+child_i*p_width;
            bit_pos_i = pos + stream.read(p, p+p_width-1)*8;//add the bit offset. now bit_pos points to child

            if (same_child) {
                //NOTE child_i=child_j is always a leaf when we enter this branch
                //we can't be following the predecessor of i because, if we were, the range i,j would be empty
                assert(!following_pred);
                auto res = scan_leaf<false>(bit_pos_i+1, i, j, symbol, node_sigma, rank_width);
                return std::make_pair(rank_i+res.first, rank_j+res.second);
            }

            const bool is_leaf = stream.read_bit(bit_pos_i++);
            rank_i += subtree_rank(bit_pos_i, i, symbol, node_sigma, rank_width, bk_sz, is_leaf, following_pred);
        }

        const size_t sp_info_j = succ_pred_info & ((1<<(child_j+1))-1);//remove right siblings of j
        if(sp_info_j!=0) {//rank for j still incomplete
            size_t pred = 63-__builtin_clzll(sp_info_j);
            size_t start = stream_type::select64(child_info, pred+1);
            size_t n_real_lsib = __builtin_ctzll(child_info>>(start+1))+1;
            bool following_pred = pred<child_j;
#ifdef VLBT_TRACE_RANK
            const size_t j_before = j;
#endif
            j = !following_pred*j + following_pred*(bk_sz*n_real_lsib-1);
            child_j = pred;

            const size_t p = parent_ptr_area+child_j*p_width;
            bit_pos_j = pos + stream.read(p, p+p_width-1)*8;//add the bit offset. now bit_pos points to child

            const bool is_leaf = stream.read_bit(bit_pos_j++);
#ifdef VLBT_TRACE_RANK
            const int64_t rank_j_before = rank_j;
#endif
            rank_j += subtree_rank(bit_pos_j, j, symbol, node_sigma, rank_width, bk_sz, is_leaf, following_pred);
#ifdef VLBT_TRACE_RANK
            fprintf(stderr, "[trace]   j side: sp_info_j=0x%llx pred=%zu start=%zu n_real_lsib=%zu "
                            "following_pred=%d  j %zu -> %zu  bk_sz=%zu is_leaf=%d  "
                            "rank_j %lld -> %lld (subtree contributed %lld)\n",
                    (unsigned long long)sp_info_j, pred, start, n_real_lsib, (int)following_pred,
                    j_before, j, bk_sz, (int)is_leaf,
                    (long long)rank_j_before, (long long)rank_j, (long long)(rank_j-rank_j_before));
#endif
        }
#ifdef VLBT_TRACE_RANK
        fprintf(stderr, "[trace] result: rank_i=%lld rank_j=%lld\n", (long long)rank_i, (long long)rank_j);
#endif
        return std::make_pair(rank_i, rank_j);
    }

    [[nodiscard]] std::tuple<uint64_t, uint64_t, bool> range_rank(size_t i, size_t j, uint8_t symbol) const {

        symbol = packed_alpha[symbol];

        size_t bit_pos_i, bit_pos_j;
        uint64_t child_i = i/block_size, child_j = j/block_size;
        find_prev_bit_pos(child_i, child_j, bit_pos_i, bit_pos_j);

        bool i_is_head = false;
        int64_t rank_i = 0, rank_j = 0;
        uint8_t rank_width = sym_width(max_freq);

        //i and j are on different trees, meaning there is no advantage, and we proceed as usual
        if(child_i != child_j) {
            const size_t prev_bit_pos_i = bit_pos_i;
            skip_ext_succ_info(bit_pos_i);
            bool is_leaf = stream.read_bit(bit_pos_i++);
            bool has_symbol = stream.read_bit(bit_pos_i+symbol);

            if(!has_symbol){
                rank_i = read_rank_from_succ(prev_bit_pos_i, i, child_i, symbol);
                i_is_head = true;
            }else {
                std::tie(rank_i, i_is_head) = subtree_rank<true>(bit_pos_i, i-child_i*block_size, symbol, sigma, rank_width, block_size, is_leaf, false);
            }

            const size_t prev_bit_pos_j = bit_pos_j;
            skip_ext_succ_info(bit_pos_j);
            is_leaf = stream.read_bit(bit_pos_j++);
            has_symbol = stream.read_bit(bit_pos_j+symbol);

            if (!has_symbol) {
                rank_j = read_rank_from_succ(prev_bit_pos_j, j, child_j, symbol);
            }else {
                rank_j = subtree_rank(bit_pos_j, j-child_j*block_size, symbol, sigma, rank_width, block_size, is_leaf, false);
            }
            return std::make_tuple(rank_i, rank_j, i_is_head);
        }

        size_t bit_pos = bit_pos_i;
        skip_ext_succ_info(bit_pos);

        //positions i and j fall in the same block, but the block does not contain the symbol
        if(!stream.read_bit(bit_pos+1+symbol)){
            return std::make_tuple(0, 0, false);
        }

        //positions i and j fall in the same block and block has the symbol, but the block is a leaf
        size_t offset = child_i*block_size;
        if(stream.read_bit(bit_pos++)){
            return scan_leaf<true>(bit_pos, i-offset, j-offset, symbol, sigma, rank_width);
        }
        //

        //positions i and j fall within the same internal node
        int64_t rank = 0;
        uint8_t node_sigma = sigma;
        size_t bk_sz = block_size;

        bool traverse_common_path = false;
        size_t child_info, succ_pred_info, pos, p_width, parent_ptr_area;
        i-=offset;
        j-=offset;

        do {
            const uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//this works only because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            const size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            node_sigma=new_sigma;

            child_info = stream.read(bit_pos, bit_pos+scale_factor-1);//children info
            bit_pos += scale_factor;

            const size_t n_children = __builtin_popcount(child_info);//number of eff children
            assert(n_children>0);
            bk_sz/=scale_factor;
            rank_width = sym_width(bk_sz*scale_factor);

            child_i = i/bk_sz;
            child_j = j/bk_sz;
            assert(child_i<scale_factor && child_j<scale_factor);

            //read the effective child for i
            size_t eff_c_info = child_info & ((1<<(child_i+1))-1);//clean the bits marking the right siblings
            child_i = __builtin_popcount(eff_c_info)-1;//eff child (zero-based)
            size_t n_real_lsib = 63-__builtin_clzll(eff_c_info);//= select_1(child_info, (eff child)+1)-1
            i-=n_real_lsib*bk_sz;//number of symbols before child i within the node

            //read the effective child for j
            eff_c_info = child_info & ((1<<(child_j+1))-1);//clean the bits marking the right siblings
            child_j = __builtin_popcount(eff_c_info)-1;//eff child (zero-based)
            n_real_lsib = 63-__builtin_clzll(eff_c_info);//= select_1(child_info, (eff child)+1)-1
            j-=n_real_lsib*bk_sz;//number of symbols before child j within the node

            //read which child_i has the symbol
            succ_pred_info = bit_pos + symbol*n_children;
            succ_pred_info = stream.read(succ_pred_info, succ_pred_info+n_children-1);
            bit_pos+=new_sigma*n_children;
            //

            //read how many bits we use to encode the pointers to the children
            p_width = stream.read(bit_pos, bit_pos+int_pt_width-1);
            bit_pos+=int_pt_width;
            //

            parent_ptr_area = bit_pos;
            pos = INT_CEIL((bit_pos+(n_children*p_width)), 8)*8;

            //read the pointer to the child for i
            const size_t p = parent_ptr_area+child_i*p_width;
            bit_pos = pos + stream.read(p, p+p_width-1)*8;//add the bit offset. now bit_pos points to child
            //

            //child_i==child_j means we descend over the same node
            //(succ_pred_inf>>child_i) & 1 means child_i has "symbol"
            //!stream.read_bit(bit_pos) means child_i is an internal node
            traverse_common_path = child_i==child_j && ((succ_pred_info >> child_i) & 1) && !stream.read_bit(bit_pos++);
            //
        } while(traverse_common_path);

        //entering this if means the range of siblings i,i+1,...,j does not contain the symbol
        if((succ_pred_info>>child_i & ((1<<(child_j-child_i+1))-1))==0) {
            return std::make_tuple(rank_i, rank_j, i_is_head);
        }
        //

        rank_i = rank;
        rank_j = rank;

        //read pred info
        child_info |= 1<<scale_factor; //avoid corner cases
        const size_t sp_info_i = succ_pred_info & ((1<<(child_i+1))-1);//remove right siblings of i
        const size_t rank_complete_i = sp_info_i==0;
        i_is_head = rank_complete_i;

        if(!rank_complete_i) {

            bool same_child = child_i==child_j;

            size_t pred = 63-__builtin_clzll(sp_info_i);
            size_t start = stream_type::select64(child_info, pred+1);
            size_t n_real_lsib = __builtin_ctzll(child_info>>(start+1))+1;
            bool following_pred = pred<child_i;
            i = !following_pred*i + following_pred*(bk_sz*n_real_lsib-1);
            child_i = pred;

            const size_t p = parent_ptr_area+child_i*p_width;
            bit_pos_i = pos + stream.read(p, p+p_width-1)*8;//add the bit offset. now bit_pos points to child

            if(same_child){
                //NOTE child_i=child_j is always a leaf when we enter this branch
                //we can't be following the predecessor of i because, if we were, the range i,j would be empty
                assert(!following_pred);
                auto res = scan_leaf<true>(bit_pos_i+1, i, j, symbol, node_sigma, rank_width);
                return std::make_tuple(rank_i+std::get<0>(res), rank_j+std::get<1>(res), std::get<2>(res));
            }

            const bool is_leaf = stream.read_bit(bit_pos_i++);
            auto res_i = subtree_rank<true>(bit_pos_i, i, symbol, node_sigma, rank_width, bk_sz, is_leaf, following_pred);
            rank_i += res_i.first;
            i_is_head = res_i.second;
        }

        const size_t sp_info_j = succ_pred_info & ((1<<(child_j+1))-1);//remove right siblings of j
        const size_t rank_complete_j = sp_info_j==0;
        if(!rank_complete_j) {
            size_t pred = 63-__builtin_clzll(sp_info_j);
            size_t start = stream_type::select64(child_info, pred+1);
            size_t n_real_lsib = __builtin_ctzll(child_info>>(start+1))+1;
            bool following_pred = pred<child_j;
            j = !following_pred*j + following_pred*(bk_sz*n_real_lsib-1);
            child_j = pred;

            const size_t p = parent_ptr_area+child_j*p_width;
            bit_pos_j = pos + stream.read(p, p+p_width-1)*8;//add the bit offset. now bit_pos points to child

            const bool is_leaf = stream.read_bit(bit_pos_j++);
            rank_j += subtree_rank(bit_pos_j, j, symbol, node_sigma, rank_width, bk_sz, is_leaf, following_pred);
        }
        return std::make_tuple(rank_i, rank_j, i_is_head);
    }

    template<bool check_head=false>
    auto subtree_rank(size_t bit_pos, size_t i, uint8_t symbol, size_t node_sigma, size_t rank_width,
                      size_t bk_sz, bool is_leaf, bool following_pred) const {

        stream.prefetch(bit_pos);

        //descend over i
        bool rank_complete=false;
        size_t i_branches[2];
        i_branches[following_pred] = i;
        int64_t rank=0;

        while(!is_leaf && !rank_complete) {

            const uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//this works only because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            const size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            node_sigma=new_sigma;

            bk_sz/=scale_factor;
            size_t child = (following_pred == 0 ? i_branches[0] : i_branches[1])/bk_sz;
            assert(child<scale_factor);

            size_t child_info = stream.read(bit_pos, bit_pos+scale_factor-1);
            bit_pos += scale_factor;

            const size_t n_children = __builtin_popcount(child_info);//number of eff children
            assert(n_children>0);
            size_t pred_info = bit_pos + symbol*n_children;
            stream.prefetch(pred_info);

            child_info &= (1<<(child+1))-1;//clean the bits marking the right siblings
            child = __builtin_popcount(child_info)-1;//eff child (zero-based)
            size_t n_real_lsib = 63-__builtin_clzll(child_info);//= select_1(child_info, (eff child)+1)-1
            i_branches[0]-=n_real_lsib*bk_sz;//number of symbols before child within the node

            //read pred info
            pred_info = stream.read(pred_info, pred_info+n_children-1);
            pred_info &=(1<<(child+1))-1;//remove right siblings
            rank_complete = pred_info==0;

            const size_t pred = 63-__builtin_clzll(pred_info|1);//the |1 is to avoid corner cases with 0
            const size_t start = stream_type::select64(child_info, pred+1);
            child_info |= 1<<scale_factor; //avoid corner cases
            n_real_lsib = __builtin_ctzll(child_info>>(start+1))+1;
            i_branches[1] = bk_sz*n_real_lsib-1;
#ifdef VLBT_TRACE_RANK
            fprintf(stderr, "[trace] subtree_rank  bk_sz=%zu  child=%zu pred=%zu start=%zu n_real_lsib=%zu "
                            "child_info=0x%llx pred_info=0x%llx  following_pred=%d->%d  "
                            "i_branches=[%zu,%zu]  rank_so_far=%lld  rank_complete=%d\n",
                    bk_sz, child, pred, start, n_real_lsib,
                    (unsigned long long)child_info, (unsigned long long)pred_info,
                    (int)following_pred, (int)(following_pred | (pred<child)),
                    i_branches[0], i_branches[1], (long long)rank, (int)rank_complete);
#endif
            following_pred |= pred<child;
            child = pred;
            bit_pos+=new_sigma*n_children;//skip int succ/pred info
            //

            //read how many bits we use to encode the pointers to the children
            const size_t p_width = stream.read(bit_pos, bit_pos+int_pt_width-1);
            bit_pos+=int_pt_width;
            //

            //read the pointer to the child
            size_t p = bit_pos+child*p_width;
            p = stream.read(p, p+p_width-1);
            //

            //skip the pointer to the children and position the bit in the next byte-aligned position
            bit_pos = INT_CEIL((bit_pos+(n_children*p_width)), 8)*8;
            //add the bit offset. now bit_pos points to child
            bit_pos+= p*8;

            rank_width = sym_width(bk_sz*scale_factor);

            //start reading the header of child (there is no ext succ/pred info)
            is_leaf = stream.read_bit(bit_pos++);
        }

        if(!rank_complete) {

            assert(is_leaf);
            i = i_branches[following_pred] + following_pred;//small hack due to the encoding of the leaves
#ifdef VLBT_TRACE_RANK
            //the leaf covers bk_sz symbols, so scanning beyond that reads whatever follows it
            if(i > bk_sz){
                fprintf(stderr, "[trace] *** leaf scan out of range: i=%zu but the leaf covers %zu symbols "
                                "(i_branches[1]=%zu, following_pred=%d)\n",
                        i, bk_sz, i_branches[1], (int)following_pred);
            }
#endif

            const uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//only works because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            const size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            const uint8_t leaf_enc = stream.read(bit_pos, bit_pos+leaf_enc_width-1);
            bit_pos+= leaf_enc_width;

            bool false_break;
            if constexpr (var==RLBWT_WITH_TOEHOLDS){
                if constexpr (check_head){
                    bit_pos+= run_byte_w + run_width;//skip the number of bytes we use to encode the runs
                    false_break = stream.read_bit(bit_pos++);//read if the leftmost run in this head is artificial
                }else{
                    bit_pos+= run_byte_w + run_width + 1;//skip the metadata of the runs
                }
            }

            const uint8_t *leaf_addr = reinterpret_cast<uint8_t *>(stream.stream)+(INT_CEIL(bit_pos, 8));
            int64_t ans;

            //LEAF encoding (bpr=bytes per run):
            //0: 1 bpr, no overflow
            //1: 1 bpr, overflow of 32 elements but not 16
            //2: 1 brp, overflow of 16 and 32 elements

            //3: 2 bpr, no_vbyte, no overflow
            //4: 2 bpr, no_vbyte, overflow 16 elements but not 8
            //5: 2 bpr, no_vbyte, overflow 8 and 16 elements
            //6: 2 bpr, vbyte, no overflow
            //7: 2 bpr, vbyte, overflow of 16 elements but not 8
            //8: 2 bpr, vbyte, overflow of 8 and 16 elements

            //9: 3 bpr, no_vbyte
            //10: 3 bpr, vbyte
            //11: 4 bpr, no_vbyte
            //12: 4 bpr, vbyte

            //13: 5 bpr, vbyte
            //14: 6 bpr, vbyte
            //15: 7 bpr, vbyte

            //scan the runs in the leaf according to the leaf encoding
            switch(leaf_enc) {
                case 0:
                    ans = RANK_8<false, false, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 1:
                    ans = RANK_8<false, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 2:
                    ans = RANK_8<true, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 2 bytes (no vbyte)
                    break;

                case 3://template param: vbyte?, overflow8?, overflow16?
                    ans = RANK_16<false, false, false, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 3 bytes (no vbyte)
                    break;
                case 4:
                    ans = RANK_16<false, false, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (no vbyte)
                    break;
                case 5:
                    ans = RANK_16<false, true, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 5 bytes (no vbyte)
                    break;
                case 6:
                    ans = RANK_16<true, false, false, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 2 bytes (vbyte)
                    break;
                case 7:
                    ans = RANK_16<true, false, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 3 bytes (vbyte)
                    break;
                case 8:
                    ans = RANK_16<true, true, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 9://template param: vbyte?, bpr
                    ans = RANK_32<false,3, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 10:
                    ans = RANK_32<true,3, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 11:
                    ans = RANK_32<false,4, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 12:
                    ans = RANK_32<true,4, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 13://template param: bpr
                    ans = RANK_64<5, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 14:
                    ans = RANK_64<6, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 15:
                    ans = RANK_64<7, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                default:
                    LOG_ERROR("Undefined leaf encoding");
                    exit(1);
            }

            if constexpr (check_head){
                //complete information
                const bool is_same_sym = ans & 1UL;//the query symbol is the same as the symbol of the run where "i" falls
                ans>>=1;
                const bool is_head = ans & 1UL;//the query symbol is the same as the symbol of the run where "i" falls and "i" is the head of that run
                ans>>=1;
                //
                return std::make_pair(rank + ans, (!is_same_sym || (is_head && i>0) || (is_head && !false_break)));
            } else {
                return rank + ans;
            }
        }

        if constexpr (check_head){
            return std::make_pair(rank, true);
        } else {
            return rank;
        }
    }

    [[nodiscard]] int64_t read_rank_from_succ(size_t succ_bit_pos, const size_t i, const size_t child, uint8_t symbol) const {

        //find successor containing sym
        bool succ_found = stream.read_bit(succ_bit_pos+symbol);

        if(!succ_found) {

            uint64_t succ_child = child;
            size_t steps = 0;

            while(!succ_found && steps < 5) {
                succ_bit_pos = find_next(++succ_child);
                skip_ext_succ_info(succ_bit_pos);
                succ_found = stream.read_bit(succ_bit_pos+1+symbol);//does the tree have the symbol?
                steps++;
            }

            if(!succ_found && stream.read_bit(symbol)) {//last opportunity: check if the node is low freq
                succ_bit_pos = find_low_freq_succ(i, symbol);
                skip_ext_succ_info(succ_bit_pos);
                succ_found = true;
            }
        } else {
            decode_ext_succ_info(succ_bit_pos, child, symbol);
            skip_ext_succ_info(succ_bit_pos);
        }

        if(!succ_found) return -1;

        ++succ_bit_pos;//the +1 is to skip the bit indicating if this node is a leaf
        symbol = stream.pop_count(succ_bit_pos, succ_bit_pos+symbol-1);
        succ_bit_pos+=sigma;
        const uint8_t r_width = sym_width(max_freq);
        const size_t r_pos = succ_bit_pos + symbol*r_width;
        const int64_t rank = stream.read(r_pos, r_pos+r_width-1);

        return rank;
    }

    template<bool check_head=false>
    [[nodiscard]] auto rank(size_t i, uint8_t symbol) const {

        //assert(i<=tot_syms);

        symbol = packed_alpha[symbol];
        // NOTE this is a partial rank, because it can sometimes answer -1 for a valid query.
        // However, rank operations in backward_search never return -1, so it is OK for pattern matching

        //initialize the block size
        size_t bk_sz = block_size;

        //get the block where index "i" lies
        uint64_t child = i/bk_sz;

        //bit-position where "child" begins in the stream
        size_t bit_pos = find_prev(child);
        size_t prev_bit_pos=bit_pos;
        skip_ext_succ_info(bit_pos);

        if(bool has_symbol = stream.read_bit(bit_pos+1+symbol); !has_symbol){
            //find successor containing sym
            size_t succ_bit_pos = prev_bit_pos;
            bool succ_found = stream.read_bit(succ_bit_pos+symbol);
            if(!succ_found) {
                uint64_t succ_child = child;
                size_t steps = 0;

                while(!succ_found && steps < 5) {
                    succ_bit_pos = find_next(++succ_child);
                    skip_ext_succ_info(succ_bit_pos);
                    succ_found = stream.read_bit(succ_bit_pos+1+symbol);//does the tree have the symbol?
                    steps++;
                }

                if(!succ_found && stream.read_bit(symbol)) {//last opportunity: check if the node is low freq
                    succ_bit_pos = find_low_freq_succ(i, symbol);
                    skip_ext_succ_info(succ_bit_pos);
                    //assert(stream.read_bit(succ_bit_pos+1+symbol));
                    succ_found = true;
                }
            } else {
                decode_ext_succ_info(succ_bit_pos, child, symbol);
                skip_ext_succ_info(succ_bit_pos);
                //assert(stream.read_bit(succ_bit_pos+1+symbol));
            }

            if(!succ_found){
                if constexpr (check_head){
                    return std::make_pair(static_cast<int64_t>(-1), true);
                }else{
                    return static_cast<int64_t>(-1);
                }
            }

            ++succ_bit_pos;//the +1 is to skip the bit indicating if this node is a leaf
            symbol = stream.pop_count(succ_bit_pos, succ_bit_pos+symbol-1);
            succ_bit_pos+=sigma;
            uint8_t r_width = sym_width(max_freq);
            size_t r_pos = succ_bit_pos + (symbol*r_width);
            int64_t rank = stream.read(r_pos, r_pos+r_width-1);

            if constexpr (check_head){
                return std::make_pair(rank, true);
            }else{
                return rank;
            }
        }
        //

        int64_t rank = 0;
        uint8_t node_sigma = sigma;
        uint8_t rank_width = sym_width(max_freq);

        //read the node header
        bool is_leaf = stream.read_bit(bit_pos++);
        bool rank_complete=false;
        bool following_pred=false;
        size_t i_branches[2]={i, 0};
        i_branches[0]-= child*bk_sz;//relative position of i within the child block

        while(!is_leaf && !rank_complete) {

            uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//this works only because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            node_sigma=new_sigma;

            bk_sz/=scale_factor;
            child = i_branches[following_pred]/bk_sz;
            assert(child<scale_factor);

            size_t child_info = stream.read(bit_pos, bit_pos+scale_factor-1);
            bit_pos += scale_factor;

            size_t n_children = __builtin_popcount(child_info);//number of eff children
            assert(n_children>0);

            child_info &= (1<<(child+1))-1;//clean the bits marking the right siblings
            child = __builtin_popcount(child_info)-1;//eff child (zero-based)
            size_t n_real_lsib = 63-__builtin_clzll(child_info);//= select_1(child_info, (eff child)+1)-1
            i_branches[0]-=n_real_lsib*bk_sz;//number of symbols before child within the node

            //read pred info
            size_t pred_info = bit_pos + (symbol*n_children);
            pred_info = stream.read(pred_info, pred_info+n_children-1);
            pred_info &=(1<<(child+1))-1;//remove right siblings
            rank_complete = pred_info==0;
            size_t pred = 64-__builtin_clzll(pred_info)-!rank_complete;
            size_t start = stream_type::select64(child_info, pred+1);
            child_info |= 1<<scale_factor; //avoid corner cases
            n_real_lsib = __builtin_ctzll(child_info>>(start+1))+1;
            i_branches[1] = (bk_sz*n_real_lsib)-1;
            following_pred |= pred<child;
            child = pred;
            bit_pos+=new_sigma*n_children;//skip int succ/pred info
            //

            //read how many bits we use to encode the pointers to the children
            size_t p_width = stream.read(bit_pos, bit_pos+int_pt_width-1);
            bit_pos+=int_pt_width;
            //

            //read the pointer to the child
            size_t p = bit_pos+(child*p_width);
            p = stream.read(p, p+p_width-1);
            //

            //skip the pointer to the children and position the bit in the next byte-aligned position
            bit_pos = INT_CEIL((bit_pos+(n_children*p_width)), 8)*8;
            //add the bit offset. now bit_pos points to child
            bit_pos+= p*8;

            rank_width = sym_width(bk_sz*scale_factor);
            //start reading the header of child (there is no ext succ/pred info)
            is_leaf = stream.read_bit(bit_pos++);
        }

        if(!rank_complete){

            assert(is_leaf);
            i = i_branches[following_pred] + following_pred;//small hack due to the encoding of the leaves

            uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//only works because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            uint8_t leaf_enc = stream.read(bit_pos, bit_pos+leaf_enc_width-1);
            bit_pos+= leaf_enc_width;

            bool false_break;
            if constexpr (var==RLBWT_WITH_TOEHOLDS){
                if constexpr (check_head){
                    bit_pos+= run_byte_w + run_width;//skip the number of bytes we use to encode the runs
                    false_break = stream.read_bit(bit_pos++);//read if the leftmost run in this head is artificial
                }else{
                    bit_pos+= run_byte_w + run_width + 1;//skip the metadata of the runs
                }
            }

            const uint8_t *leaf_addr = reinterpret_cast<uint8_t *>(stream.stream)+(INT_CEIL(bit_pos, 8));
            int64_t ans;

            //LEAF encoding (bpr=bytes per run):
            //0: 1 bpr, no overflow
            //1: 1 bpr, overflow of 32 elements but not 16
            //2: 1 brp, overflow of 16 and 32 elements

            //3: 2 bpr, no_vbyte, no overflow
            //4: 2 bpr, no_vbyte, overflow 16 elements but not 8
            //5: 2 bpr, no_vbyte, overflow 8 and 16 elements
            //6: 2 bpr, vbyte, no overflow
            //7: 2 bpr, vbyte, overflow of 16 elements but not 8
            //8: 2 bpr, vbyte, overflow of 8 and 16 elements

            //9: 3 bpr, no_vbyte
            //10: 3 bpr, vbyte
            //11: 4 bpr, no_vbyte
            //12: 4 bpr, vbyte

            //13: 5 bpr, vbyte
            //14: 6 bpr, vbyte
            //15: 7 bpr, vbyte

            //scan the runs in the leaf according to the leaf encoding
            switch(leaf_enc) {
                case 0:
                    ans = RANK_8<false, false, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 1:
                    ans = RANK_8<false, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 2:
                    ans = RANK_8<true, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 2 bytes (no vbyte)
                    break;

                case 3://template param: vbyte?, overflow8?, overflow16?
                    ans = RANK_16<false, false, false, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 3 bytes (no vbyte)
                    break;
                case 4:
                    ans = RANK_16<false, false, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (no vbyte)
                    break;
                case 5:
                    ans = RANK_16<false, true, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 5 bytes (no vbyte)
                    break;
                case 6:
                    ans = RANK_16<true, false, false, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 2 bytes (vbyte)
                    break;
                case 7:
                    ans = RANK_16<true, false, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 3 bytes (vbyte)
                    break;
                case 8:
                    ans = RANK_16<true, true, true, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 9://template param: vbyte?, bpr
                    ans = RANK_32<false,3, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 10:
                    ans = RANK_32<true,3, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 11:
                    ans = RANK_32<false,4, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 12:
                    ans = RANK_32<true,4, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 13://template param: bpr
                    ans = RANK_64<5, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 14:
                    ans = RANK_64<6, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 15:
                    ans = RANK_64<7, check_head>(&leaf_addr, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                default:
                    LOG_ERROR("Undefined leaf encoding");
                    exit(1);
            }

            if constexpr (check_head){
                //complete information
                bool is_same_sym = ans & 1UL;//the query symbol is the same as the symbol of the run where "i" falls
                ans>>=1;
                bool is_head = ans & 1UL;//the query symbol is the same as the symbol of the run where "i" falls and "i" is the head of that run
                ans>>=1;
                //
                return std::make_pair(rank + ans, (!is_same_sym || (is_head && i>0) || (is_head && !false_break)));
            }else{
                return rank + ans;
            }
        }

        if constexpr (check_head){
            return std::make_pair(rank, true);
        } else {
            return rank;
        }
    }

    int64_t get_sa_from_leftmost_leaf(succ_info& si, uint8_t pck_sym) const {

        bool is_leaf = stream.read_bit(si.bit_pos++);

        uint64_t bit_pos = si.bit_pos;
        uint8_t symbol = si.symbol;
        uint8_t node_sigma = si.node_sigma;
        uint8_t rank_width = si.r_width;
        uint64_t rank = si.rank;
        size_t bk_sz = si.bk_sz;

        while(!is_leaf){

            uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//this works only because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            node_sigma=new_sigma;

            bk_sz/=scale_factor;

            size_t child_info = stream.read(bit_pos, bit_pos+scale_factor-1);
            bit_pos += scale_factor;

            size_t n_children = __builtin_popcount(child_info);//number of eff children
            assert(n_children>0);

            size_t succ_info = bit_pos + (symbol*n_children);
            succ_info = stream.read(succ_info, succ_info+n_children-1);
            assert(succ_info>0);

            //succ_info>0=true, meaning there is a right sibling containing the symbol
            size_t succ_child = __builtin_ctz(succ_info);

            bit_pos+=new_sigma*n_children;//skip int succ/pred info

            //read how many bits we use to encode the pointers to the children
            size_t p_width = stream.read(bit_pos, bit_pos+int_pt_width-1);
            bit_pos+=int_pt_width;
            //

            //read the pointer to the next sibling of "child" containing the symbol
            size_t p_succ = bit_pos+(succ_child*p_width);
            p_succ = stream.read(p_succ, p_succ+p_width-1);
            //

            //skip the pointer to the children and position the bit in the next byte-aligned position
            bit_pos = INT_CEIL((bit_pos+(n_children*p_width)), 8)*8;
            rank_width = sym_width(bk_sz*scale_factor);

            //add the bit offset. now bit_pos points to child
            bit_pos+= p_succ*8;

            //start reading the header of child (there is no ext succ/pred info)
            is_leaf = stream.read_bit(bit_pos++);
            assert(stream.read_bit(bit_pos+symbol));
        }

        uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
        symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//only works because bit_stream[bit_pos+symbol] is true
        bit_pos+=node_sigma;
        size_t r_pos = bit_pos + symbol*rank_width;
        rank+=stream.read(r_pos, r_pos+rank_width-1);
        bit_pos+=new_sigma*rank_width;

        const uint8_t leaf_enc = stream.read(bit_pos, bit_pos+leaf_enc_width-1);
        bit_pos+= leaf_enc_width;

        const size_t run_bytes= stream.read(bit_pos, bit_pos+run_byte_w-1);
        bit_pos+= run_byte_w;
        const size_t n_runs = stream.read(bit_pos, bit_pos+run_width-1)+1;
        bit_pos+= run_width+1;//+1 to skip the bit indicating if the head of the leftmost run is fake
        size_t byte_pos = INT_CEIL(bit_pos, 8);

        const uint8_t *leaf_addr = reinterpret_cast<uint8_t *>(stream.stream)+byte_pos;

        size_t run_idx;

        //scan the runs in the leaf according to the leaf encoding
        switch(leaf_enc) {
            case 0:
            case 1:
            case 2:
                run_idx = FIRST_RUN_8(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 2 bytes (no vbyte)
                break;

            case 3://template param: vbyte?
            case 4:
            case 5:
                run_idx = FIRST_RUN_16<false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 5 bytes (no vbyte)
                break;
            case 6:
            case 7:
            case 8:
                run_idx = FIRST_RUN_16<true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 4 bytes (vbyte)
                break;

            case 9://template param: vbyte?, bpr
                run_idx = FIRST_RUN_32<false,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 4 bytes (vbyte)
                break;
            case 10:
                run_idx = FIRST_RUN_32<true,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 4 bytes (vbyte)
                break;
            case 11:
                run_idx = FIRST_RUN_32<false,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 4 bytes (vbyte)
                break;
            case 12:
                run_idx = FIRST_RUN_32<true,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 4 bytes (vbyte)
                break;

            case 13://template param: bpr
                run_idx = FIRST_RUN_64<5>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 5 bytes (vbyte)
                break;
            case 14:
                run_idx = FIRST_RUN_64<6>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 5 bytes (vbyte)
                break;
            case 15:
                run_idx = FIRST_RUN_64<7>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, symbol);//runs use 5 bytes (vbyte)
                break;
            default:
                LOG_ERROR("Undefined leaf encoding");
                exit(1);
        }

        bit_pos = (byte_pos + run_bytes)*8;

        bool has_sa_sample = stream.read_bit(bit_pos+int_pt_width+run_idx);
        int64_t sa_samp=0;
        if(has_sa_sample){
            sa_samp = decode_sa_value(bit_pos, run_idx, n_runs);
            return sa_samp;
        }

        size_t n_steps = 0;

        symbol = pck_sym;
        while(!has_sa_sample){
            size_t lf = C[symbol] + rank;
            auto res = inverse_select_with_sa(lf);
            symbol = res.sym;
            rank = res.rank;
            sa_samp =  res.sa_samp;
            has_sa_sample = sa_samp>=0;
            n_steps++;
        }
        assert(n_steps<=subsamp_step);
        return sa_samp+n_steps;
    }

    [[nodiscard]] int64_t decode_sa(size_t bwt_pos) const {
        int64_t sa_samp=std::numeric_limits<int64_t>::min();
        uint64_t n_steps=0;
        while(sa_samp<0 && n_steps<subsamp_step){
            auto res = inverse_select_with_sa(bwt_pos);
            sa_samp =  res.sa_samp;
            bwt_pos = C[res.sym] + res.rank;
            n_steps++;
        }
        return sa_samp+static_cast<int64_t>(n_steps)-1;
    }

    [[nodiscard]] size_t subsampling_value() const {
        return subsamp_step;
    }

    [[nodiscard]] int64_t sa_samp_of_succ_head(size_t i, uint8_t symbol) const {

        static_assert(var==RLBWT_WITH_TOEHOLDS);
        uint8_t pck_sym = packed_alpha[symbol];
        symbol = pck_sym;
        // NOTE this is a partial successor, because it can sometimes answer -1 for a valid query.
        // However, it will never return -1 for a query coming from a pattern that exists in the text

        //initialize the block size
        size_t bk_sz = block_size;

        //succ info
        //s_info[1] contains the successor information
        //s_info[2] is a temporary value discarded when there is no successor
        //we use two variables to avoid branching as we descend over the tree
        succ_info s_info[2];

        //get the block where index "i" lies
        uint64_t child = i/bk_sz, succ_child;
        size_t bit_pos = find_prev(child);//bit-position where "child" begins in the stream

        //we back up this information in case we have to follow a successor tree
        size_t bit_pos_root = bit_pos;
        size_t child_root = child;
        size_t root_i = i;
        //

        skip_ext_succ_info(bit_pos);

        int64_t rank = 0;
        uint8_t node_sigma = sigma;
        uint8_t rank_width = sym_width(max_freq);

        //read the node header
        bool is_leaf = stream.read_bit(bit_pos++);
        bool has_symbol = stream.read_bit(bit_pos+symbol);
        bool succ_found;
        i-= child*bk_sz;//relative position of i within the child block

        while(!is_leaf && has_symbol) {

            uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//this works only because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            node_sigma=new_sigma;

            bk_sz/=scale_factor;
            child = i/bk_sz;
            assert(child<scale_factor);

            size_t child_info = stream.read(bit_pos, bit_pos+scale_factor-1);
            bit_pos += scale_factor;

            size_t n_children = __builtin_popcount(child_info);//number of eff children
            assert(n_children>0);

            child_info &= (1<<(child+1))-1;//clean the bits marking the right siblings
            child = __builtin_popcount(child_info)-1;//eff child (zero-based)
            size_t n_real_lsib = 63-__builtin_clzll(child_info);//= select_1(child_info, (eff child)+1)-1
            i-=n_real_lsib*bk_sz;//number of symbols before child within the node

            //read succ info
            size_t succ_info = bit_pos + (symbol*n_children);
            succ_info = stream.read(succ_info, succ_info+n_children-1);
            //eg: succ_info=1010 for symbol with child=1
            succ_info >>=child+1;//remove "child" and its left siblings
            //succ_info = 10
            succ_found = succ_info>0;
            //succ_info>0=true, meaning there is a right sibling containing the symbol
            succ_child = child+((__builtin_ctz(succ_info)+1)*succ_found);
            //__builtin_ctz(succ_info=10)=1
            //thus, succ_child = child + 1 + 1 = 3
            bit_pos+=new_sigma*n_children;//skip int succ/pred info
            //

            //read how many bits we use to encode the pointers to the children
            size_t p_width = stream.read(bit_pos, bit_pos+int_pt_width-1);
            bit_pos+=int_pt_width;
            //

            //read the pointer to the child
            size_t p = bit_pos+(child*p_width);
            p = stream.read(p, p+p_width-1);
            //

            //read the pointer to the next sibling of "child" containing the symbol
            size_t p_succ = bit_pos+(succ_child*p_width);
            p_succ = stream.read(p_succ, p_succ+p_width-1);
            //

            //skip the pointer to the children and position the bit in the next byte-aligned position
            bit_pos = INT_CEIL((bit_pos+(n_children*p_width)), 8)*8;
            rank_width = sym_width(bk_sz*scale_factor);

            //conditionally calculate successor information
            s_info[succ_found].bit_pos = bit_pos + (p_succ*8);
            s_info[succ_found].symbol = symbol;
            s_info[succ_found].node_sigma = node_sigma;
            s_info[succ_found].rank = rank;
            s_info[succ_found].r_width = rank_width;
            s_info[succ_found].bk_sz = bk_sz;
            //

            //add the bit offset. now bit_pos points to child
            bit_pos+= p*8;

            //start reading the header of child (there is no ext succ/pred info)
            is_leaf = stream.read_bit(bit_pos++);
            has_symbol = stream.read_bit(bit_pos+symbol);
        }

        if(is_leaf && has_symbol) {

            uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//only works because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            uint8_t leaf_enc = stream.read(bit_pos, bit_pos+leaf_enc_width-1);
            bit_pos+= leaf_enc_width;

            size_t run_bytes= stream.read(bit_pos, bit_pos+run_byte_w-1);
            bit_pos+= run_byte_w;
            size_t n_runs = stream.read(bit_pos, bit_pos+run_width-1)+1;
            bit_pos+= run_width+1;//we also skip the bit indicating if the head of the leftmost run is fake
            size_t byte_pos = INT_CEIL(bit_pos, 8);

            const uint8_t *leaf_addr = ((uint8_t *)stream.stream)+(byte_pos);

            std::pair<int64_t, int64_t> ans;//(run_id, rank)

            //LEAF encoding (bpr=bytes per run):
            //0: 1 bpr, no overflow
            //1: 1 bpr, overflow of 32 elements but not 16
            //2: 1 brp, overflow of 16 and 32 elements

            //3: 2 bpr, no_vbyte, no overflow
            //4: 2 bpr, no_vbyte, overflow 16 elements but not 8
            //5: 2 bpr, no_vbyte, overflow 8 and 16 elements
            //6: 2 bpr, vbyte, no overflow
            //7: 2 bpr, vbyte, overflow of 16 elements but not 8
            //8: 2 bpr, vbyte, overflow of 8 and 16 elements

            //9: 3 bpr, no_vbyte
            //10: 3 bpr, vbyte
            //11: 4 bpr, no_vbyte
            //12: 4 bpr, vbyte

            //13: 5 bpr, vbyte
            //14: 6 bpr, vbyte
            //15: 7 bpr, vbyte

            //scan the runs in the leaf according to the leaf encoding
            switch(leaf_enc) {
                case 0:
                    ans = SUCC_8<false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 1:
                    ans = SUCC_8<false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 2:
                    ans = SUCC_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 2 bytes (no vbyte)
                    break;

                case 3://template param: vbyte?, overflow8?, overflow16?
                    ans = SUCC_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 3 bytes (no vbyte)
                    break;
                case 4:
                    ans = SUCC_16<false, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 4 bytes (no vbyte)
                    break;
                case 5:
                    ans = SUCC_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 5 bytes (no vbyte)
                    break;
                case 6:
                    ans = SUCC_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 2 bytes (vbyte)
                    break;
                case 7:
                    ans = SUCC_16<true, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 3 bytes (vbyte)
                    break;
                case 8:
                    ans = SUCC_16<true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 9://template param: vbyte?, bpr
                    ans = SUCC_32<false,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 10:
                    ans = SUCC_32<true,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 11:
                    ans = SUCC_32<false,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 12:
                    ans = SUCC_32<true,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 13://template param: bpr
                    ans = SUCC_64<5>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 14:
                    ans = SUCC_64<6>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 15:
                    ans = SUCC_64<7>(reinterpret_cast<const uint8_t **>(&leaf_addr), n_runs, new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                default:
                    LOG_ERROR("Undefined leaf encoding");
                    exit(1);
            }

            if(ans.first>=0){//it means we found a successor
                bit_pos = (byte_pos + run_bytes)*8;
                bool has_sa_sample = stream.read_bit(bit_pos+int_pt_width+ans.first);

                int64_t sa_samp=0;
                if(has_sa_sample){
                    sa_samp = decode_sa_value(bit_pos, ans.first, n_runs);
                    return sa_samp;
                }

                size_t n_steps = 0;
                rank += ans.second;
                symbol = pck_sym;
                while(!has_sa_sample){
                    size_t lf = C[symbol] + rank;
                    auto res = inverse_select_with_sa(lf);
                    symbol = res.sym;
                    rank = res.rank;
                    sa_samp =  res.sa_samp;
                    has_sa_sample = sa_samp>=0;
                    n_steps++;
                }
                assert(n_steps<=subsamp_step);
                return sa_samp+(int64_t)n_steps;
            }
        }

        //we did not find any successor for "symbol"
        if(s_info[1].bit_pos==0xffffffffffffffff){

            //find the successor tree containing sym
            symbol = pck_sym;
            succ_found = stream.read_bit(bit_pos_root+symbol);
            if(!succ_found) {
                succ_child = child_root;
                size_t steps = 0;

                while(!succ_found && steps < 5) {
                    s_info[1].bit_pos = find_next(++succ_child);
                    skip_ext_succ_info(s_info[1].bit_pos);
                    succ_found = stream.read_bit( s_info[1].bit_pos+1+symbol);//does the tree have the symbol?
                    steps++;
                }

                if(!succ_found && stream.read_bit(symbol)) {//last opportunity: check if the node is low freq
                    s_info[1].bit_pos = find_low_freq_succ(root_i, symbol);
                    skip_ext_succ_info(s_info[1].bit_pos);
                    succ_found = true;
                    //assert(stream.read_bit(succ_bit_pos+1+symbol));
                }
            } else {
                s_info[1].bit_pos = bit_pos_root;
                decode_ext_succ_info(s_info[1].bit_pos, child_root, symbol);
                skip_ext_succ_info(s_info[1].bit_pos);
                //assert(stream.read_bit(succ_bit_pos+1+symbol));
            }

            if(!succ_found) return -1;
            s_info[1].node_sigma = sigma;
            s_info[1].symbol = pck_sym;
            s_info[1].r_width = sym_width(max_freq);
            s_info[1].bk_sz = block_size;
        }

        return get_sa_from_leftmost_leaf(s_info[1], pck_sym);
    }

    void find_path_to_leaf(tree_path_type& path, size_t& i) const {

        //initialize the block size
        size_t bk_sz = block_size;

        //get the block where index "i" lies
        uint64_t child = i/bk_sz;

        //skip the bits with the ext. succ info of low-freq symbols
        path.bit_pos = lfs_bits;

        //get the effective block where "i" lies and its byte offset within the stream
        //[p..p+ext_pt_width-1] is the area where the pointer information of bk lies in the stream
        size_t p = path.bit_pos+(ext_pt_width*child);
        p = stream.read(p, p+ext_pt_width-1);
        size_t offset = (p >> (run_width+1)) * ((p & 1)>0);
        child = child-offset;//eff child in the representation where "i" lies
        p = path.bit_pos + (ext_pt_width*child);
        p = stream.read(p, p+ext_pt_width-1)>>1;

        path.bit_pos =  (header_bytes + p)*8;//bit-position where "child" begins in the stream

        //skip ext succ. information
        const size_t n_samps = stream.pop_count(path.bit_pos, path.bit_pos+sigma-1);
        path.bit_pos+=sigma;
        const uint8_t w = stream.read(path.bit_pos, path.bit_pos+mtd_bits-1);//number bits we use to encode the tree distances for "child"
        path.bit_pos+=mtd_bits;
        path.bit_pos+=n_samps*w;//skip the n_samp tree distances
        path.bit_pos= INT_CEIL(path.bit_pos, 8)*8;//next byte-aligned position
        //
        path.node_sigma[path.lvl] = sigma;

        //read the node header
        ++path.lvl;
        bool is_leaf = stream.read_bit(path.bit_pos++);
        path.node_sigma[path.lvl] = stream.pop_count(path.bit_pos, path.bit_pos+path.node_sigma[path.lvl-1]-1);
        path.sigma_pos[path.lvl]=path.bit_pos;
        path.bit_pos+=path.node_sigma[path.lvl-1];
        path.rank_pos[path.lvl] = path.bit_pos;
        path.rank_width[path.lvl] = sym_width(max_freq);
        path.bit_pos+=path.rank_width[path.lvl]*path.node_sigma[path.lvl];//skip rank information

        i-= child*bk_sz;//relative position of i within the child block

        while(!is_leaf){
            bk_sz/=scale_factor;
            child = i/bk_sz;
            assert(child<scale_factor);

            size_t child_info = stream.read(path.bit_pos, path.bit_pos+scale_factor-1);
            path.bit_pos+=scale_factor;

            size_t n_children = __builtin_popcount(child_info);//number of eff children
            assert(n_children>0);

            child_info &= (1<<(child+1))-1;//clean the bits marking the right siblings
            child = __builtin_popcount(child_info)-1;//eff child (zero-based)
            size_t n_real_lsib = 63-__builtin_clzll(child_info);//= select_1(child_info, (eff child)+1)-1
            i-=n_real_lsib*bk_sz;//number of symbols before child within the node

            path.bit_pos+=path.node_sigma[path.lvl]*n_children;//skip int succ/pred info

            //read how many bits we use to encode the pointers to the children
            size_t p_width = stream.read(path.bit_pos, path.bit_pos+int_pt_width-1);
            path.bit_pos+=int_pt_width;

            p = path.bit_pos+(child*p_width);
            p = stream.read(p, p+p_width-1);

            //skip the pointer to the children and position the bit in the next byte-aligned position
            path.bit_pos = INT_CEIL((path.bit_pos+(n_children*p_width)), 8)*8;
            //add the bit offset. now bit_pos points to child
            path.bit_pos+= p*8;

            //start reading the header of child (there is no ext succ/pred info)
            ++path.lvl;
            is_leaf = stream.read_bit(path.bit_pos++);
            path.node_sigma[path.lvl] = stream.pop_count(path.bit_pos, path.bit_pos+path.node_sigma[path.lvl-1]-1);
            path.sigma_pos[path.lvl]=path.bit_pos;
            path.bit_pos+=path.node_sigma[path.lvl-1];
            path.rank_pos[path.lvl] = path.bit_pos;
            path.rank_width[path.lvl] = sym_width(bk_sz*scale_factor);
            path.bit_pos+=path.rank_width[path.lvl]*path.node_sigma[path.lvl];//skip rank information
        }

        path.leaf_enc = stream.read(path.bit_pos, path.bit_pos+leaf_enc_width-1);
        path.bit_pos+= leaf_enc_width;
        if constexpr (var==RLBWT_WITH_TOEHOLDS){
            path.bit_pos+= run_byte_w+run_width+1;//skip the bits encoding the number of bytes we use to store the runs in their encoding
        }
    }

    [[nodiscard]] std::pair<uint64_t, uint8_t> inverse_select(size_t i) const {

        tree_path_type p;
        find_path_to_leaf(p, i);
        const uint8_t *leaf_addr = ((uint8_t *)stream.stream)+(INT_CEIL(p.bit_pos, 8));

        std::pair<uint64_t, uint8_t> rank_answer;

        //LEAF encoding (bpr=bytes per run):
        //0: 1 bpr, no overflow
        //1: 1 bpr, overflow of 32 elements but not 16
        //2: 1 brp, overflow of 16 and 32 elements

        //3: 2 bpr, no_vbyte, no overflow
        //4: 2 bpr, no_vbyte, overflow 16 elements but not 8
        //5: 2 bpr, no_vbyte, overflow 8 and 16 elements
        //6: 2 bpr, vbyte, no overflow
        //7: 2 bpr, vbyte, overflow of 16 elements but not 8
        //8: 2 bpr, vbyte, overflow of 8 and 16 elements

        //9: 3 bpr, no_vbyte
        //10: 3 bpr, vbyte
        //11: 4 bpr, no_vbyte
        //12: 4 bpr, vbyte

        //13: 5 bpr, vbyte
        //14: 6 bpr, vbyte
        //15: 7 bpr, vbyte

        //scan the runs in the leaf according to the leaf encoding
        switch(p.leaf_enc) {
            case 0:
                rank_answer = INV_SELECT_8<false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 1:
                rank_answer = INV_SELECT_8<false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 2:
                rank_answer = INV_SELECT_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 2 bytes (no vbyte)
                break;

            case 3://template param: vbyte?, overflow8?, overflow16?
                rank_answer = INV_SELECT_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 3 bytes (no vbyte)
                break;
            case 4:
                rank_answer = INV_SELECT_16<false, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (no vbyte)
                break;
            case 5:
                rank_answer = INV_SELECT_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (no vbyte)
                break;
            case 6:
                rank_answer = INV_SELECT_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 2 bytes (vbyte)
                break;
            case 7:
                rank_answer = INV_SELECT_16<true, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 3 bytes (vbyte)
                break;
            case 8:
                rank_answer = INV_SELECT_16<true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;

            case 9://template param: vbyte?, bpr
                rank_answer = INV_SELECT_32<false,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 10:
                rank_answer = INV_SELECT_32<true,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 11:
                rank_answer = INV_SELECT_32<false,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 12:
                rank_answer = INV_SELECT_32<true,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;

            case 13://template param: bpr
                rank_answer = INV_SELECT_64<5>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            case 14:
                rank_answer = INV_SELECT_64<6>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            case 15:
                rank_answer = INV_SELECT_64<7>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            default:
                LOG_ERROR("Undefined leaf encoding");
                exit(1);
        }

        //go back in the path to compute the rank information and the original symbol
        p.bit_pos = p.rank_pos[p.lvl]+(rank_answer.second*p.rank_width[p.lvl]);
        rank_answer.first += stream.read(p.bit_pos, p.bit_pos+p.rank_width[p.lvl]-1);//add rank information
        rank_answer.second = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, rank_answer.second+1);//update the symbol
        --p.lvl;

        while(p.lvl>0){
            p.bit_pos = p.rank_pos[p.lvl]+(rank_answer.second*p.rank_width[p.lvl]);
            rank_answer.first += stream.read(p.bit_pos, p.bit_pos+p.rank_width[p.lvl]-1);
            rank_answer.second = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, rank_answer.second+1);
            --p.lvl;
        }

        return rank_answer;
    }

    [[nodiscard]] uint64_t decode_sa_value(size_t bit_pos, size_t run_id, size_t n_runs) const {
        //NOTE: this function does ont check if run_id is valid
        const size_t sa_width = stream.read(bit_pos, bit_pos + int_pt_width - 1);//read the width
        bit_pos += int_pt_width;
        //NOTE: this operation assumes stream[bit_pos+run_id] is true
        const size_t pos = stream.pop_count(bit_pos, bit_pos + run_id-1);//position of the sample in the encoding
        bit_pos += n_runs;//move to the area where the SA values lie
        bit_pos += sa_width*pos;//move to the area where the SA for run_id lies
        return stream.read(bit_pos, bit_pos + sa_width - 1);
    }

    //this function returns a tuple (symbol, rank, sa_sample)
    //symbol: symbol at index BWT[i]
    //rank: number of occurrences of BWT[i] in the prefix BWT[0..i-1];
    //sa_sample: SA value for BWT[i] *iff* BWT[i] is the head of its run *and* it was subsampled,
    // otherwise sa_sample=-1
    [[nodiscard]] inv_sel_sa_ans inverse_select_with_sa(size_t i) const {

        tree_path_type p;
        find_path_to_leaf(p, i);
        const size_t byte_pos = INT_CEIL(p.bit_pos, 8);

        const uint8_t *leaf_addr = reinterpret_cast<uint8_t *>(stream.stream)+byte_pos;
        inv_sel_sa_ans ans;

        //scan the runs in the leaf according to the leaf encoding
        //NOTE: in this context, INV_SELECT_X returns a tuple (rank (uint64), sym (uint8_t), run_id (int64_t))
        //"sym" is the symbol at position i in the block of runs encoded by the leaf
        //"rank" is number of occurrences of sym in the prefix before the run where "i" falls
        //"run_id" is the id of the run where "i" falls. this value is -1 if "i" is not the head of that run
        switch(p.leaf_enc) {
            case 0:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_8<false, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                            p.node_sigma[p.lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 1:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_8<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                           p.node_sigma[p.lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 2:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_8<true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                          p.node_sigma[p.lvl], i);//runs use 2 bytes (no vbyte)
                break;

            case 3://template param: vbyte?, overflow8?, overflow16?
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_16<false, false, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                                    p.node_sigma[p.lvl], i);//runs use 3 bytes (no vbyte)
                break;
            case 4:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_16<false, false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                                   p.node_sigma[p.lvl], i);//runs use 4 bytes (no vbyte)
                break;
            case 5:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_16<false, true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                                  p.node_sigma[p.lvl], i);//runs use 5 bytes (no vbyte)
                break;
            case 6:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_16<true, false, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                                   p.node_sigma[p.lvl], i);//runs use 2 bytes (vbyte)
                break;
            case 7:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_16<true, false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                                  p.node_sigma[p.lvl], i);//runs use 3 bytes (vbyte)
                break;
            case 8:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_16<true, true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                                 p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;

            case 9://template param: vbyte?, bpr
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_32<false,3, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                        p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 10:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_32<true,3, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                       p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 11:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_32<false,4, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                        p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 12:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_32<true,4, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                       p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;

            case 13://template param: bpr
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_64<5, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                  p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            case 14:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_64<6, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                  p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            case 15:
                std::tie(ans.rank, ans.sym, ans.sa_samp) = INV_SELECT_64<7, true>(reinterpret_cast<const uint8_t **>(&leaf_addr),
                                                                                  p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            default:
                LOG_ERROR("Undefined leaf encoding");
                exit(1);
        }

        //recover the sa sample (if any)
        if(ans.sa_samp>=0) {//"i" is the head of its run run_id=ans.sa_samp
            size_t bit_pos = p.bit_pos - (run_byte_w+run_width+1);
            const size_t run_bytes = stream.read(bit_pos, bit_pos + run_byte_w - 1);//read the number of bytes for the runs
            bit_pos+=run_byte_w;
            const size_t n_runs = stream.read(bit_pos, bit_pos + run_width - 1)+1;
            bit_pos = (byte_pos + run_bytes) * 8;

            //we are now at the start of the area for the SA values. it contains:
            //  the width of the SA values and the number of runs (int_pt_width + run_width)
            //  a bit stream marking each run with a sampled SA value,
            //  the sampled SA values
            if(stream.read_bit(bit_pos + int_pt_width +  ans.sa_samp)) {
                ans.sa_samp = decode_sa_value(bit_pos, ans.sa_samp, n_runs);//ans.sa_samp is the run_id where i lies in the runs
            }else{
                ans.sa_samp = -1;
            }
        }
        //

        //go back in the path to compute the rank information and the original symbol
        p.bit_pos = p.rank_pos[p.lvl]+(ans.sym*p.rank_width[p.lvl]);
        ans.rank += stream.read(p.bit_pos, p.bit_pos+p.rank_width[p.lvl]-1);//add rank information
        ans.sym = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, ans.sym+1);//update the symbol
        --p.lvl;

        while(p.lvl>0){
            p.bit_pos = p.rank_pos[p.lvl]+(ans.sym*p.rank_width[p.lvl]);
            ans.rank += stream.read(p.bit_pos, p.bit_pos+p.rank_width[p.lvl]-1);
            ans.sym = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, ans.sym+1);
            --p.lvl;
        }
        return ans;
    }

    [[nodiscard]] inline uint8_t operator[](size_t i) const {
        assert(i<tot_syms);
        tree_path_type p;
        find_path_to_leaf(p, i);
        const uint8_t *leaf_addr = ((uint8_t *)stream.stream)+(INT_CEIL(p.bit_pos, 8));

        uint8_t symbol;

        //LEAF encoding (bpr=bytes per run):
        //0: 1 bpr, no overflow
        //1: 1 bpr, overflow of 32 elements but not 16
        //2: 1 brp, overflow of 16 and 32 elements

        //3: 2 bpr, no_vbyte, no overflow
        //4: 2 bpr, no_vbyte, overflow 16 elements but not 8
        //5: 2 bpr, no_vbyte, overflow 8 and 16 elements
        //6: 2 bpr, vbyte, no overflow
        //7: 2 bpr, vbyte, overflow of 16 elements but not 8
        //8: 2 bpr, vbyte, overflow of 8 and 16 elements

        //9: 3 bpr, no_vbyte
        //10: 3 bpr, vbyte
        //11: 4 bpr, no_vbyte
        //12: 4 bpr, vbyte

        //13: 5 bpr, vbyte
        //14: 6 bpr, vbyte
        //15: 7 bpr, vbyte

        //scan the runs in the leaf according to the leaf encoding
        switch(p.leaf_enc) {
            case 0:
                symbol = ACCESS_8<false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 1:
                symbol = ACCESS_8<false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 2:
                symbol = ACCESS_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 2 bytes (no vbyte)
                break;

            case 3://template param: vbyte?, overflow8?, overflow16?
                symbol = ACCESS_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 3 bytes (no vbyte)
                break;
            case 4:
                symbol = ACCESS_16<false, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (no vbyte)
                break;
            case 5:
                symbol = ACCESS_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (no vbyte)
                break;
            case 6:
                symbol = ACCESS_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 2 bytes (vbyte)
                break;
            case 7:
                symbol = ACCESS_16<true, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 3 bytes (vbyte)
                break;
            case 8:
                symbol = ACCESS_16<true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;

            case 9://template param: vbyte?, bpr
                symbol = ACCESS_32<false,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 10:
                symbol = ACCESS_32<true,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 11:
                symbol = ACCESS_32<false,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;
            case 12:
                symbol = ACCESS_32<true,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (vbyte)
                break;

            case 13://template param: bpr
                symbol = ACCESS_64<5>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            case 14:
                symbol = ACCESS_64<6>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            case 15:
                symbol = ACCESS_64<7>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (vbyte)
                break;
            default:
                LOG_ERROR("Undefined leaf encoding");
                exit(1);
        }

        //go up in the tree to compute the symbol
        //go back in the path to compute the rank information and the original symbol
        p.bit_pos = p.rank_pos[p.lvl]+(symbol*p.rank_width[p.lvl]);
        symbol = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, symbol+1);//update the symbol
        --p.lvl;

        while(p.lvl>0){
            p.bit_pos = p.rank_pos[p.lvl]+(symbol*p.rank_width[p.lvl]);
            symbol = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, symbol+1);
            --p.lvl;
        }

        return unpacked_alpha[symbol];
    }

    [[nodiscard]] std::pair<uint64_t, uint64_t> count(const std::string &pat) const {
        int64_t l=0, r=static_cast<int64_t>(size())-1;
        size_t j=pat.size();
        while(j-->0 && l<=r){
            const uint8_t cc = packed_alpha[static_cast<uint8_t>(pat[j])];
            //l = C[cc] + rank(l, pat[j]); // count c in bwt[0..l-1]
            //r = C[cc] + rank(r+1, pat[j]) - 1; // count c in bwt[0..r]
            auto [fst, snd] = range_rank_no_sa_head(l, r+1, pat[j]);

            l = C[cc] + fst;
            r = C[cc] + snd-1;
        }
        if(l>r) return {1, 0};
        return {static_cast<uint64_t>(l), static_cast<uint64_t>(r)};
    }

    [[nodiscard]] std::tuple<uint64_t, uint64_t, uint64_t> count_with_head(const std::string &pat) const {
        size_t j=pat.size();
        if (j==0) return {1, 0, 0};
        int64_t l=0, r=static_cast<int64_t>(size())-1;
        std::pair<uint64_t, uint64_t> head[2]={{0,0}, {j-1, l}};

        while(j-->0 && l<=r){

            const uint8_t cc = packed_alpha[static_cast<uint8_t>(pat[j])];

            auto tmp = range_rank(l, r+1, pat[j]);
            head[std::get<2>(tmp)] = {j, l};
            l = C[cc] + std::get<0>(tmp); // count c in bwt[0..l-1]
            r = C[cc] + std::get<1>(tmp) - 1; // count c in bwt[0..r]
        }

        if (l>r) return {1, 0, 0};

        const int64_t sa_samp = sa_samp_of_succ_head(head[1].second, pat[head[1].first]);
        return {static_cast<uint64_t>(l), static_cast<uint64_t>(r), sa_samp-head[1].first-1};
    }

    [[nodiscard]] uint8_t eff2byte(uint8_t eff_sym) const {
        assert(eff_sym<sigma);
        return unpacked_alpha[eff_sym];
    }

    [[nodiscard]] uint64_t size() const {
        return tot_syms;
    }

    [[nodiscard]] size_t alphabet_size() const {
        return sigma;
    }
};

template<size_t b_size>
using vlbt_rlbwt=vlbt_bwt<RLBWT, b_size>;

template<size_t b_size>
using vlbt_rlbwt_th=vlbt_bwt<RLBWT_WITH_TOEHOLDS, b_size>;
#endif //VLBT_BWT_TH_H