//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef VLBT_BWT_TH_H
#define VLBT_BWT_TH_H

#include <cmath>
#include <vector>
#include "bitstream.h"
#include "def_scan.h"

template<size_t b_size, size_t b_runs, size_t s_factor>
struct vlbt_bwt_th {

    static constexpr size_t block_size = b_size;
    static constexpr size_t scale_factor = s_factor;
    static constexpr size_t max_block_runs = b_runs;
    static constexpr uint8_t int_pt_width=6;//number of bits we use to encode the number of bits we use to encode pointers
    static constexpr uint8_t run_width = (sizeof(unsigned long)*8) - __builtin_clzl(b_runs-1);
    static constexpr uint8_t leaf_enc_width=4;
    //number of bits to encode metadata about the sequence of runs:
    // number of bits we use to encode the number of bytes that the runs use in a leaf
    // one extra bit that indicates if the head of the first run is real or artificial
    static constexpr uint8_t runs_mt_bits = ((sizeof(unsigned long)*8) - __builtin_clzl((b_runs*8) + (b_runs/8)))+1;
    typedef bitstream<size_t> stream_type;

    struct tree_path_type{
        uint8_t lvl=0;
        uint64_t bit_pos=0;
        uint64_t sigma_pos[10]={0};
        uint64_t rank_pos[10]={0};
        uint8_t node_sigma[10]={0};
        uint8_t rank_width[10]={0};
        uint8_t leaf_enc=0;
    };

    struct succ_info{
        size_t bit_pos=0xffffffffffffffff;
        uint64_t rank=0;
        uint8_t symbol=0;
        uint8_t r_width=0;
        uint8_t node_sigma=0;

        inline void get_rank(const stream_type& st) {
            bit_pos++;//the +1 is to skip the bit indicating if this node is a leaf
            symbol = st.pop_count(bit_pos, bit_pos+symbol-1);
            bit_pos+=node_sigma;
            size_t r_pos = bit_pos + (symbol*r_width);
            rank += st.read(r_pos, r_pos+r_width-1);
        }
    };

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

    vlbt_bwt_th(): levels(size_t(ceil(log(b_size) / log(s_factor)) - ceil(log(b_runs) / log(s_factor))) + 1){
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

    const std::vector<uint8_t>& get_packed_alpha(){
        return packed_alpha;
    }

    const std::vector<uint8_t>& get_unpacked_alpha(){
        return unpacked_alpha;
    }

    inline size_t find_prev(uint64_t& child) const {
        //get the effective block where i lies and its byte offset within the stream
        //[p..p+ext_pt_width-1] is the area where the pointer information of bk lies in the stream
        size_t p = lfs_bits + (ext_pt_width*child);
        p = stream.read(p, p+ext_pt_width-1);
        size_t offset = (p >> (run_width+1)) * ((p & 1)>0);
        child -= offset;//eff child in the representation where i lies
        size_t child_off = lfs_bits + (ext_pt_width*child);
        p = stream.read(child_off, child_off+ext_pt_width-1)>>1;
        return (header_bytes + p)*8;//bit position where child begins in the stream
    }

    inline size_t find_next(uint64_t& child) const {
        size_t p = lfs_bits + (ext_pt_width*child);
        p = stream.read(p, p+ext_pt_width-1);
        //std::cout<<p<<" "<<int(run_width)<<" "<<" "<<(p>>1)<<" "<<((p >> 1) & run_width)<<std::endl;
        size_t offset = ((p >> 1) & ((1<<run_width)-1)) * ((p & 1)>0);
        child += offset;//eff child in the representation where i lies
        size_t child_off = lfs_bits + (ext_pt_width * child);
        p = stream.read(child_off, child_off + ext_pt_width - 1) >> 1;
        return (header_bytes + p) * 8;
    }

    inline void skip_ext_succ_info(size_t& bit_pos) const {
        //skip ext succ/pred information
        size_t n_samps = stream.pop_count(bit_pos, bit_pos+sigma-1);
        bit_pos+=sigma;
        uint8_t w = stream.read(bit_pos, bit_pos+mtd_bits-1);//number bits we use to encode the tree distances for child
        bit_pos+=mtd_bits;
        bit_pos+=n_samps*w;//skip the n_samp tree distances
        bit_pos= INT_CEIL(bit_pos, 8)*8;//next byte-aligned position (trees are byte-aligned)
        //
    }

    inline void decode_ext_succ_info(size_t& bit_pos, size_t child, uint8_t symbol) const {
        //the position where the offset for the successor is located
        size_t succ_pos = stream.pop_count(bit_pos, bit_pos+symbol)-1;//works only because bit_stream[bit_pos+symbol] is true
        bit_pos+=sigma;
        uint8_t w = stream.read(bit_pos, bit_pos+mtd_bits-1);//number bits we use to encode the tree distances for child
        bit_pos+=mtd_bits + (succ_pos*w);
        size_t succ_child = child + stream.read(bit_pos, bit_pos+w-1);//read the offset of the successor
        size_t p = lfs_bits + (ext_pt_width*succ_child);
        p = stream.read(p, p+ext_pt_width-1)>>1;
        bit_pos = (header_bytes+p)*8;
    }

    inline size_t find_low_freq_succ(size_t i, uint8_t symbol) const {

        symbol = stream.pop_count(0, symbol)-1;//this works because stream[symbol] is true
        size_t c_bits = sigma;//c_bits + (n_symbol+1)*40 contains pointers to the areas where the info lies
        uint8_t w = sym_width(INT_CEIL(tot_syms, block_size)*block_size);

        //read the area of the stream where the info of symbol lies
        size_t ptr = c_bits + symbol*40;
        int64_t first = stream.read(ptr, ptr+39);
        ptr+=40;
        int64_t last = stream.read(ptr, ptr+39)-w;
        size_t n = ((last-first)/w)+1;

        size_t idx;
        while(first<=last){
            int64_t mid = first + int64_t((n/2)*w);
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

    [[nodiscard]] inline int64_t successor(size_t i, uint8_t symbol) const {

        symbol = packed_alpha[symbol];
        // NOTE this is a partial rank, because it can sometimes answer -1 for a valid query.
        // However, rank operations in backwardsearch never return -1, so it is OK for pattern matching

        //initialize the block size
        size_t bk_sz = block_size;

        //succ info
        //s_info[1] contains the successor information
        //s_info[2] is a temporary value that is discarded when there is no successor
        //we use two variables to avoid branching as we descend over the tree
        succ_info s_info[2];

        //get the block where index i lies
        uint64_t child = i/bk_sz, succ_child;
        size_t bit_pos = find_prev(child);//bit position where child begins in the stream

        //find successor containing sym
        bool succ_found = stream.read_bit(bit_pos+symbol);
        size_t succ_bit_pos=0xffffffffffffffff;
        if(!succ_found) {
            succ_child = child;
            size_t steps = 0;
            while(!succ_found && steps < 5) {
                succ_bit_pos = find_next(++succ_child);
                skip_ext_succ_info(succ_bit_pos);
                succ_found = stream.read_bit(succ_bit_pos+1+symbol);//does the tree have the symbol?
                steps++;
            }

            if(!succ_found && stream.read_bit(symbol)) {//check if the node is low-freq
                succ_bit_pos = find_low_freq_succ(i, symbol);
                skip_ext_succ_info(succ_bit_pos);
                succ_found = true;
                //assert(stream.read_bit(succ_bit_pos+1+symbol));
            }
        } else {
            succ_bit_pos = bit_pos;
            decode_ext_succ_info(succ_bit_pos, child, symbol);
            skip_ext_succ_info(succ_bit_pos);
            //assert(stream.read_bit(succ_bit_pos+1+symbol));
        }

        skip_ext_succ_info(bit_pos);

        int64_t rank = 0;
        uint8_t node_sigma = sigma;
        uint8_t rank_width = sym_width(max_freq);

        s_info[succ_found].bit_pos = succ_bit_pos;
        s_info[succ_found].node_sigma = node_sigma;
        s_info[succ_found].symbol = symbol;
        s_info[succ_found].r_width = rank_width;

        //read the node header
        bool is_leaf = stream.read_bit(bit_pos++);
        bool has_symbol = stream.read_bit(bit_pos+symbol);

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
            succ_info >>=child+1;//remove child and its left siblings
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

            //read the pointer to the next sibling of child containing the symbol
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
            //

            //add the bit offset. Now bit_pos points to child
            bit_pos+= p*8;

            //start reading the header of child (there is no ext succ/pred info)
            is_leaf = stream.read_bit(bit_pos++);
            has_symbol = stream.read_bit(bit_pos+symbol);
        }

        if(is_leaf && has_symbol){
            uint8_t new_sigma = stream.pop_count(bit_pos, bit_pos+node_sigma-1);//node_sigma is always >0
            symbol = stream.pop_count(bit_pos, bit_pos+symbol)-1;//only works because bit_stream[bit_pos+symbol] is true
            bit_pos+=node_sigma;
            size_t r_pos = bit_pos + symbol*rank_width;
            rank+=stream.read(r_pos, r_pos+rank_width-1);
            bit_pos+=new_sigma*rank_width;

            uint8_t leaf_enc = stream.read(bit_pos, bit_pos+leaf_enc_width-1);
            bit_pos+= leaf_enc_width;
            bit_pos+= runs_mt_bits;//skip the bits storing the number of bytes used by the runs

            const uint8_t *leaf_addr = ((uint8_t *)stream.stream)+(INT_CEIL(bit_pos, 8));

            //scan the runs in the leaf according to the leaf encoding
            switch(leaf_enc) {
                case 0:
                    rank += RANK_8<false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 1:
                    rank += RANK_8<true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 2:
                    rank += RANK_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 2 bytes (no vbyte)
                    break;

                case 3://template param: vbyte?, overflow8?, overflow16?
                    rank += RANK_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 3 bytes (no vbyte)
                    break;
                case 4:
                    rank += RANK_16<false, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (no vbyte)
                    break;
                case 5:
                    rank += RANK_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (no vbyte)
                    break;
                case 6:
                    rank += RANK_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 2 bytes (vbyte)
                    break;
                case 7:
                    rank += RANK_16<true, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 3 bytes (vbyte)
                    break;
                case 8:
                    rank += RANK_16<true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 9://template param: vbyte?, bpr
                    rank += RANK_32<false,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 10:
                    rank += RANK_32<true,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 11:
                    rank += RANK_32<false,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 12:
                    rank += RANK_32<true,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 13://template param: bpr
                    rank += RANK_64<5>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 14:
                    rank += RANK_64<6>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 15:
                    rank += RANK_64<7>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                default:
                    std::cout<<"Undefined encoding"<<std::endl;
                    exit(1);
            }
            return rank;
        }

        if(s_info[1].bit_pos==0xffffffffffffffff){
            return -1;
        }

        s_info[1].get_rank(stream);
        return s_info[1].rank;
    }

    inline std::pair<int64_t, bool> rank_and_succ(size_t i, uint8_t symbol) {

        symbol = packed_alpha[symbol];
        // NOTE this is a partial rank, because it can sometimes answer -1 for a valid query.
        // However, rank operations in backwardsearch never return -1, so it is OK for pattern matching

        //initialize the block size
        size_t bk_sz = block_size;

        //get the block where index i lies
        uint64_t child = i/bk_sz;

        //bit position where child begins in the stream
        size_t bit_pos = find_prev(child);
        size_t prev_bit_pos=bit_pos;
        skip_ext_succ_info(bit_pos);

        bool has_symbol = stream.read_bit(bit_pos+1+symbol);//does the tree have the symbol?

        if(!has_symbol){
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

                if(!succ_found && stream.read_bit(symbol)) {//last opportunity: check if the node is low-freq
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
                return {-1,true};
            }

            succ_bit_pos++;//the +1 is to skip the bit indicating if this node is a leaf
            symbol = stream.pop_count(succ_bit_pos, succ_bit_pos+symbol-1);
            succ_bit_pos+=sigma;
            uint8_t r_width = sym_width(max_freq);
            size_t r_pos = succ_bit_pos + (symbol*r_width);
            int64_t rank = stream.read(r_pos, r_pos+r_width-1);

            return {rank, true};
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
            rank_width = sym_width(bk_sz*scale_factor);

            //add the bit offset. Now bit_pos points to child
            bit_pos+= p*8;

            //start reading the header of child (there is no ext succ/pred info)
            is_leaf = stream.read_bit(bit_pos++);
            //has_symbol = stream.read_bit(bit_pos+symbol);
        }

        bool succ_is_head = true;

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
            bit_pos+= runs_mt_bits;

            const uint8_t *leaf_addr = ((uint8_t *)stream.stream)+(INT_CEIL(bit_pos, 8));

            //scan the runs in the leaf according to the leaf encoding
            switch(leaf_enc) {
                case 0:
                    rank += RANK_8<false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 1:
                    rank += RANK_8<true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 2:
                    rank += RANK_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 2 bytes (no vbyte)
                    break;

                case 3://template param: vbyte?, overflow8?, overflow16?
                    rank += RANK_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 3 bytes (no vbyte)
                    break;
                case 4:
                    rank += RANK_16<false, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (no vbyte)
                    break;
                case 5:
                    rank += RANK_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (no vbyte)
                    break;
                case 6:
                    rank += RANK_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 2 bytes (vbyte)
                    break;
                case 7:
                    rank += RANK_16<true, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 3 bytes (vbyte)
                    break;
                case 8:
                    rank += RANK_16<true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 9://template param: vbyte?, bpr
                    rank += RANK_32<false,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 10:
                    rank += RANK_32<true,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 11:
                    rank += RANK_32<false,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 12:
                    rank += RANK_32<true,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 13://template param: bpr
                    rank += RANK_64<5>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 14:
                    rank += RANK_64<6>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 15:
                    rank += RANK_64<7>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                default:
                    std::cout<<"Undefined encoding"<<std::endl;
                    exit(1);
            }
        }

        return {rank, succ_is_head};
    }

    inline int64_t rank(size_t i, uint8_t symbol) {

        symbol = packed_alpha[symbol];
        // NOTE this is a partial rank, because it can sometimes answer -1 for a valid query.
        // However, rank operations in backwardsearch never return -1, so it is OK for pattern matching

        //initialize the block size
        size_t bk_sz = block_size;

        //get the block where index i lies
        uint64_t child = i/bk_sz;

        //bit position where child begins in the stream
        size_t bit_pos = find_prev(child);
        size_t prev_bit_pos=bit_pos;
        skip_ext_succ_info(bit_pos);

        bool has_symbol = stream.read_bit(bit_pos+1+symbol);//does the tree have the symbol?

        if(!has_symbol){
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

                if(!succ_found && stream.read_bit(symbol)) {//last opportunity: check if the node is low-freq
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
                return -1;
            }

            succ_bit_pos++;//the +1 is to skip the bit indicating if this node is a leaf
            symbol = stream.pop_count(succ_bit_pos, succ_bit_pos+symbol-1);
            succ_bit_pos+=sigma;
            uint8_t r_width = sym_width(max_freq);
            size_t r_pos = succ_bit_pos + (symbol*r_width);
            int64_t rank = stream.read(r_pos, r_pos+r_width-1);

            return rank;
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
            rank_width = sym_width(bk_sz*scale_factor);

            //add the bit offset. Now bit_pos points to child
            bit_pos+= p*8;

            //start reading the header of child (there is no ext succ/pred info)
            is_leaf = stream.read_bit(bit_pos++);
            //has_symbol = stream.read_bit(bit_pos+symbol);
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
            bit_pos+= runs_mt_bits;

            const uint8_t *leaf_addr = ((uint8_t *)stream.stream)+(INT_CEIL(bit_pos, 8));

            //scan the runs in the leaf according to the leaf encoding
            switch(leaf_enc) {
                case 0:
                    rank += RANK_8<false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 1:
                    rank += RANK_8<true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 1 byte (no vbyte)
                    break;
                case 2:
                    rank += RANK_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 2 bytes (no vbyte)
                    break;

                case 3://template param: vbyte?, overflow8?, overflow16?
                    rank += RANK_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 3 bytes (no vbyte)
                    break;
                case 4:
                    rank += RANK_16<false, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (no vbyte)
                    break;
                case 5:
                    rank += RANK_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (no vbyte)
                    break;
                case 6:
                    rank += RANK_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 2 bytes (vbyte)
                    break;
                case 7:
                    rank += RANK_16<true, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 3 bytes (vbyte)
                    break;
                case 8:
                    rank += RANK_16<true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 9://template param: vbyte?, bpr
                    rank += RANK_32<false,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 10:
                    rank += RANK_32<true,3>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 11:
                    rank += RANK_32<false,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;
                case 12:
                    rank += RANK_32<true,4>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 4 bytes (vbyte)
                    break;

                case 13://template param: bpr
                    rank += RANK_64<5>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 14:
                    rank += RANK_64<6>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                case 15:
                    rank += RANK_64<7>(reinterpret_cast<const uint8_t **>(&leaf_addr), new_sigma, i, symbol);//runs use 5 bytes (vbyte)
                    break;
                default:
                    std::cout<<"Undefined encoding"<<std::endl;
                    exit(1);
            }
        }
        return rank;
    }

    inline void find_path_to_leaf(tree_path_type& path, size_t& i) const {

        //initialize the block size
        size_t bk_sz = block_size;

        //get the block where index i lies
        uint64_t child = i/bk_sz;

        //skip the bits with the ext. succ info of low-freq symbols
        path.bit_pos = lfs_bits;

        //get the effective block where i lies and its byte offset within the stream
        //[p..p+ext_pt_width-1] is the area where the pointer information of bk lies in the stream
        size_t p = path.bit_pos+(ext_pt_width*child);
        p = stream.read(p, p+ext_pt_width-1);
        size_t offset = (p >> (run_width+1)) * ((p & 1)>0);
        //std::cout<<p<<" "<<int(run_mask)<<" "<<int(run_width)<<" byte_pos"<<int(p>>run_width)<<" offset:"<<int(p & run_mask)<<std::endl;
        //child = child-(p & run_mask);//eff child in the representation where i lies
        child = child-offset;//eff child in the representation where i lies
        p = path.bit_pos + (ext_pt_width*child);
        p = stream.read(p, p+ext_pt_width-1)>>1;

        path.bit_pos =  (header_bytes + p)*8;//bit position where child begins in the stream

        //skip ext succ. information
        size_t n_samps = stream.pop_count(path.bit_pos, path.bit_pos+sigma-1);
        path.bit_pos+=sigma;
        uint8_t w = stream.read(path.bit_pos, path.bit_pos+mtd_bits-1);//number bits we use to encode the tree distances for child
        path.bit_pos+=mtd_bits;
        path.bit_pos+=n_samps*w;//skip the n_samp tree distances
        path.bit_pos= INT_CEIL(path.bit_pos, 8)*8;//next byte-aligned position
        //
        path.node_sigma[path.lvl] = sigma;

        //read the node header
        path.lvl++;
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
            //add the bit offset. Now bit_pos points to child
            path.bit_pos+= p*8;

            //start reading the header of child (there is no ext succ/pred info)
            path.lvl++;
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
        path.bit_pos+= runs_mt_bits;//skip the bits encoding the number of bytes we use to store the runs in their encoding
    }

    [[nodiscard]] inline std::pair<uint64_t, uint8_t> inverse_select(size_t i) const {

        tree_path_type p;
        find_path_to_leaf(p, i);
        const uint8_t *leaf_addr = ((uint8_t *)stream.stream)+(INT_CEIL(p.bit_pos, 8));

        std::pair<uint64_t, uint8_t> rank_answer;

        //0 : 1 bpr, no overflow
        //1 : 1 bpr, overflow window=16
        //2 : 1 brp, overflow window=32

        //3 : 2 bpr, no_vbyte, no overflow
        //4 : 2 bpr, no_vbyte, overflow window=8
        //5 : 2 bpr, no_vbyte, overflow window=16
        //6 : 2 bpr, vbyte, no overflow
        //7 : 2 bpr, vbyte, overflow window=8
        //8 : 2 bpr, vbyte, overflow window=16

        //9 : 3 bpr, no_vbyte
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
                rank_answer = INV_SELECT_8<true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 2:
                rank_answer = INV_SELECT_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 2 bytes (no vbyte)
                break;

            case 3://template param: vbyte?, overflow8?, overflow16?
                rank_answer = INV_SELECT_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 3 bytes (no vbyte)
                break;
            case 4:
                rank_answer = INV_SELECT_16<false, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (no vbyte)
                break;
            case 5:
                rank_answer = INV_SELECT_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (no vbyte)
                break;
            case 6:
                rank_answer = INV_SELECT_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 2 bytes (vbyte)
                break;
            case 7:
                rank_answer = INV_SELECT_16<true, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 3 bytes (vbyte)
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
                std::cout<<"Undefined encoding"<<std::endl;
                exit(1);
        }

        //go back in the path to obtain the rank information and the original symbol
        p.bit_pos = p.rank_pos[p.lvl]+(rank_answer.second*p.rank_width[p.lvl]);
        rank_answer.first += stream.read(p.bit_pos, p.bit_pos+p.rank_width[p.lvl]-1);//add rank information
        rank_answer.second = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, rank_answer.second+1);//update the symbol
        p.lvl--;

        while(p.lvl>0){
            p.bit_pos = p.rank_pos[p.lvl]+(rank_answer.second*p.rank_width[p.lvl]);
            rank_answer.first += stream.read(p.bit_pos, p.bit_pos+p.rank_width[p.lvl]-1);
            rank_answer.second = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, rank_answer.second+1);
            p.lvl--;
        }

        return rank_answer;
    }

    [[nodiscard]] inline uint8_t operator[](size_t i) const {
        assert(i<tot_syms);
        tree_path_type p;
        find_path_to_leaf(p, i);
        const uint8_t *leaf_addr = ((uint8_t *)stream.stream)+(INT_CEIL(p.bit_pos, 8));

        uint8_t symbol;
        //LEAF encoding (bpr=bytes per run):

        //0 : 1 bpr, no overflow
        //1 : 1 bpr, overflow window=16
        //2 : 1 brp, overflow window=32

        //3 : 2 bpr, no_vbyte, no overflow
        //4 : 2 bpr, no_vbyte, overflow window=8
        //5 : 2 bpr, no_vbyte, overflow window=16
        //6 : 2 bpr, vbyte, no overflow
        //7 : 2 bpr, vbyte, overflow window=8
        //8 : 2 bpr, vbyte, overflow window=16

        //9 : 3 bpr, no_vbyte
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
                symbol = ACCESS_8<true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 1 byte (no vbyte)
                break;
            case 2:
                symbol = ACCESS_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 2 bytes (no vbyte)
                break;

            case 3://template param: vbyte?, overflow8?, overflow16?
                symbol = ACCESS_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 3 bytes (no vbyte)
                break;
            case 4:
                symbol = ACCESS_16<false, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 4 bytes (no vbyte)
                break;
            case 5:
                symbol = ACCESS_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 5 bytes (no vbyte)
                break;
            case 6:
                symbol = ACCESS_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 2 bytes (vbyte)
                break;
            case 7:
                symbol = ACCESS_16<true, true, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), p.node_sigma[p.lvl], i);//runs use 3 bytes (vbyte)
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
                std::cout<<"Undefined encoding"<<std::endl;
                exit(1);
        }

        //go up in the tree to compute the symbol
        //go back in the path to obtain the rank information and the original symbol
        p.bit_pos = p.rank_pos[p.lvl]+(symbol*p.rank_width[p.lvl]);
        symbol = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, symbol+1);//update the symbol
        p.lvl--;

        while(p.lvl>0){
            p.bit_pos = p.rank_pos[p.lvl]+(symbol*p.rank_width[p.lvl]);
            symbol = stream.select(p.sigma_pos[p.lvl], p.sigma_pos[p.lvl]+p.node_sigma[p.lvl-1]-1, symbol+1);
            p.lvl--;
        }

        return unpacked_alpha[symbol];
    }

    [[nodiscard]] inline uint8_t eff2byte(uint8_t eff_sym) const {
        assert(eff_sym<sigma);
        return unpacked_alpha[eff_sym];
    }

    [[nodiscard]] inline uint64_t size() const {
        return tot_syms;
    }
};
#endif //VLBT_BWT_TH_H