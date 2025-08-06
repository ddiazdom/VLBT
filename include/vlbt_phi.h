//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef VLBT_PHI_H
#define VLBT_PHI_H

#include <cmath>
#include "bit_stream.h"

enum phi_variant {
    NO_VALID_AREA = 0,
    WITH_VALID_AREA = 1
};

template<phi_variant type, size_t b_size, size_t b_runs, size_t s_factor>
struct vlbt_phi {
    typedef bit_stream<size_t> stream_type;
    static constexpr size_t block_size = b_size;
    static constexpr size_t scale_factor = s_factor;
    static constexpr size_t max_block_runs = b_runs;
    static constexpr phi_variant variant = type;

    static constexpr uint8_t int_pt_width = 7;
    //number of bits we use to encode the number of bits we use to encode pointers
    static constexpr uint8_t run_width = (sizeof(unsigned long) * 8) - __builtin_clzl(b_runs - 1);
    //number of bits we use to encode the nearest valid tree in representation
    static constexpr uint8_t leaf_enc_width = 4; //number of bits to encode the leaf encoding
    static constexpr uint8_t run_bytes = (sizeof(unsigned long) * 8) - __builtin_clzl((b_runs * 8) + (b_runs / 8));
    //number of bits we use to encode the number of bytes that the run lengths use in a leaf

    uint64_t tot_syms = 0; //total symbols in the text
    uint16_t ext_pt_width = 0; //number of bits we use store the pointers to the trees
    uint64_t orig_runs = 0; //original number of runs in the BWT
    uint64_t eff_runs = 0; //number of runs in the representation
    uint64_t header_bytes = 0; //bytes used for the header
    uint64_t subsamp_step = 0;//the maximum subsamp step we support is 2^{15}-1= 0x7FFF
    uint8_t levels = 0; //maximum number of levels
    stream_type stream; //stream with the data

    vlbt_phi(): levels(size_t(ceil(log(b_size) / log(s_factor)) - ceil(log(b_runs) / log(s_factor))) + 1) {
        //TODO static asserts in block_size, scale_factor, and b_runs
        // logarithm function to calculate value
        float lg = log(b_size) / log(s_factor);
        assert(lg==floor(lg));
        float lg2 = log(b_runs) / log(s_factor);
        assert(lg2==floor(lg2));
    }

    size_t serialize(std::ostream &ofs) const {
        size_t written_bytes = 0;
        written_bytes += serialize_elm(ofs, tot_syms);
        written_bytes += serialize_elm(ofs, ext_pt_width);
        written_bytes += serialize_elm(ofs, orig_runs);
        written_bytes += serialize_elm(ofs, eff_runs);
        written_bytes += serialize_elm(ofs, header_bytes);
        written_bytes += serialize_elm(ofs, subsamp_step);
        written_bytes += serialize_elm(ofs, levels);
        written_bytes += stream.serialize(ofs);
        return written_bytes;
    }

    void load(std::istream &ifs) {
        load_elm(ifs, tot_syms);
        load_elm(ifs, ext_pt_width);
        load_elm(ifs, orig_runs);
        load_elm(ifs, eff_runs);
        load_elm(ifs, header_bytes);
        load_elm(ifs, subsamp_step);
        load_elm(ifs, levels);
        stream.load(ifs);
    }

    void find_path_to_leaf(uint64_t &bit_pos, size_t& i) const {
        bit_pos = 0;

        //initialize the block size
        size_t bk_sz = block_size;

        //get the block where index "i" lies
        uint64_t child = i / bk_sz;

        //get the effective block where "i" lies and its byte offset within the stream
        //[p..p+ext_pt_width-1] is the area where the pointer information of bk lies in the stream
        size_t p = bit_pos + (ext_pt_width * child);
        p = stream.read(p, p + ext_pt_width - 1);
        const size_t offset = (p >> (run_width + 1)) * ((p & 1) > 0);

        child = child - offset; //eff child in the representation where "i" lies
        p = bit_pos + (ext_pt_width * child);
        p = stream.read(p, p + ext_pt_width - 1) >> 1;

        bit_pos = (header_bytes + p) * 8; //bit-position where "child" begins in the stream

        //read the node header
        bool is_leaf = stream.read_bit(bit_pos++);
        i -= child * bk_sz; //relative position of i within the child block

        while (!is_leaf) {
            bk_sz /= scale_factor;
            child = i / bk_sz;
            assert(child<scale_factor);

            size_t child_info = stream.read(bit_pos, bit_pos + scale_factor - 1);
            bit_pos += scale_factor;

            const size_t n_children = __builtin_popcount(child_info); //number of eff children
            assert(n_children>0);

            child_info &= (1 << (child + 1)) - 1; //clean the bits marking the right siblings
            child = __builtin_popcount(child_info) - 1; //eff child (zero-based)
            const size_t n_real_lsib = 63 - __builtin_clzll(child_info); //= select_1(child_info, (eff child)+1)-1
            i -= n_real_lsib * bk_sz; //number of symbols before child within the node

            //read how many bits we use to encode the pointers to the children
            const size_t p_width = stream.read(bit_pos, bit_pos + int_pt_width - 1);
            bit_pos += int_pt_width;

            p = bit_pos + (child * p_width);
            p = stream.read(p, p + p_width - 1);

            //skip the pointer to the children and position the bit in the next byte-aligned position
            bit_pos = INT_CEIL((bit_pos+(n_children*p_width)), 8) * 8;
            //add the bit offset. now bit_pos points to child
            bit_pos += p * 8;

            is_leaf = stream.read_bit(bit_pos++);
        }
    }


    int64_t operator()(size_t i) const {
        assert(i<tot_syms);
        uint64_t bit_pos;
        const size_t bck_i = i;

        find_path_to_leaf(bit_pos, i);

        const uint8_t leaf_enc = stream.read(bit_pos, bit_pos + leaf_enc_width - 1);
        bit_pos += leaf_enc_width;

        const size_t r_bytes = stream.read(bit_pos, bit_pos + run_bytes - 1);//value of run_width/8
        bit_pos += run_bytes;

        const size_t diff_width = stream.read(bit_pos, bit_pos + int_pt_width - 1); //width we use to store the offsets
        bit_pos += int_pt_width;
        bit_pos = INT_CEIL(bit_pos, 8)*8;//byte aligned

        const uint8_t *leaf_addr = reinterpret_cast<uint8_t *>(stream.stream) + bit_pos/8;

        std::pair<uint64_t, uint64_t> run; // id and offset for the run where "i" falls

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
        switch (leaf_enc) {
            case 0:
                run = GET_PHI_RUN_8<false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 1 byte (no vbyte)
                break;
            case 1:
                run = GET_PHI_RUN_8<false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 1 byte overflow in the sum of 16 symbols
                break;
            case 2:
                //runs use 1 byte overflow in a sum of 32 elements
                run = GET_PHI_RUN_8<true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                break;

            case 3: //template param: vbyte?, overflow8?, overflow16?
                run = GET_PHI_RUN_16<false, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 3 bytes (no vbyte)
                break;
            case 4:
                run = GET_PHI_RUN_16<false, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 4 bytes (no vbyte)
                break;
            case 5:
                run = GET_PHI_RUN_16<false, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 5 bytes (no vbyte)
                break;
            case 6:
                run = GET_PHI_RUN_16<true, false, false>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 2 bytes (vbyte)
                break;
            case 7:
                run = GET_PHI_RUN_16<true, false, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 3 bytes (vbyte)
                break;
            case 8:
                run = GET_PHI_RUN_16<true, true, true>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 4 bytes (vbyte)
                break;

            case 9: //template param: vbyte?, bpr
                run = GET_PHI_RUN_32<false, 3>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 4 bytes (vbyte)
                break;
            case 10:
                run = GET_PHI_RUN_32<true, 3>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 4 bytes (vbyte)
                break;
            case 11:
                run = GET_PHI_RUN_32<false, 4>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 4 bytes (vbyte)
                break;
            case 12:
                run = GET_PHI_RUN_32<true, 4>(reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 4 bytes (vbyte)
                break;

            case 13: //template param: bpr
                run = GET_PHI_RUN_64 <5> (reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 5 bytes (vbyte)
                break;
            case 14:
                run = GET_PHI_RUN_64 <6> (reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 6 bytes (vbyte)
                break;
            case 15:
                run = GET_PHI_RUN_64 <7> (reinterpret_cast<const uint8_t **>(&leaf_addr), i);
                //runs use 7 bytes (vbyte)
                break;
            default:
                std::cout << "Undefined encoding" << std::endl;
                exit(1);
        }

        bit_pos += r_bytes*8;

        if constexpr (variant == NO_VALID_AREA) {
            bit_pos+=run_width;

            bit_pos+=run.first*diff_width;
            const uint64_t val = stream.read(bit_pos, bit_pos + diff_width - 1);
            const uint64_t sa_val = val & 1 ? bck_i-(val>>1UL) : bck_i+(val>>1UL);//the first bit indicates if the different is negative or positive
            return static_cast<int64_t>(sa_val);
        } else {
            const size_t n_runs = stream.read(bit_pos, bit_pos + run_width - 1);
            bit_pos += run_width;
            const size_t val_bit_pos = bit_pos + run.first*diff_width;
            const uint64_t val = stream.read(val_bit_pos, val_bit_pos + diff_width - 1);
            uint64_t sa_val = val & 1 ? bck_i-(val>>1UL) : bck_i+(val>>1UL);//the first bit indicates if the different is negative or positive
            bit_pos += n_runs*diff_width;
            if(!stream.pop_count(bit_pos, bit_pos + run.first)) {//the whole block is a valid area
                return static_cast<int64_t>(sa_val);
            }
            const size_t pos = stream.pop_count(bit_pos, bit_pos + run.first);
            bit_pos += n_runs;
            const size_t w = sym_width(subsamp_step-1);
            bit_pos += pos*w;
            const size_t valid_area = stream.read(bit_pos, bit_pos + w - 1);
            sa_val = valid_area > run.second ?  sa_val : -1;
            return static_cast<int64_t>(sa_val);
        }
    }
};

#endif //VLBT_PHI_H