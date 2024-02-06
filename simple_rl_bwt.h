//
// Created by Diaz, Diego on 19.1.2024.
//

#ifndef SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
#define SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
#include "bitstream.h"
#include "bwt_io.h"
#include <vector>

#define ONE_BYTE_ENC 0
#define TWO_BYTE_ENC 1

struct simple_rl_bwt{

    typedef uint8_t sym_type;

    static const size_t b_size = 4096; //block size
    static const size_t mb_size = 512; //mini block size
    static constexpr uint8_t mb_width = 9; //log2(mb_size)

    static const size_t max_runs_per_block = 128; //threshold to split a block into mini blocks
    static constexpr size_t max_run_len = 2048;//max value we can encode in 11 bits
    static constexpr size_t n_mini_blocks = INT_CEIL(b_size, mb_size);//number of mini blocks of a block
    static constexpr uint8_t mb_header_widths[16] = {0, 12, 24, 36, 48, 60, 72, 84, 96,
                                                     108, 120, 132, 144, 156, 168, 180};//the cumulative bits used by the mini block rank samples

    size_t alphabet=0; //text alphabet
    size_t b_header_bits=0; //number of bits used by the block header
    size_t n_symbols=0; //number of text symbols
    size_t mb_header_bytes=0; //number of bytes used by a mini block header
    size_t tot_runs=0; //number of BWT runs after cutting the text into blocks and mini blocks
    size_t orig_n_runs=0; //original number of BWT runs (before creating the blocks)
    size_t n_sampled_blocks=0; //number of blocks with mini blocks
    uint8_t *data_pointer= nullptr; //pointer to the bwt stream

    std::vector<uint8_t> sym_map; //map symbols to their compacted alphabet
    std::vector<uint8_t> sym_inv_map; //map compacted symbols to their original values
    std::vector<uint16_t> b_header_widths; //cumulative bits for the ranks in the block header
    std::vector<size_t> block_pointers; //position of each block within the BWT stream
    bitstream<size_t> bwt; //BWT stream

    //run encoded in two bytes
    inline size_t insert_run(size_t sym, size_t len, size_t& bwt_pos) {
        assert(len>0);
        uint16_t run = ((len-1) << 4) | sym;
        bwt.write(bwt_pos, bwt_pos + 16 -1, run);
        bwt_pos +=16;
        tot_runs++;
        return 2;
    }

    //run encoded in one byte
    inline size_t insert_mini_run(size_t sym, size_t len, size_t& bwt_pos) {
        //mini run encoded in one byte
        uint8_t run = ((len-1)<< 4) | sym;
        data_pointer[bwt_pos>>3] = run;
        bwt_pos +=8;
        tot_runs++;
        return 1;
    }

    //break a block into mini blocks
    inline void subsample_block(std::vector<std::pair<uint8_t, uint16_t>>& block_runs, size_t& bwt_pos) {

        // store the start of the block just to assert
        // the construction's correctness
        size_t block_start = bwt_pos;

        //metadata within the block
        size_t mb_metadata_offset = bwt_pos;

        //offset of the runs within the block
        bwt_pos +=mb_header_bytes*8;

        //start of the runs of mini blocks within the block
        size_t runs_byte_pos = 0;

        std::vector<std::pair<uint8_t, uint16_t>> sub_block_runs;
        sub_block_runs.reserve(block_runs.size());

        std::vector<size_t> b_sym_freqs(alphabet+1, 0);
        std::vector<uint16_t> mb_rank_widths(alphabet+2, 0);
        for(size_t i=0;i<=alphabet+1;i++) mb_rank_widths[i] = i*12;

        insert_block_header(b_sym_freqs, mb_rank_widths, mb_header_widths[alphabet+1], mb_metadata_offset);

        size_t acc_block=0, broken_run_len, one_byte_run=0;
        for(auto const& run : block_runs) {

            if((acc_block+run.second)>mb_size) {

                //last run of the previous mini block
                broken_run_len = mb_size-acc_block;
                if(broken_run_len>0){
                    one_byte_run += INT_CEIL(broken_run_len, 16);
                    sub_block_runs.emplace_back(run.first, broken_run_len);
                }

                //check which encoding uses fewer bytes : one-byte or two-bytes
                if(one_byte_run<=(sub_block_runs.size()*2)){
                    bwt.write(mb_metadata_offset, mb_metadata_offset, ONE_BYTE_ENC);//mark the mini block as using one-byte encoding
                    for(auto& mini_run : sub_block_runs){
                        while(mini_run.second>16){
                            runs_byte_pos+=insert_mini_run(mini_run.first, 16, bwt_pos);
                            b_sym_freqs[mini_run.first]+= 16;
                            mini_run.second-=16;
                        }
                        if(mini_run.second>0){
                            runs_byte_pos+=insert_mini_run(mini_run.first, mini_run.second, bwt_pos);
                            b_sym_freqs[mini_run.first]+= mini_run.second;
                        }
                    }
                }else{
                    bwt.write(mb_metadata_offset, mb_metadata_offset, TWO_BYTE_ENC);//mark the mini block as using two-byte encoding

                    //a mini block with two-byte encoding has to be two-byte aligned,
                    // so we have to update the pointer in the header if it is not two-byte aligned
                    bool two_byte_unaligned = reinterpret_cast<uintptr_t>(&data_pointer[bwt_pos>>3]) & 1;
                    if(two_byte_unaligned){
                        assert(runs_byte_pos==bwt.read(mb_metadata_offset-12, mb_metadata_offset-1));
                        bwt_pos+=8;
                        runs_byte_pos++;
                    }

                    for(auto const& mini_run : sub_block_runs) {
                        runs_byte_pos+=insert_run(mini_run.first, mini_run.second, bwt_pos);
                        b_sym_freqs[mini_run.first]+= mini_run.second;
                    }
                }
                mb_metadata_offset++;

                //break rules into mini blocks as long as they are
                // bigger than the mini block size
                broken_run_len = (acc_block+run.second)-mb_size;
                while(broken_run_len>mb_size){
                    b_sym_freqs[alphabet] = runs_byte_pos;
                    insert_block_header(b_sym_freqs, mb_rank_widths, mb_header_widths[alphabet+1], mb_metadata_offset);
                    bwt.write(mb_metadata_offset, mb_metadata_offset, TWO_BYTE_ENC);//mark the mini block as using two-byte encoding

                    //a mini block with two-byte encoding has to be two-byte aligned,
                    // so we have to update the pointer it is not
                    bool two_byte_unaligned = reinterpret_cast<uintptr_t>(&data_pointer[bwt_pos>>3]) & 1;
                    if(two_byte_unaligned){
                        bwt_pos+=8;
                        runs_byte_pos++;
                    }

                    runs_byte_pos+=insert_run(run.first, mb_size, bwt_pos);
                    b_sym_freqs[run.first]+=mb_size;
                    broken_run_len-=mb_size;
                    mb_metadata_offset++;
                }

                sub_block_runs.clear();
                one_byte_run=0;
                //insert the header of the first mini run
                b_sym_freqs[alphabet] = runs_byte_pos;
                insert_block_header(b_sym_freqs, mb_rank_widths, mb_header_widths[alphabet+1], mb_metadata_offset);

                //store the first mini run
                sub_block_runs.emplace_back(run.first, broken_run_len);
                one_byte_run+=INT_CEIL(broken_run_len, 16);
                acc_block = broken_run_len;

            } else {
                //the run fits the mini block size
                sub_block_runs.push_back(run);
                one_byte_run+=INT_CEIL(run.second, 16);
                acc_block+=run.second;
            }
        }

        if(one_byte_run<=(sub_block_runs.size()*2)){
            bwt.write(mb_metadata_offset, mb_metadata_offset, ONE_BYTE_ENC);//mark the mini block as using one-byte encoding
            for(auto& mini_run : sub_block_runs){
                while(mini_run.second>16){
                    runs_byte_pos+=insert_mini_run(mini_run.first, 16, bwt_pos);
                    mini_run.second-=16;
                }
                if(mini_run.second>0){
                    runs_byte_pos+=insert_mini_run(mini_run.first, mini_run.second, bwt_pos);
                }
            }
        }else{
            //a mini block with two-byte encoding has to be two-byte aligned,
            // so we have to update the pointer it is not
            bool two_byte_unaligned = reinterpret_cast<uintptr_t>(&data_pointer[bwt_pos>>3]) & 1;
            if(two_byte_unaligned){
                assert(runs_byte_pos==bwt.read(mb_metadata_offset-12, mb_metadata_offset-1));
                bwt_pos+=8;
                runs_byte_pos++;
            }
            bwt.write(mb_metadata_offset, mb_metadata_offset, TWO_BYTE_ENC);//mark the mini block as using one-byte encoding

            for(auto const& mini_run : sub_block_runs) {
                runs_byte_pos+=insert_run(mini_run.first, mini_run.second, bwt_pos);
            }
        }

        //a small assert to check everything is in order
        mb_metadata_offset++;
        assert((mb_metadata_offset-block_start)<=(mb_header_bytes*8));
    }

    inline void insert_block_header(std::vector<size_t>& sym_freq, std::vector<uint16_t>& rank_widths, size_t h_width, size_t& bwt_pos) {
        for(size_t i=0;i<sym_freq.size();i++){
            bwt.write(bwt_pos+rank_widths[i], bwt_pos+rank_widths[i+1]-1, sym_freq[i]);
        }
        bwt_pos+=h_width;
    }

    explicit simple_rl_bwt(std::string& plain_rl_bwt) {

        bwt_buff_reader bwt_buff(plain_rl_bwt);
        size_t n_runs = bwt_buff.size();

        std::vector<size_t> sym_freqs(256, 0);
        size_t sym, len;
        size_t acc_block=0, broken_run_len, eff_runs=0, runs_in_block=0;
        for(size_t i=0;i<n_runs;i++){
            bwt_buff.read_run(i, sym, len);

            if((acc_block+len)>b_size){
                //last run of the previous block
                broken_run_len = b_size-acc_block;
                if(broken_run_len>0){//corner case : the current run is also the start of a block
                    sym_freqs[sym]+=broken_run_len;
                    eff_runs += INT_CEIL(broken_run_len, max_run_len);
                    runs_in_block += INT_CEIL(broken_run_len, max_run_len);
                }

                if(runs_in_block>=max_runs_per_block){
                    n_sampled_blocks++;
                    eff_runs += n_mini_blocks;//each mini blocks adds at most one extra run to the BWT
                }

                //break rules into blocks as long as they are bigger than the block size
                broken_run_len = (acc_block+len)-b_size;
                while(broken_run_len>b_size){
                    eff_runs += INT_CEIL(b_size, max_run_len);
                    sym_freqs[sym]+=b_size;
                    broken_run_len-=b_size;
                }

                //insert the header and first run in the next block
                acc_block = broken_run_len;
                sym_freqs[sym] += broken_run_len;
                eff_runs += INT_CEIL(broken_run_len, max_run_len);
                runs_in_block= INT_CEIL(broken_run_len, max_run_len);
            }else{
                //the run fits the block size
                sym_freqs[sym]+=len;
                acc_block+=len;
                eff_runs+= INT_CEIL(len, max_run_len);
                runs_in_block+= INT_CEIL(len, max_run_len);
            }
        }

        size_t n_syms = 0, u=0;
        sym_map.resize(256);
        sym_inv_map.resize(256);
        b_header_widths.resize(256);
        for(unsigned long sym_freq : sym_freqs){
            n_syms+=sym_freq;
            if(sym_freq>0){
                b_header_widths[alphabet] = sym_width(sym_freq);
                sym_map[u] = alphabet;
                sym_inv_map[alphabet] = u;
                alphabet++;
                n_symbols+=sym_freq;
            }
            u++;
        }

        sym_inv_map.resize(alphabet);
        sym_map.resize(sym_inv_map.back()+1);

        size_t acc=0, tmp;
        b_header_widths.resize(alphabet+1);
        for(size_t i=0;i<alphabet;i++){
            tmp = b_header_widths[i];
            b_header_widths[i] = acc;
            acc +=tmp;
        }
        b_header_widths[alphabet] = acc;
        b_header_bits = INT_CEIL(acc, 8)*8;

        orig_n_runs = bwt_buff.size();

        size_t n_blocks = INT_CEIL(n_syms, b_size);
        block_pointers.resize(n_blocks+1, 0);

        // number of bytes used by the concatenated rank samples of the mini blocks within a block.
        // the alphabet+1 position within the mini block header stores the byte position of the mini block within
        // the block. This position requires 12 bits: in the extreme case that all the runs use one byte, then
        // the highest value is 4096 bytes. The other extreme case is that all the runs use two bytes: then we have
        // at most ceil(4096/9)=456 runs and thus 456*2=912 bytes
        // we add extra bit to indicate if each mini block uses one or two-byte encoding
        mb_header_bytes = INT_CEIL((n_mini_blocks*(mb_header_widths[alphabet+1]+1)), 8);//align metadata to bytes

        //estimate the number of bits in the BWT
        size_t bwt_size_bits = n_blocks*(b_header_bits+8) + n_sampled_blocks*(mb_header_bytes*8) + eff_runs*16;
        bwt.stream_size = INT_CEIL(bwt_size_bits, (sizeof(size_t)*8));
        bwt.stream = (size_t *) malloc(bwt.stream_size*sizeof(size_t));
        data_pointer = (uint8_t *)bwt.stream;

        size_t bwt_pos = 0;
        size_t idx_block=0;

        std::vector<size_t> acc_ranks(alphabet, 0);

        block_pointers[idx_block++] = bwt_pos;
        insert_block_header(acc_ranks, b_header_widths, b_header_bits, bwt_pos);

        acc_block=0;
        std::vector<std::pair<uint8_t, uint16_t>> block_runs;
        block_runs.reserve(b_size);

        for(size_t k=0;k<n_runs;k++){

            bwt_buff.read_run(k, sym, len);
            sym = sym_map[sym];

            if((acc_block+len)>b_size) {

                //last run of the previous block
                broken_run_len = b_size-acc_block;
                if(broken_run_len>0){//corner case : the current run is also the start of a block
                    acc_ranks[sym] += broken_run_len;
                    assert(acc_block+broken_run_len==b_size);
                    while(broken_run_len>max_run_len){
                        block_runs.emplace_back(sym, broken_run_len);
                        broken_run_len-=max_run_len;
                    }
                    if(broken_run_len>0) block_runs.emplace_back(sym, broken_run_len);
                }

                if(block_runs.size()>=max_runs_per_block){
                    //mark the block as subsampled (i.e., it has mini blocks)
                    bwt.write(bwt_pos, bwt_pos+8-1, 1);
                    bwt_pos+=8;
                    subsample_block(block_runs, bwt_pos);
                }else{
                    //mark the block as not subsampled (i.e., it does not have mini blocks)
                    bwt.write(bwt_pos, bwt_pos+8-1, 0);
                    bwt_pos+=8;

                    //make the block two-byte aligned
                    bool two_byte_unaligned = reinterpret_cast<std::uintptr_t>(&data_pointer[bwt_pos>>3]) & 1;
                    bwt_pos +=8*two_byte_unaligned;
                    assert((reinterpret_cast<std::uintptr_t>(&data_pointer[bwt_pos>>3]) & 1)==0);

                    for(auto const& run : block_runs){
                        insert_run(run.first, run.second, bwt_pos);
                    }
                }

                //break rules into blocks as long as they are bigger than the block size
                broken_run_len = (acc_block+len)-b_size;
                while(broken_run_len>b_size){
                    block_pointers[idx_block++] = bwt_pos;
                    insert_block_header(acc_ranks, b_header_widths, b_header_bits, bwt_pos);

                    //mark as not sub sampled
                    bwt.write(bwt_pos, bwt_pos+8-1, 0);
                    bwt_pos+=8;

                    //make the block two-byte aligned
                    bool two_byte_unaligned = reinterpret_cast<std::uintptr_t>(&data_pointer[bwt_pos>>3]) & 1;
                    bwt_pos+=8*two_byte_unaligned;
                    assert((reinterpret_cast<std::uintptr_t>(&data_pointer[bwt_pos>>3]) & 1)==0);

                    acc_ranks[sym]+=b_size;
                    broken_run_len-=b_size;

                    //break the run to fit the max run len
                    len = b_size;
                    while(len>max_run_len){
                        insert_run(sym, max_run_len, bwt_pos);
                        len-=max_run_len;
                    }
                    if(len>0) insert_run(sym, max_run_len, bwt_pos);
                }

                //insert the block header
                block_runs.clear();
                block_pointers[idx_block++] = bwt_pos;
                insert_block_header(acc_ranks, b_header_widths, b_header_bits, bwt_pos);

                //insert the first run
                acc_ranks[sym] += broken_run_len;
                acc_block = broken_run_len;

                //break the run if it exceeds the max run len
                while(broken_run_len>max_run_len){
                    block_runs.emplace_back(sym, max_run_len);
                    broken_run_len-=max_run_len;
                }
                if(broken_run_len>0) block_runs.emplace_back(sym, broken_run_len);

            }else{//the run fits the block size
                acc_ranks[sym]+=len;
                acc_block+=len;

                //break the run if it exceeds the max run len
                while(len>max_run_len){
                    block_runs.emplace_back(sym, max_run_len);
                    len-=max_run_len;
                }
                if(len>0){
                    block_runs.emplace_back(sym, len);
                }
            }
        }

        //insert the last run
        assert(acc_block<=b_size);
        if(block_runs.size()>=max_runs_per_block){
            //mark the block as subsampled (i.e., it has mini blocks)
            bwt.write(bwt_pos, bwt_pos+8-1, 1);
            bwt_pos+=8;
            subsample_block(block_runs, bwt_pos);
        }else{
            //mark the block as not subsampled (i.e., it does not have mini blocks)
            bwt.write(bwt_pos, bwt_pos+8-1, 0);
            bwt_pos+=8;

            //make the block two-byte aligned
            bool two_byte_unaligned = reinterpret_cast<std::uintptr_t>(&data_pointer[bwt_pos>>3]) & 1;
            bwt_pos+=8*two_byte_unaligned;
            assert((reinterpret_cast<std::uintptr_t>(&data_pointer[bwt_pos>>3]) & 1)==0);

            for(auto const& run : block_runs){
                insert_run(run.first, run.second, bwt_pos);
            }
        }

        //appending the full ranks at the end of the encoding
        // to scan the blocks backwards
        block_pointers[idx_block] = bwt_pos;
        insert_block_header(acc_ranks, b_header_widths, b_header_bits, bwt_pos);

        assert(idx_block == n_blocks);
        assert(bwt_pos<=bwt_size_bits);

        //shrink to fit
        bwt.stream_size = INT_CEIL(bwt_pos, (sizeof(size_t)*8));
        bwt.stream = (size_t *) realloc(bwt.stream, bwt.stream_size*sizeof(size_t));
        data_pointer = (uint8_t *)bwt.stream;
    }

    [[nodiscard]] inline std::pair<size_t, sym_type> inverse_select(size_t idx) const  {

        uint16_t b_freq[16] = {0};
        uint8_t symbol;

        size_t rank=0;
        size_t block = idx>>12;
        size_t block_pos = block_pointers[block];

        size_t b_start =  block_pos;
        block_pos += b_header_bits;

        auto *bwt_ptr = data_pointer + (block_pos>>3);
        bool has_mini_blocks = *bwt_ptr;
        block_pos+=8;
        bwt_ptr++;
        size_t tmp_idx = block<<12;

        //the current position indicates if the block has mini blocks or not
        if(has_mini_blocks) {
            //idx of the mini block within the block
            size_t mini_block = (idx-tmp_idx) >> mb_width;

            //read metadata of the mini block
            size_t mb_start = block_pos + mini_block*(mb_header_widths[alphabet+1]+1);

            //read the mini block position within the block and the mini block encoding
            size_t mini_block_pos = bwt.read(mb_start + mb_header_widths[alphabet], mb_start + mb_header_widths[alphabet+1]);
            bool one_byte_encoding = !(mini_block_pos & 4096);
            mini_block_pos &= 4095;//clear the bit indicating the mini block encoding

            tmp_idx += mini_block<< mb_width;

            if((idx-tmp_idx)>(mb_size>>1)){

                size_t next_mb_mt_start;
                if(mini_block<(n_mini_blocks-1)){
                    next_mb_mt_start = block_pos + (mini_block+1)*(mb_header_widths[alphabet+1]+1);
                    mini_block_pos = bwt.read(next_mb_mt_start + mb_header_widths[alphabet],
                                              next_mb_mt_start + mb_header_widths[alphabet+1]-1)-1;
                }else{
                    next_mb_mt_start = block_pointers[block+1];
                    mini_block_pos = next_mb_mt_start-block_pointers[block];
                    mini_block_pos -= b_header_bits + (mb_header_bytes<<3)+8+8;//subtract the header bits
                    mini_block_pos >>=3;//in bytes
                }

                bwt_ptr+= mb_header_bytes + mini_block_pos;
                if(one_byte_encoding){
                    b_scan<true, true>(idx, tmp_idx+mb_size, bwt_ptr, b_freq, symbol);
                }else{
                    b_scan<false, true>(idx, tmp_idx+mb_size, bwt_ptr, b_freq, symbol);
                }

                if(mini_block<(n_mini_blocks-1)){
                    rank = bwt.read(b_start + b_header_widths[symbol], b_start + b_header_widths[symbol+1] -1);
                    rank += bwt.read(next_mb_mt_start + mb_header_widths[symbol], next_mb_mt_start + mb_header_widths[symbol+1]-1);
                }else{
                    rank = bwt.read(next_mb_mt_start + b_header_widths[symbol], next_mb_mt_start + b_header_widths[symbol+1] -1);
                }
                rank-=b_freq[symbol];
            }else{
                bwt_ptr+= mb_header_bytes + mini_block_pos;
                if(one_byte_encoding){
                    f_scan<true, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);
                }else{
                    bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
                    f_scan<false, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);
                }
                rank = bwt.read(b_start + b_header_widths[symbol], b_start + b_header_widths[symbol+1] -1);
                rank+= bwt.read(mb_start + mb_header_widths[symbol], mb_start + mb_header_widths[symbol+1] -1);
                rank+=b_freq[symbol];
            }
        } else {
            if((idx-tmp_idx)>(b_size>>1)){
                block_pos = block_pointers[block+1]-block_pointers[block];
                block_pos -= b_header_bits + 16;//subtract the header bits 8+8=16 consider the mini block flag and moving one position back from the end
                block_pos >>=3;//in bytes
                bwt_ptr+=block_pos;

                b_scan<false, true>(idx, tmp_idx+b_size, bwt_ptr, b_freq, symbol);
                b_start = block_pointers[block+1];
                rank = bwt.read(b_start + b_header_widths[symbol], b_start + b_header_widths[symbol+1] -1);
                rank-=b_freq[symbol];
            }else{
                bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
                f_scan<false, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);
                rank = bwt.read(b_start + b_header_widths[symbol], b_start + b_header_widths[symbol+1] -1);
                rank+=b_freq[symbol];
            }
        }
        return {rank, sym_inv_map[symbol]};
    }

    template<bool one_byte_encoding, bool perform_count>
    static inline void f_scan(size_t idx, size_t tmp_idx, uint8_t * bwt_ptr, uint16_t* b_freq, uint8_t& symbol){
        uint16_t data;
        if constexpr (one_byte_encoding){
            while(tmp_idx<=idx){
                //get the run symbol
                data = *bwt_ptr;
                symbol = data & 15;

                //get the run len
                data>>=4;
                data++;

                if constexpr (perform_count){
                    b_freq[symbol]+=data;
                }

                //move to the next position
                tmp_idx+=data;
                bwt_ptr++;
            }

            if constexpr (perform_count){
                b_freq[symbol]-=tmp_idx-idx;
            }
        }else{
            auto *run = (uint16_t *) bwt_ptr;
            while(tmp_idx<=idx) {
                //I assume the compiler will unroll this loop
                for(size_t i=0;i<4;i++){
                    symbol = run[i] & 15;
                    data = (run[i]>>4) + 1;
                    if constexpr (perform_count){
                        b_freq[symbol]+=data;
                    }
                    tmp_idx+=data;
                }
                run+=4;
            }

            run-=4;
            size_t pos = 3;
            while(tmp_idx>idx){
                symbol = run[pos] & 15;
                data = (run[pos]>>4) + 1;

                if constexpr (perform_count){
                    b_freq[symbol]-=data;
                }

                tmp_idx-=data;
                pos--;
            }

            if constexpr (perform_count){
                b_freq[symbol]+=idx-tmp_idx;
            }
        }
    }

    //TODO pad the BWT stream with zeroes
    template<bool one_byte_encoding, bool perform_count>
    static inline void b_scan(size_t idx, size_t tmp_idx, uint8_t * bwt_ptr, uint16_t* b_freq,
                              uint8_t& symbol){
        uint16_t data;
        if constexpr (one_byte_encoding) {
            while(tmp_idx>idx){
                //get the run symbol
                data = *bwt_ptr;
                symbol = data & 15;

                //get the run len
                data>>=4;
                data++;

                if constexpr (perform_count){
                    b_freq[symbol]+=data;
                }

                //move to the next position
                tmp_idx-=data;
                bwt_ptr--;
            }
            if constexpr (perform_count){
                b_freq[symbol]-=idx-tmp_idx;
            }
        }else{
            bwt_ptr--;
            auto *run = (uint16_t *) bwt_ptr;

            while(tmp_idx>idx) {
                //I assume the compiler will unroll this loop
                for(size_t i=0;i<4;i++){
                    symbol = run[0] & 15;
                    data = (run[0]>>4) + 1;
                    if constexpr (perform_count){
                        b_freq[symbol]+=data;
                    }
                    tmp_idx-=data;
                    run--;
                }
            }

            run++;
            size_t pos = 0;
            while(tmp_idx<=idx){
                symbol = run[pos] & 15;
                data = (run[pos] >> 4) + 1;
                if constexpr (perform_count){
                    b_freq[symbol]-=data;
                }
                tmp_idx+=data;
                pos++;
            }
            if constexpr (perform_count){
                b_freq[symbol]+=tmp_idx-idx;
            }
        }
    }

    [[nodiscard]] inline size_t rank(size_t idx, sym_type symbol) const {

        uint16_t b_freq[16] = {0};
        uint8_t q_symbol = sym_map[symbol];

        size_t rank;
        size_t block = idx>>12;
        size_t block_pos = block_pointers[block];

        size_t b_start =  block_pos;
        block_pos += b_header_bits;

        auto *bwt_ptr = data_pointer + (block_pos>>3);
        bool has_mini_blocks = *bwt_ptr;
        block_pos+=8;
        bwt_ptr++;
        size_t tmp_idx = block<<12;

        //the current position indicates if the block has mini blocks or not
        if(has_mini_blocks) {
            //idx of the mini block within the block
            size_t mini_block = (idx-tmp_idx) >> mb_width;

            //read metadata of the mini block
            size_t mb_start = block_pos + mini_block*(mb_header_widths[alphabet+1]+1);

            //read the mini block position within the block and the mini block encoding
            size_t mini_block_pos = bwt.read(mb_start + mb_header_widths[alphabet], mb_start + mb_header_widths[alphabet+1]);
            bool one_byte_encoding = !(mini_block_pos & 4096);
            mini_block_pos &= 4095;//clear the bit indicating the mini block encoding

            tmp_idx += mini_block<< mb_width;

            if((idx-tmp_idx)>(mb_size>>1)){

                size_t next_mb_mt_start;
                if(mini_block<(n_mini_blocks-1)){
                    next_mb_mt_start = block_pos + (mini_block+1)*(mb_header_widths[alphabet+1]+1);
                    mini_block_pos = bwt.read(next_mb_mt_start + mb_header_widths[alphabet],
                                              next_mb_mt_start + mb_header_widths[alphabet+1]-1)-1;
                }else{
                    next_mb_mt_start = block_pointers[block+1];
                    mini_block_pos = next_mb_mt_start-block_pointers[block];
                    mini_block_pos -= b_header_bits + (mb_header_bytes<<3)+8+8;//subtract the header bits
                    mini_block_pos >>=3;//in bytes
                }

                bwt_ptr+= mb_header_bytes + mini_block_pos;
                if(one_byte_encoding){
                    b_scan<true, true>(idx, tmp_idx+mb_size, bwt_ptr, b_freq, symbol);
                }else{
                    b_scan<false, true>(idx, tmp_idx+mb_size, bwt_ptr, b_freq, symbol);
                }

                if(mini_block<(n_mini_blocks-1)){
                    rank = bwt.read(b_start + b_header_widths[q_symbol], b_start + b_header_widths[q_symbol+1] -1);
                    rank += bwt.read(next_mb_mt_start + mb_header_widths[q_symbol], next_mb_mt_start + mb_header_widths[q_symbol+1]-1);
                }else{
                    rank = bwt.read(next_mb_mt_start + b_header_widths[q_symbol], next_mb_mt_start + b_header_widths[q_symbol+1] -1);
                }
                rank-=b_freq[q_symbol];
            }else{
                bwt_ptr+= mb_header_bytes + mini_block_pos;
                if(one_byte_encoding){
                    f_scan<true, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);
                }else{
                    bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
                    f_scan<false, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);
                }
                rank = bwt.read(b_start + b_header_widths[q_symbol], b_start + b_header_widths[q_symbol+1] -1);
                rank+= bwt.read(mb_start + mb_header_widths[q_symbol], mb_start + mb_header_widths[q_symbol+1] -1);
                rank+=b_freq[q_symbol];
            }
        } else {
            if((idx-tmp_idx)>(b_size>>1)){
                block_pos = block_pointers[block+1]-block_pointers[block];
                block_pos -= b_header_bits + 16;//subtract the header bits 8+8=16 consider the mini block flag and moving one position back from the end
                block_pos >>=3;//in bytes
                bwt_ptr+=block_pos;

                b_scan<false, true>(idx, tmp_idx+b_size, bwt_ptr, b_freq, symbol);
                b_start = block_pointers[block+1];
                rank = bwt.read(b_start + b_header_widths[q_symbol], b_start + b_header_widths[q_symbol+1] -1);
                rank-=b_freq[q_symbol];
            }else{
                bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
                f_scan<false, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);
                rank = bwt.read(b_start + b_header_widths[q_symbol], b_start + b_header_widths[q_symbol+1] -1);
                rank+=b_freq[q_symbol];
            }
        }
        return rank;
    }

    inline void count_in_block(size_t tmp_idx, size_t idx, size_t block_pos,
                               size_t block, uint8_t* bwt_ptr,
                               std::vector<size_t>& ranks) const {

        sym_type symbol;
        uint16_t b_freq[16]={0};

        if((idx-tmp_idx)>(b_size>>1)) {
            block_pos = block_pointers[block+1]-block_pointers[block];
            block_pos -= b_header_bits + 16;//subtract the header bits 8+8=16 consider the mini block flag and moving one position back from the end
            block_pos >>=3;//in bytes
            bwt_ptr+= block_pos;

            b_scan<false, true>(idx, tmp_idx+b_size, bwt_ptr, b_freq, symbol);
            size_t b_start = block_pointers[block+1];
            for(size_t u=0;u<alphabet;u++){
                ranks[u] = bwt.read(b_start + b_header_widths[u], b_start + b_header_widths[u+1] -1);
                ranks[u]-= b_freq[u];
            }
        } else{
            bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
            f_scan<false, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);

            size_t b_start = block_pointers[block];
            for(size_t u=0;u<alphabet;u++){
                ranks[u] = bwt.read(b_start + b_header_widths[u], b_start + b_header_widths[u+1] -1);
                ranks[u]+= b_freq[u];
            }
        }
    }

    inline void count_in_mini_block(size_t tmp_idx, size_t idx, size_t block_pos,
                                    size_t block, uint8_t* bwt_ptr,
                                    std::vector<size_t>& ranks) const {

        sym_type symbol;
        uint16_t b_freq[16]={0};

        size_t mini_block = (idx-tmp_idx) >> mb_width;

        //read metadata of the mini block
        size_t mb_start = block_pos + mini_block*(mb_header_widths[alphabet+1]+1);

        //read the mini block position within the block and the mini block encoding
        size_t mini_block_pos = bwt.read(mb_start + mb_header_widths[alphabet], mb_start + mb_header_widths[alphabet+1]);
        bool one_byte_encoding = !(mini_block_pos & 4096);
        mini_block_pos &= 4095;//clear the bit indicating the mini block encoding
        tmp_idx += mini_block<< mb_width;

        if((idx-tmp_idx)>(mb_size>>1)) {

            size_t next_mb_mt_start;
            if(mini_block<(n_mini_blocks-1)){
                next_mb_mt_start = block_pos + (mini_block+1)*(mb_header_widths[alphabet+1]+1);
                mini_block_pos = bwt.read(next_mb_mt_start + mb_header_widths[alphabet],
                                          next_mb_mt_start + mb_header_widths[alphabet+1]-1)-1;
            }else{
                next_mb_mt_start = block_pointers[block+1];
                mini_block_pos = next_mb_mt_start-block_pointers[block];
                mini_block_pos -= b_header_bits + (mb_header_bytes<<3)+8+8;//subtract the header bits
                mini_block_pos >>=3;//in bytes
            }

            bwt_ptr+= mb_header_bytes + mini_block_pos;
            if(one_byte_encoding){
                b_scan<true, true>(idx, tmp_idx+mb_size, bwt_ptr, b_freq, symbol);
            }else{
                b_scan<false, true>(idx, tmp_idx+mb_size, bwt_ptr, b_freq, symbol);
            }

            if(mini_block<(n_mini_blocks-1)){
                //size_t b_start = block_pointers[block];
                for(size_t u=0;u<alphabet;u++){
                    ranks[u]  = bwt.read(block_pointers[block] + b_header_widths[u], block_pointers[block] + b_header_widths[u+1] -1);
                    ranks[u] += bwt.read(next_mb_mt_start + mb_header_widths[u], next_mb_mt_start + mb_header_widths[u+1]-1);
                    ranks[u] -= b_freq[u];
                }
            }else{
                for(size_t u=0;u<alphabet;u++){
                    ranks[u] = bwt.read(next_mb_mt_start + b_header_widths[u], next_mb_mt_start + b_header_widths[u+1] -1);
                    ranks[u] -= b_freq[u];
                }
            }
        }else{
            bwt_ptr+= mb_header_bytes + mini_block_pos;
            if(one_byte_encoding){
                f_scan<true, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);
            }else{
                bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
                f_scan<false, true>(idx, tmp_idx, bwt_ptr, b_freq, symbol);
            }

            //size_t b_start = block_pointers[block];
            for(size_t u=0;u<alphabet;u++){
                ranks[u] = bwt.read(block_pointers[block] + b_header_widths[u], block_pointers[block] + b_header_widths[u+1] -1);
                ranks[u]+= bwt.read(mb_start + mb_header_widths[u], mb_start + mb_header_widths[u+1] -1);
                ranks[u]+= b_freq[u];
            }
        }
    }

    inline void interval_symbols(size_t i, size_t j, size_t& k,
                                 std::vector<sym_type>& cs,
                                 std::vector<size_t>& rank_c_i,
                                 std::vector<size_t>& rank_c_j)  {

        size_t i_block = i>>12;
        size_t i_block_pos = block_pointers[i_block]+b_header_bits;
        uint8_t * i_bwt_ptr = data_pointer + (i_block_pos>>3);
        bool i_has_mini_block = *i_bwt_ptr;
        i_bwt_ptr++;
        i_block_pos+=8;
        size_t tmp_i = i_block<<12;

        size_t j_block = j>>12;

        if(i_block==j_block) {//i and j-1 are in the same block b
            if(i_has_mini_block){ //b has mini blocks
                count_in_mini_block(tmp_i, i, i_block_pos, i_block, i_bwt_ptr, rank_c_i);
                count_in_mini_block(tmp_i, j, i_block_pos, i_block, i_bwt_ptr, rank_c_j);
            } else { //b does not have mini blocks
                count_in_block(tmp_i, i, i_block_pos, i_block, i_bwt_ptr, rank_c_i);
                count_in_block(tmp_i, j, i_block_pos, i_block, i_bwt_ptr, rank_c_j);
            }
        } else {// i and j-1 are in different blocks

            if(i_has_mini_block) { //i's block has mini blocks
                //idx of the mini block within the block
                count_in_mini_block(tmp_i, i, i_block_pos, i_block, i_bwt_ptr, rank_c_i);
            }else{ //i's block does not have mini blocks
                count_in_block(tmp_i, i, i_block_pos, i_block, i_bwt_ptr, rank_c_i);
            }

            //process j-1
            size_t j_block_pos = block_pointers[j_block]+b_header_bits;
            uint8_t * j_bwt_ptr = data_pointer + (j_block_pos>>3);
            bool j_has_mini_blocks = *j_bwt_ptr;
            j_bwt_ptr+=1;
            j_block_pos+=8;
            size_t tmp_j = j_block<<12;

            if(j_has_mini_blocks){ //j's block has mini blocks
                count_in_mini_block(tmp_j, j, j_block_pos, j_block, j_bwt_ptr, rank_c_j);
            }else{ //j's block does not have mini blocks
                count_in_block(tmp_j, j, j_block_pos, j_block, j_bwt_ptr, rank_c_j);
            }
        }

        /*size_t i_block = i>>12;
        j--;
        size_t j_block = j>>12;

        size_t i_block_pos = block_pointers[i_block];
        size_t j_block_pos;

        //sample the ranks up to the respective blocks
        if(j_block==i_block){
            //copy the ranks if i and j-1 are withing the same block
            j_block_pos = i_block_pos;
            for(size_t u=0;u<alphabet;u++){
                rank_c_i[u] = bwt.read(i_block_pos+b_header_widths[u], i_block_pos+b_header_widths[u+1]-1);
            }
            memcpy(rank_c_j.data(), rank_c_i.data(), alphabet*sizeof(size_t));
        }else{
            j_block_pos = block_pointers[j_block];
            for(size_t u=0;u<alphabet;u++){
                rank_c_i[u] = bwt.read(i_block_pos+b_header_widths[u], i_block_pos+b_header_widths[u+1]-1);
                rank_c_j[u] = bwt.read(j_block_pos+b_header_widths[u], j_block_pos+b_header_widths[u+1]-1);
            }
        }
        i_block_pos += b_header_bits;
        j_block_pos += b_header_bits;

        //scan the block for i
        uint8_t * i_bwt_ptr = data_pointer + (i_block_pos>>3);
        i_block_pos+=8;
        size_t tmp_i = i_block<<12;
        size_t i_mini_block=0;

        //the current position indicates if the block was sub sampled or not
        if(*i_bwt_ptr){
            //idx of the mini block within the block
            i_mini_block = (i-tmp_i) >> mb_width;

            //metadata of the mini block
            size_t i_mb_start = i_block_pos + i_mini_block*mb_header_widths[alphabet+1];

            //read the ranks
            size_t tmp_val, off_set=0, pos;
            for(size_t u=0;u<alphabet;u++){
                pos = i_mb_start+off_set;
                memcpy(&tmp_val, &data_pointer[pos>>3], 2);
                rank_c_i[u] += (tmp_val>> (pos & 7)) & 4095;
                //rank_c_i[u] += bwt.read(i_mb_start + off_set, i_mb_start + off_set+11);
                off_set+=12;
            }

            //read the mini block position within the block
            i_bwt_ptr += mb_header_bytes;//+ mini_block_pos;
            i_bwt_ptr += bwt.read(i_mb_start + mb_header_widths[alphabet], i_mb_start + mb_header_widths[alphabet+1] -1);
            tmp_i += i_mini_block<< mb_width;
        }

        //scan the block for j
        uint8_t *j_bwt_ptr = data_pointer + (j_block_pos>>3);
        j_block_pos+=8;
        size_t tmp_j = j_block<<12;
        bool same_mini_block=false;

        //the current position indicates if the block was sub sampled
        if(*j_bwt_ptr){
            //idx of the mini block within the block
            size_t j_mini_block = (j-tmp_j) >> mb_width;

            //i and j-1 are within the same mini block, so copy the information
            same_mini_block = i_block==j_block && i_mini_block==j_mini_block;
            if(same_mini_block){
                tmp_j = tmp_i;
                j_bwt_ptr = i_bwt_ptr;
            }else{//not in the same mini block

                //the bit position within the BWT for the mini block metadata
                size_t j_mb_start = j_block_pos + j_mini_block*mb_header_widths[alphabet+1];

                //read all the precomputed ranks from the mini block metadata
                size_t tmp_val, off_set=0, pos;
                for(size_t u=0;u<alphabet;u++){
                    pos = j_mb_start+off_set;
                    memcpy(&tmp_val, &data_pointer[pos>>3], 2);
                    rank_c_j[u] += (tmp_val>> (pos & 7)) & 4095;
                    //rank_c_j[u] += bwt.read(j_mb_start + off_set, j_mb_start + off_set + 11);
                    off_set+=12;
                }

                //read the mini block byte position within the BWT block
                j_bwt_ptr += mb_header_bytes;
                j_bwt_ptr += bwt.read(j_mb_start + mb_header_widths[alphabet], j_mb_start + mb_header_widths[alphabet+1] -1);
                tmp_j += j_mini_block<< mb_width;
            }
        }

        i_bwt_ptr++;
        j_bwt_ptr++;

        uint16_t data;
        uint8_t long_run;
        sym_type symbol;

        //size_t n_runs=0;
        while(tmp_i<i){
            data = *i_bwt_ptr;
            i_bwt_ptr++;
            long_run = data & 1;//is a long or short run?
            data>>=1;

            symbol = (data & 15);//get the run symbol
            data>>=4;//get the run len

            data |= (-long_run & ((*i_bwt_ptr)<<3));
            i_bwt_ptr+=long_run;
            data++;

            rank_c_i[symbol]+=data;
            tmp_i+=data;
            //n_runs++;
        }
        rank_c_i[symbol]-=tmp_i-i;


        if(same_mini_block){
            memcpy(rank_c_j.data(), rank_c_i.data(), alphabet*sizeof(size_t));
            rank_c_j[symbol]+=tmp_i-i;
            tmp_j=tmp_i;
            j_bwt_ptr=i_bwt_ptr;
        }

        while(tmp_j<=j) {
            data = *j_bwt_ptr;
            j_bwt_ptr++;
            long_run = data & 1;//is a long or short run?
            data>>=1;

            symbol = (data & 15);//get the run symbol
            data>>=4;//get the run len

            data |= (-long_run & ((*j_bwt_ptr)<<3));
            j_bwt_ptr+=long_run;
            data++;

            rank_c_j[symbol]+=data;
            tmp_j+=data;
            //n_runs++;
        }
        rank_c_j[symbol]-=tmp_j-j-1;*/

        k=0;
        for(size_t u=0;u<alphabet;u++){
            if(rank_c_i[u]<rank_c_j[u]){
                rank_c_i[k] = rank_c_i[u];
                rank_c_j[k] = rank_c_j[u];
                cs[k] = sym_inv_map[u];
                k++;
            }
        }
    }

    [[nodiscard]] inline size_t select(size_t rank, sym_type symbol) const {
        assert(rank>0);

        size_t select_ans=0;
        symbol = sym_map[symbol];
        auto n_blocks = (long int)blocks();

        long int l =0, r = n_blocks-1, m;
        size_t tmp_a, tmp_b;

        while(l<=r){
            m = l + ((r - l)>>1);

            tmp_a = bwt.read(block_pointers[m] + b_header_widths[symbol],
                             block_pointers[m] + b_header_widths[symbol+1] -1);
            tmp_b = bwt.read(block_pointers[m+1] + b_header_widths[symbol],
                             block_pointers[m+1] + b_header_widths[symbol+1] -1);

            if(tmp_a<rank && rank<=tmp_b) break;

            int mask = (tmp_a<rank == 0) - 1;
            l = (l & ~mask) | ((m + 1) & mask);
            r = ((m - 1) & ~mask) | (r & mask);
        }

        //corner case: 'symbol' has less than 'rank' occurrences
        if(l>r) return n_symbols;

        size_t block_pos =  block_pointers[m]+b_header_bits;
        auto *bwt_ptr = data_pointer + (block_pos>>3);
        bool has_mini_blocks = *bwt_ptr;
        block_pos+=8;
        bwt_ptr++;
        select_ans += m<<12;

        rank-=tmp_a;

        if(has_mini_blocks) {

            auto n_m_blocks = (long int) (m<int64_t(blocks()-1)? n_mini_blocks : INT_CEIL((n_symbols-select_ans), mb_size));

            l = 0, r = n_m_blocks-1;
            size_t mb_start_a, mb_start_b;

            while(l<=r){
                m = l + ((r - l)>>1);

                mb_start_a = block_pos + m*(mb_header_widths[alphabet+1]+1);
                tmp_a = bwt.read(mb_start_a + mb_header_widths[symbol], mb_start_a + mb_header_widths[symbol+1] -1);

                mb_start_b = block_pos + (m+1)*(mb_header_widths[alphabet+1]+1);
                tmp_b = bwt.read(mb_start_b + mb_header_widths[symbol], mb_start_b + mb_header_widths[symbol+1] -1);

                if(tmp_a<rank && rank<=tmp_b) break;

                int mask = (tmp_a<rank == 0) - 1;
                l = (l & ~mask) | ((m + 1) & mask);
                r = ((m - 1) & ~mask) | (r & mask);
            }

            //read the mini block position within the block and the mini block encoding
            size_t mini_block_pos = bwt.read(mb_start_a + mb_header_widths[alphabet],
                                             mb_start_a + mb_header_widths[alphabet+1]);

            bool one_byte_encoding = !(mini_block_pos & 4096);
            mini_block_pos &= 4095;//clear the bit indicating the mini block encoding
            bwt_ptr+= mb_header_bytes + mini_block_pos;
            select_ans += m<<mb_width;

            rank-=tmp_a;

            uint16_t b_freq[16]={0};
            if(one_byte_encoding){
                uint8_t  len;
                while(b_freq[symbol]<rank){
                    len = (*bwt_ptr>>4)+1;
                    b_freq[*bwt_ptr & 15]+=len;
                    select_ans+=len;
                    bwt_ptr++;
                }
            }else{
                bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
                auto *run = (uint16_t *) bwt_ptr;
                uint16_t len;
                while(b_freq[symbol]<rank){
                    len = (*run>>4) + 1;
                    b_freq[*run & 15]+=len;
                    select_ans+=len;
                    run++;
                }
            }
            select_ans-=b_freq[symbol]-rank+1;
        }else{
            bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
            auto *run = (uint16_t *) bwt_ptr;
            uint16_t b_freq[16]={0};
            uint16_t len;
            while(b_freq[symbol]<rank){
                len = (*run>>4) + 1;
                b_freq[*run & 15]+=len;
                select_ans+=len;
                run++;
            }
            select_ans-=b_freq[symbol]-rank+1;
        }
        return select_ans;
    }

    inline sym_type operator[](size_t idx) const {

        sym_type symbol;
        size_t block = idx>>12;
        size_t block_pos = block_pointers[block];
        block_pos += b_header_bits;

        auto *bwt_ptr = data_pointer + (block_pos>>3);
        bool has_mini_blocks = *bwt_ptr;
        block_pos+=8;
        bwt_ptr++;
        size_t tmp_idx = block<<12;

        //the current position indicates if the block has mini blocks or not
        if(has_mini_blocks) {
            //idx of the mini block within the block
            size_t mini_block = (idx-tmp_idx) >> mb_width;

            //read metadata of the mini block
            size_t mb_start = block_pos + mini_block*(mb_header_widths[alphabet+1]+1);

            //read the mini block position within the block and the mini block encoding
            size_t mini_block_pos = bwt.read(mb_start + mb_header_widths[alphabet], mb_start + mb_header_widths[alphabet+1]);
            bool one_byte_encoding = !(mini_block_pos & 4096);
            mini_block_pos &= 4095;//clear the bit indicating the mini block encoding

            tmp_idx += mini_block<< mb_width;

            if((idx-tmp_idx)>(mb_size>>1)){

                if(mini_block<(n_mini_blocks-1)){
                    size_t next_mb_mt_start = block_pos + (mini_block+1)*(mb_header_widths[alphabet+1]+1);
                    mini_block_pos = bwt.read(next_mb_mt_start + mb_header_widths[alphabet],
                                              next_mb_mt_start + mb_header_widths[alphabet+1]-1)-1;
                }else{
                    mini_block_pos = block_pointers[block+1]-block_pointers[block];
                    mini_block_pos -= b_header_bits + (mb_header_bytes<<3)+8+8;//subtract the header bits
                    mini_block_pos >>=3;//in bytes
                }

                bwt_ptr+= mb_header_bytes + mini_block_pos;
                if(one_byte_encoding){
                    b_scan<true, false>(idx, tmp_idx+mb_size, bwt_ptr, nullptr, symbol);
                }else{
                    b_scan<false, false>(idx, tmp_idx+mb_size, bwt_ptr, nullptr, symbol);
                }

            }else{
                bwt_ptr+= mb_header_bytes + mini_block_pos;
                if(one_byte_encoding){
                    f_scan<true, false>(idx, tmp_idx, bwt_ptr, nullptr, symbol);
                }else{
                    bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
                    f_scan<false, false>(idx, tmp_idx, bwt_ptr, nullptr, symbol);
                }
            }
        } else {
            if((idx-tmp_idx)>(b_size>>1)){
                block_pos = block_pointers[block+1]-block_pointers[block];
                block_pos -= b_header_bits + 16;//subtract the header bits 8+8=16 consider the mini block flag and moving one position back from the end
                block_pos >>=3;//transform to byte position
                bwt_ptr+=block_pos;
                b_scan<false, false>(idx, tmp_idx+b_size, bwt_ptr, nullptr, symbol);
            }else{
                bwt_ptr+=reinterpret_cast<uintptr_t>(bwt_ptr) & 1;//move to the next two-byte-aligned position
                f_scan<false, false>(idx, tmp_idx, bwt_ptr, nullptr, symbol);
            }
        }

        return sym_inv_map[symbol];
    }

    void stats() const {
        std::cout<<"Number of runs before:        "<<orig_n_runs<<std::endl;
        std::cout<<"Number of runs now:           "<<tot_runs<<" ("<<100*((double(tot_runs)/double(orig_n_runs))-1)<<"% increase)"<<std::endl;
        std::cout<<"Number of blocks:             "<<blocks()<<std::endl;
        std::cout<<"Block space overhead:         "<<double(INT_CEIL(blocks()*b_header_bits, 8))/1000000<<" Mb"<<std::endl;
        std::cout<<"% of blocks with mini blocks: "<<100*(double(n_sampled_blocks)/double(blocks()))<<"%"<<std::endl;
        std::cout<<"Total number of mini blocks:  "<<mini_blocks()<<std::endl;
        std::cout<<"Mini blocks space overhead:   "<<double(n_sampled_blocks*mb_header_bytes)/1000000<<" Mb"<<std::endl;
        std::cout<<"BWT space usage:              "<<double(bwt.stream_size*sizeof(size_t))/1000000<<" Mb"<<std::endl;
    }

    [[nodiscard]] inline size_t size() const {
        return n_symbols;
    }

    [[nodiscard]] inline size_t mini_blocks() const {
        return n_sampled_blocks*n_mini_blocks;
    }

    [[nodiscard]] inline size_t blocks() const {
        return INT_CEIL(n_symbols, b_size);
    }

    size_t serialize(std::ofstream & ofs){
        size_t written_bytes = 0;
        written_bytes += serialize_elm(ofs, alphabet);
        written_bytes += serialize_elm(ofs, b_header_bits);
        written_bytes += serialize_elm(ofs, n_symbols);
        written_bytes += serialize_elm(ofs, mb_header_bytes);
        written_bytes += serialize_elm(ofs, tot_runs);
        written_bytes += serialize_elm(ofs, orig_n_runs);
        written_bytes += serialize_elm(ofs, n_sampled_blocks);

        written_bytes += serialize_plain_vector(ofs, sym_map);
        written_bytes += serialize_plain_vector(ofs, sym_inv_map);
        written_bytes += serialize_plain_vector(ofs, b_header_widths);
        written_bytes += serialize_plain_vector(ofs, block_pointers);
        written_bytes += bwt.serialize(ofs);

        return  written_bytes;
    }

    void load(std::ifstream & ifs){
        load_elm(ifs, alphabet);
        load_elm(ifs, b_header_bits);
        load_elm(ifs, n_symbols);
        load_elm(ifs, mb_header_bytes);
        load_elm(ifs, tot_runs);
        load_elm(ifs, orig_n_runs);
        load_elm(ifs, n_sampled_blocks);

        load_plain_vector(ifs, sym_map);
        load_plain_vector(ifs, sym_inv_map);
        load_plain_vector(ifs, b_header_widths);
        load_plain_vector(ifs, block_pointers);
        bwt.load(ifs);

        data_pointer = (uint8_t *)bwt.stream;
    }
};
#endif //SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
