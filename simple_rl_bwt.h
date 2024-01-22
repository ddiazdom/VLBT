//
// Created by Diaz, Diego on 19.1.2024.
//

#ifndef SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
#define SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
#include "bitstream.h"
#include "bwt_io.h"
#include <cstdlib>
#include <vector>
#include <queue>

struct simple_rl_bwt{

    typedef uint8_t sym_type;

    static const size_t b_size = 4096;
    static const size_t max_runs_per_block = 256;
    static const size_t mb_size = 512;
    static constexpr size_t max_run_len = 2048;
    static const size_t n_mini_blocks = INT_CEIL(b_size, mb_size);

    size_t n_blocks=0;
    size_t alphabet=0;
    size_t header_width=0;
    size_t n_symbols=0;
    size_t mb_header_bytes=0;

    std::vector<uint8_t> sym_map;
    std::vector<uint8_t> sym_inv_map;
    std::vector<size_t> freq_widths;
    std::vector<size_t> block_pointers;
    bitstream<size_t> bwt;

    inline size_t insert_run(size_t sym, size_t len, size_t& bwt_pos) {

        //short run encoded in one byte
        size_t run_byte_pos=0;
        if(len<=8){
            uint8_t run = ((len-1)<< 4) | sym;
            run = (run<<1) | 0;
            bwt.write(bwt_pos, bwt_pos + 8 -1, run);
            bwt_pos +=8;
            run_byte_pos+=1;
        }else{
            //long run encoded in two bytes
            uint16_t run = ((len-1) << 4) | sym;
            run = (run<<1) | 1;
            bwt.write(bwt_pos, bwt_pos + 16 -1, run);
            bwt_pos +=16;
            run_byte_pos+=2;
        }
        return run_byte_pos;
    }

    inline void subsample_block(std::vector<std::pair<uint8_t, uint16_t>>& block_runs, size_t& bwt_pos) {

        //metadata within the block
        size_t mb_rank_offset = bwt_pos;

        //offset of the runs within the block
        bwt_pos +=mb_header_bytes*8;

        //start of the runs of mini blocks within the block
        size_t runs_byte_pos = 0;

        std::vector<size_t> b_sym_freqs(alphabet+1, 0);
        std::vector<size_t> mb_rank_widths(alphabet+2, 0);
        for(size_t i=0;i<=alphabet+1;i++) mb_rank_widths[i] = i*12;

        insert_block_header(b_sym_freqs, mb_rank_widths, 12*(alphabet+1), mb_rank_offset);

        size_t acc_block=0, broken_run_len;
        for(auto const& run : block_runs) {

            if((acc_block+run.second)>mb_size) {

                //last run of the previous mini block
                broken_run_len = mb_size-acc_block;
                runs_byte_pos +=insert_run(run.first, broken_run_len, bwt_pos);
                b_sym_freqs[run.first]+=broken_run_len;

                //break rules into mini blocks as long as they are
                // bigger than the mini block size
                broken_run_len = (acc_block+run.second)-mb_size;
                while(broken_run_len>mb_size){

                    b_sym_freqs[alphabet] = runs_byte_pos;
                    insert_block_header(b_sym_freqs, mb_rank_widths, 12*(alphabet+1), mb_rank_offset);

                    runs_byte_pos+=insert_run(run.first, mb_size, bwt_pos);
                    b_sym_freqs[run.first]+=mb_size;
                    broken_run_len-=mb_size;
                }

                //insert the header of the first mini run
                b_sym_freqs[alphabet] = runs_byte_pos;
                insert_block_header(b_sym_freqs, mb_rank_widths, 12*(alphabet+1), mb_rank_offset);

                //insert the header of the first mini run
                runs_byte_pos+=insert_run(run.first, broken_run_len, bwt_pos);
                b_sym_freqs[run.first]+= broken_run_len;
                acc_block = broken_run_len;

            } else {
                //the run fits the mini block size
                runs_byte_pos+=insert_run(run.first, run.second, bwt_pos);
                b_sym_freqs[run.first]+=run.second;
                acc_block+=run.second;
            }
        }
        assert(acc_block==mb_size);
    }

    inline void insert_block_header(std::vector<size_t>& sym_freq, std::vector<size_t>& rank_widths, size_t h_width, size_t& bwt_pos) {
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
        size_t acc_block=0, broken_run_len, eff_runs=0, runs_in_block=0, mb_overhead=0;
        for(size_t i=0;i<n_runs;i++){
            bwt_buff.read_run(i, sym, len);

            if((acc_block+len)>b_size){
                //last run of the previous block
                broken_run_len = b_size-acc_block;
                sym_freqs[sym]+=broken_run_len;
                eff_runs += INT_CEIL(broken_run_len, max_run_len);
                runs_in_block += INT_CEIL(broken_run_len, max_run_len);

                if(runs_in_block>=max_runs_per_block){
                    mb_overhead += INT_CEIL(n_mini_blocks*12*(alphabet+1),8)*8;//overhead of the block metadata
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
        freq_widths.resize(256);
        for(unsigned long sym_freq : sym_freqs){
            n_syms+=sym_freq;
            if(sym_freq>0){
                freq_widths[alphabet] = sym_width(sym_freq);
                sym_map[u] = alphabet;
                sym_inv_map[alphabet] = u;
                alphabet++;
                n_symbols+=sym_freq;
            }
            u++;
        }

        sym_inv_map.resize(alphabet);
        sym_map.resize(sym_inv_map.back());

        size_t acc=0, tmp;
        freq_widths.resize(alphabet+1);
        for(size_t i=0;i<alphabet;i++){
            tmp = freq_widths[i];
            freq_widths[i] = acc;
            acc +=tmp;
        }
        freq_widths[alphabet] = acc;
        header_width = INT_CEIL(acc, 8)*8;

        n_blocks = INT_CEIL(n_syms, b_size);
        block_pointers.resize(n_blocks, 0);

        //estimate the number of bits in the BWT
        size_t bwt_size_bits = n_blocks*(freq_widths.back()+20) + eff_runs*16 + mb_overhead;
        bwt.stream_size = INT_CEIL(bwt_size_bits, (sizeof(size_t)*8));
        bwt.stream = (size_t *) malloc(bwt.stream_size*sizeof(size_t));

        size_t bwt_pos = 0;
        size_t idx_block=0;

        std::vector<size_t> acc_ranks(alphabet, 0);

        block_pointers[idx_block++] = bwt_pos;
        insert_block_header(acc_ranks, freq_widths, header_width, bwt_pos);

        acc_block=0;
        std::vector<std::pair<uint8_t, uint16_t>> block_runs;
        block_runs.reserve(b_size);

        // number of bytes for the mini blocks' header
        size_t bits_mb_metadata = 12 * (alphabet+1) * n_mini_blocks;
        mb_header_bytes = INT_CEIL(bits_mb_metadata, 8);//align metadata to bytes

        for(size_t k=0;k<n_runs;k++){

            bwt_buff.read_run(k, sym, len);

            sym = sym_map[sym];
            if((acc_block+len)>b_size){

                //last run of the previous block
                broken_run_len = b_size-acc_block;
                acc_ranks[sym] += broken_run_len;
                assert(acc_block+broken_run_len==b_size);

                while(broken_run_len>max_run_len){
                    block_runs.emplace_back(sym, broken_run_len);
                    broken_run_len-=max_run_len;
                }
                if(broken_run_len>0) block_runs.emplace_back(sym, broken_run_len);

                if(block_runs.size()>=max_runs_per_block){
                    //mark the block as subsampled (i.e., it has mini blocks)
                    bwt.write(bwt_pos, bwt_pos+8-1, 1);
                    bwt_pos+=8;
                    subsample_block(block_runs, bwt_pos);
                }else{
                    //mark the block as not subsampled (i.e., it does not have mini blocks)
                    bwt.write(bwt_pos, bwt_pos+8-1, 0);
                    bwt_pos+=8;
                    for(auto const& run : block_runs){
                        insert_run(run.first, run.second, bwt_pos);
                    }
                }

                //break rules into blocks as long as they are bigger than the block size
                broken_run_len = (acc_block+len)-b_size;
                while(broken_run_len>b_size){
                    block_pointers[idx_block++] = bwt_pos;
                    insert_block_header(acc_ranks, freq_widths, header_width, bwt_pos);
                    //mark as not sub sampled
                    bwt.write(bwt_pos, bwt_pos+8-1, 0);
                    bwt_pos+=8;

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
                insert_block_header(acc_ranks, freq_widths, header_width, bwt_pos);

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
                if(len>0) block_runs.emplace_back(sym, len);
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
            for(auto const& run : block_runs){
                insert_run(run.first, run.second, bwt_pos);
            }
        }

        assert(idx_block == n_blocks);
        assert(bwt_pos<=bwt_size_bits);

        //shrink to fit
        bwt.stream_size = INT_CEIL(bwt_pos, (sizeof(size_t)*8));
        bwt.stream = (size_t *) realloc(bwt.stream, bwt.stream_size*sizeof(size_t));

        std::vector<size_t> tmp_ranks(alphabet, 0);
        size_t pos=0;
        for(size_t i=0;i<bwt_buff.size();i++){
            bwt_buff.read_run(i, sym, len);

            if((i % 100000)==0){
                std::cout<<" -> "<<sym<<" "<<len<<" "<<pos<<" "<<n_symbols<<std::endl;
            }


            for(size_t k=0;k<len;k++){
                bool equal = rank(pos, sym)==tmp_ranks[sym_map[sym]];
                if(!equal){
                    std::cout<<"? "<<pos<<","<<sym<<" -> "<<rank(pos, sym)<<" "<<tmp_ranks[sym_map[sym]]<<std::endl;
                    rank(pos, sym);
                }
                assert(equal);
                tmp_ranks[sym_map[sym]]++;
                pos++;
            }
        }
    }

    std::pair<size_t, uint8_t> inverse_select(size_t& idx){

    }

    [[nodiscard]] inline size_t rank(size_t idx, sym_type symbol) const {

        symbol = sym_map[symbol];
        size_t block = idx>>12;
        size_t block_pos = block_pointers[block];

        size_t rank = bwt.read(block_pos+freq_widths[symbol], block_pos+freq_widths[symbol+1]-1);
        block_pos += header_width;

        //the header is a byte aligned
        assert((block_pos & 3) ==0);

        auto *bwt_ptr = (uint8_t *)bwt.stream;
        bwt_ptr += (block_pos>>3);
        block_pos+=8;

        //bool sub_sampled_block = *bwt_ptr;
        size_t tmp_idx = block<<12;

        //the current position indicates if the block was sub sampled or not
        if(*bwt_ptr){

            //idx of the mini block within the block
            size_t mini_block = (idx-tmp_idx) >> 9;

            //read metadata of the mini block
            size_t mt_start = block_pos + mini_block*(alphabet+1)*12;

            //read the rank queried symbol
            rank+= bwt.read(mt_start + symbol*12, mt_start + (symbol+1)*12 -1);

            //read the mini block position within the block
            size_t mini_block_pos = bwt.read(mt_start + alphabet*12, mt_start + (alphabet+1)*12 -1);
            bwt_ptr += mb_header_bytes + mini_block_pos;
            tmp_idx += mini_block<<9;
            //assert(idx-tmp_idx<=mb_size);
        }
        bwt_ptr++;

        uint16_t data;
        bool long_run;
        int equal;
        while(tmp_idx<idx){

            data = *bwt_ptr;
            bwt_ptr++;
            long_run = data & 1;//is a long or short run?
            data>>=1;

            equal = (data & 15)==symbol;//run symbol matches query symbol?
            data>>=4;//get the run len

            if(long_run){//get the missing bits of the len
                data |=(*bwt_ptr)<<3;
                bwt_ptr++;
            }
            data++;

            rank+= (-equal & data) | (-!equal & 0);//add len only if it matches
            tmp_idx+=data;
        }

        size_t excess = tmp_idx-idx;
        rank-= (-equal & excess) | (-!equal & 0);

        return rank;
    }

    void interval_symbols(size_t& idx){

    }

    sym_type operator[](size_t idx){
        return 'A';
    }

    size_t serialize(std::ofstream & ofs){
        size_t written_bytes = 0;
        written_bytes += serialize_elm(ofs, n_blocks);
        written_bytes += serialize_elm(ofs, alphabet);
        written_bytes += serialize_elm(ofs, b_size);
        written_bytes += serialize_plain_vector(ofs, sym_map);
        written_bytes += serialize_plain_vector(ofs, sym_inv_map);
        written_bytes += serialize_plain_vector(ofs, freq_widths);
        written_bytes += serialize_plain_vector(ofs, block_pointers);
        written_bytes += bwt.serialize(ofs);
        return  written_bytes;
    }

    void load(std::ifstream & ifs){
        load_elm(ifs, n_blocks);
        load_elm(ifs, alphabet);
        load_elm(ifs, b_size);
        load_plain_vector(ifs, sym_map);
        load_plain_vector(ifs, sym_inv_map);
        load_plain_vector(ifs, freq_widths);
        load_plain_vector(ifs, block_pointers);
        bwt.load(ifs);
    }
};
#endif //SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
