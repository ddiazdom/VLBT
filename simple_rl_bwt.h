//
// Created by Diaz, Diego on 19.1.2024.
//

#ifndef SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
#define SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
#include "bitstream.h"
#include "bwt_io.h"
#include <vector>
#include <queue>

struct simple_rl_bwt{

    typedef uint8_t sym_type;

    static const size_t b_size = 4096;
    static const size_t mb_size = 256;
    static constexpr uint8_t mb_width = 8;

    static const size_t max_runs_per_block = 128;
    static constexpr size_t max_run_len = 2048;//max value we can encode in 11 bits
    static const size_t n_mini_blocks = INT_CEIL(b_size, mb_size);
    static constexpr uint8_t mb_header_widths[16] = {0, 12, 24, 36, 48, 60, 72, 84, 96,
                                                     108, 120, 132, 144, 156, 168, 180};

    size_t alphabet=0;
    size_t b_header_bits=0;
    size_t n_symbols=0;
    size_t mb_header_bytes=0;
    size_t tot_runs=0;
    size_t orig_n_runs=0;
    size_t n_sampled_blocks=0;
    uint8_t *data_pointer;

    std::vector<uint8_t> sym_map;
    std::vector<uint8_t> sym_inv_map;
    std::vector<size_t> freq_widths;
    std::vector<size_t> block_pointers;
    bitstream<size_t> bwt;

    size_t len_hist[4096]={0};

    inline size_t insert_run(size_t sym, size_t len, size_t& bwt_pos) {

        //short run encoded in one byte
        assert(len>0);
        len_hist[len]++;

        size_t run_byte_pos=0;
        if(len<=8){
            uint8_t run = ((len-1)<< 4) | sym;
            run = (run<<1) | 0;
            //bwt.write(bwt_pos, bwt_pos + 8 -1, run);
            data_pointer[bwt_pos>>3] = run;
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
        tot_runs++;
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
                if(broken_run_len>0){
                    runs_byte_pos +=insert_run(run.first, broken_run_len, bwt_pos);
                    b_sym_freqs[run.first]+=broken_run_len;
                }

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
        sym_map.resize(sym_inv_map.back()+1);

        size_t acc=0, tmp;
        freq_widths.resize(alphabet+1);
        for(size_t i=0;i<alphabet;i++){
            tmp = freq_widths[i];
            freq_widths[i] = acc;
            acc +=tmp;
        }
        freq_widths[alphabet] = acc;
        b_header_bits = INT_CEIL(acc, 8)*8;

        orig_n_runs = bwt_buff.size();

        size_t n_blocks = INT_CEIL(n_syms, b_size);
        block_pointers.resize(n_blocks, 0);

        // number of bytes used by the concatenated rank samples of the mini blocks within a block.
        // the alphabet+1 position within the mini block header stores the byte position of the mini block within
        // the block. This position requires 12 bits: in the extreme case that all the runs use one byte, then
        // the highest value is 4096 bytes. The other extreme case is that all the runs use two bytes: then we have
        // at most ceil(4096/9)=456 runs and thus 456*2=912 bytes
        mb_header_bytes = INT_CEIL((n_mini_blocks*mb_header_widths[alphabet+1]), 8);//align metadata to bytes

        //estimate the number of bits in the BWT
        size_t bwt_size_bits = n_blocks*(b_header_bits+8) + n_sampled_blocks*(mb_header_bytes*8) + eff_runs*16;
        bwt.stream_size = INT_CEIL(bwt_size_bits, (sizeof(size_t)*8));
        bwt.stream = (size_t *) malloc(bwt.stream_size*sizeof(size_t));
        data_pointer = (uint8_t *)bwt.stream;

        size_t bwt_pos = 0;
        size_t idx_block=0;

        std::vector<size_t> acc_ranks(alphabet, 0);

        block_pointers[idx_block++] = bwt_pos;
        insert_block_header(acc_ranks, freq_widths, b_header_bits, bwt_pos);

        acc_block=0;
        std::vector<std::pair<uint8_t, uint16_t>> block_runs;
        block_runs.reserve(b_size);

        for(size_t k=0;k<n_runs;k++){

            bwt_buff.read_run(k, sym, len);

            sym = sym_map[sym];
            if((acc_block+len)>b_size){

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

                    //bool alph[16]={false};
                    for(auto const& run : block_runs){
                        insert_run(run.first, run.second, bwt_pos);
                        //alph[run.first] = true;
                    }

                    /*size_t eff_alph=0;
                    for(size_t i=0;i<16;i++){
                        eff_alph+=alph[i];
                    }
                    size_t alph_bits = sym_width(eff_alph);
                    size_t max_len = 1 << (8-alph_bits);
                    size_t comp_runs=0;
                    for(auto const& run : block_runs){
                        std::cout<<run.second<<" "<<eff_alph<<" "<<alph_bits<<" "<<(run.second<=max_len && run.second>8)<<std::endl;
                        comp_runs+=run.second<=max_len && run.second>8;
                    }
                    std::cout<<double(comp_runs)/double(block_runs.size())<<" "<<comp_runs<<" out of "<<block_runs.size()<<" "<<eff_alph<<std::endl;*/
                }

                //break rules into blocks as long as they are bigger than the block size
                broken_run_len = (acc_block+len)-b_size;
                while(broken_run_len>b_size){
                    block_pointers[idx_block++] = bwt_pos;
                    insert_block_header(acc_ranks, freq_widths, b_header_bits, bwt_pos);
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
                insert_block_header(acc_ranks, freq_widths, b_header_bits, bwt_pos);

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
        data_pointer = (uint8_t *)bwt.stream;

        /*size_t acc_l=0;
        for(size_t i=0;i<4096;i++){
            if(len_hist[i]!=0){
                acc_l+=len_hist[i];
                std::cout<<i<<" "<<len_hist[i]<<" "<<double(acc_l)/double(tot_runs)<<std::endl;
            }
        }*/
    }

    [[nodiscard]] inline std::pair<size_t, sym_type> inverse_select(size_t idx) const {

        uint16_t b_freq[16] = {0};

        size_t block = idx>>12;
        size_t block_pos = block_pointers[block];

        size_t b_start =  block_pos;
        size_t mb_start;

        block_pos += b_header_bits;

        auto *bwt_ptr = data_pointer + (block_pos>>3);
        bool has_mini_blocks = *bwt_ptr;
        block_pos+=8;
        bwt_ptr++;

        size_t tmp_idx = block<<12;
        size_t mini_block=0;

        //the current position indicates if the block was sub sampled or not
        if(has_mini_blocks){
            //idx of the mini block within the block
            mini_block = (idx-tmp_idx) >> mb_width;

            //read metadata of the mini block
            mb_start = block_pos + mini_block*mb_header_widths[alphabet+1];

            //read the mini block position within the block
            bwt_ptr += mb_header_bytes + bwt.read(mb_start + mb_header_widths[alphabet], mb_start + mb_header_widths[alphabet+1] -1);
            tmp_idx += mini_block<< mb_width;
        }

        uint16_t data;
        uint8_t long_run, symbol;

        while(tmp_idx<=idx){

            data = *bwt_ptr;
            bwt_ptr++;

            long_run = data & 1;//is a long or short run?
            data>>=1;

            symbol = (data & 15);//run symbol
            data>>=4;//run len

            //visit the next bit if the run uses two bytes
            data |= (-long_run & ((*bwt_ptr)<<3));
            bwt_ptr+=long_run;

            //the len is in data (adding one due to encoding)
            data++;

            b_freq[symbol]+=data;
            tmp_idx+=data;
        }

        size_t rank = bwt.read(b_start + freq_widths[symbol], b_start + freq_widths[symbol+1] -1);
        if(mini_block>0){
            rank+= bwt.read(mb_start + mb_header_widths[symbol], mb_start + mb_header_widths[symbol+1] -1);
        }
        rank+=b_freq[symbol];

        return {rank-(tmp_idx-idx), sym_inv_map[symbol]};
    }

    [[nodiscard]] inline size_t rank(size_t idx, sym_type symbol) const {

        symbol = sym_map[symbol];
        size_t block = idx>>12;
        size_t block_pos = block_pointers[block];

        size_t rank = bwt.read(block_pos+freq_widths[symbol], block_pos+freq_widths[symbol+1]-1);
        block_pos += b_header_bits;

        auto *bwt_ptr = (uint8_t *)bwt.stream;
        bwt_ptr += (block_pos>>3);
        block_pos+=8;

        //bool sub_sampled_block = *bwt_ptr;
        size_t tmp_idx = block<<12;

        //the current position indicates if the block was sub sampled or not
        if(*bwt_ptr){
            //idx of the mini block within the block
            size_t mini_block = (idx-tmp_idx) >> mb_width;

            //read metadata of the mini block
            size_t mt_start = block_pos + mini_block*mb_header_widths[alphabet+1];

            //read the rank of the queried symbol
            rank+= bwt.read(mt_start + mb_header_widths[symbol], mt_start + mb_header_widths[symbol+1] -1);

            //read the mini block position within the block
            size_t mini_block_pos = bwt.read(mt_start + mb_header_widths[alphabet], mt_start + mb_header_widths[alphabet+1] -1);
            bwt_ptr += mb_header_bytes + mini_block_pos;
            tmp_idx += mini_block<< mb_width;
        }
        bwt_ptr++;

        uint16_t data;
        uint8_t long_run;
        int equal;
        while(tmp_idx<idx){
            data = *bwt_ptr;
            bwt_ptr++;
            long_run = data & 1;//is a long or short run?
            data>>=1;

            equal = (data & 15)==symbol;//run symbol matches query symbol?
            data>>=4;//get the run len

            data |= (-long_run & ((*bwt_ptr)<<3));
            bwt_ptr+=long_run;
            data++;

            rank+= (-equal & data);
            tmp_idx+=data;
        }
        return rank - (-equal & (tmp_idx-idx));
    }

    inline void interval_symbols(size_t i, size_t j, size_t& k,
                                 std::vector<sym_type>& cs,
                                 std::vector<size_t>& rank_c_i,
                                 std::vector<size_t>& rank_c_j) const {

        if(j-i==0){
            k=0;
        }else if(j-i==1){
            auto res = inverse_select(i);
            k=1;
            cs[0] = res.second;
            rank_c_i[0] = res.first;
            rank_c_j[0] = res.first+1;
        }else{
            //auto t1 = std::chrono::high_resolution_clock::now();

            size_t i_block = i>>12;
            j--;
            size_t j_block = j>>12;

            size_t i_block_pos = block_pointers[i_block];
            size_t j_block_pos;

            //sample the ranks up to the respective blocks
            if(j_block==i_block){
                //copy the ranks if i and j-1 are withing the same block
                j_block_pos = i_block_pos;
                for(size_t u=0;u<alphabet;u++){
                    rank_c_i[u] = bwt.read(i_block_pos+freq_widths[u], i_block_pos+freq_widths[u+1]-1);
                }
                memcpy(rank_c_j.data(), rank_c_i.data(), alphabet*sizeof(size_t));
            }else{
                j_block_pos = block_pointers[j_block];
                for(size_t u=0;u<alphabet;u++){
                    rank_c_i[u] = bwt.read(i_block_pos+freq_widths[u], i_block_pos+freq_widths[u+1]-1);
                    rank_c_j[u] = bwt.read(j_block_pos+freq_widths[u], j_block_pos+freq_widths[u+1]-1);
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

            //auto t2 = std::chrono::high_resolution_clock::now();
            //std::cout<<"The first part: "<<std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count()<<" now we will scan "<<i-tmp_i<<" symbols from "<<tmp_i<<" to "<<i<<" block: "<<i_block<<" and mini block: "<<i_mini_block<<std::endl;


            //auto t1 = std::chrono::high_resolution_clock::now();

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
            rank_c_j[symbol]-=tmp_j-j-1;

            //auto t2 = std::chrono::high_resolution_clock::now();
            //std::cout<<"The second part: "<<std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count()<<" "<<i<<" "<<j<<" "<<n_runs<<" "<<same_mini_block<<std::endl;

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
    }

    sym_type operator[](size_t idx){

        sym_type symbol;
        size_t block = idx>>12;
        size_t block_pos = block_pointers[block] + b_header_bits;

        auto *bwt_ptr = (uint8_t *)bwt.stream;
        bwt_ptr += (block_pos>>3);
        block_pos+=8;

        size_t tmp_idx = block<<12;

        //the current position indicates if the block was sub sampled or not
        if(*bwt_ptr){
            //idx of the mini block within the block
            size_t mini_block = (idx-tmp_idx) >> mb_width;

            //metadata of the mini block
            size_t mt_start = block_pos + mini_block*mb_header_widths[alphabet+1];

            //read the mini block position within the block
            size_t mini_block_pos = bwt.read(mt_start + mb_header_widths[alphabet], mt_start + mb_header_widths[alphabet+1] -1);
            bwt_ptr += mb_header_bytes + mini_block_pos;
            tmp_idx += mini_block<< mb_width;
        }
        bwt_ptr++;

        uint16_t data;
        uint8_t long_run;
        while(tmp_idx<=idx){
            data = *bwt_ptr;
            bwt_ptr++;
            long_run = data & 1;//is a long or short run?
            data>>=1;

            symbol = (data & 15);//run symbol matches query symbol?
            data>>=4;//get the run len

            data |= (-long_run & ((*bwt_ptr)<<3));
            bwt_ptr+=long_run;
            data++;

            tmp_idx+=data;
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

    [[nodiscard]] size_t size() const {
        return n_symbols;
    }

    [[nodiscard]] size_t mini_blocks() const {
        return n_sampled_blocks*n_mini_blocks;
    }

    [[nodiscard]] size_t blocks() const {
        return INT_CEIL(n_symbols, b_size);
    }

    size_t serialize(std::ofstream & ofs){
        size_t written_bytes = 0;
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
        //load_elm(ifs, n_blocks);
        load_elm(ifs, alphabet);
        load_elm(ifs, b_size);
        load_plain_vector(ifs, sym_map);
        load_plain_vector(ifs, sym_inv_map);
        load_plain_vector(ifs, freq_widths);
        load_plain_vector(ifs, block_pointers);
        bwt.load(ifs);
        data_pointer = (uint8_t *)bwt.stream;
    }
};
#endif //SIMPLE_RL_BWT_2024_SIMPLE_RL_BWT_H
