//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef BWT_DTS_BENCHMARKS_DY_BL_RLBWT_H
#define BWT_DTS_BENCHMARKS_DY_BL_RLBWT_H

using block_type = std::vector<std::pair<uint32_t, uint64_t>>;

enum IFORMAT{
   GRL_BWT=0,
   RL_PLAIN=1,
   PLAIN=2
};

template<class bwt_dt_type>
void create_subtree(block_type& block, size_t b_size, size_t s_factor, size_t b_runs,
                    std::vector<uint64_t>& sym_pointers, bwt_dt_type& bwt_rep);

template<class bwt_dt_type>
void create_leaf(std::vector<block_type>& block_sequence, size_t n_blocks, std::vector<uint64_t>& sym_ptr, bwt_dt_type& bwt_rep);

template<class bwt_dt_type>
struct c_state{//state of the compression

    size_t bk_len=0;//length of the active block
    size_t bk_id=0;//id of the active block
    size_t acc_runs=0;//sum of the runs in the active blocks

    const size_t b_size;//size of each block
    const size_t s_factor;//shrinking factor for further subdivision
    const size_t b_runs;//maximum number of runs in a sequence of blocks
    std::vector<block_type>  active_blocks;//the run-length compressed c_blocks
    bwt_dt_type& bwt_rep;//data structure encoding the representation
    std::vector<uint64_t>& sym_pointers;

    explicit c_state(size_t _b_size, size_t _s_factor, size_t _b_runs,
                     std::vector<uint64_t>& _sym_pointers, bwt_dt_type& _bwt_rep): b_size(_b_size),
                                                                                   s_factor(_s_factor),
                                                                                   b_runs(_b_runs),
                                                                                   active_blocks(b_runs),
                                                                                   sym_pointers(_sym_pointers),
                                                                                   bwt_rep(_bwt_rep){}

    inline void process_run(const size_t& sym, size_t& len){

        while(len>0){
            if((bk_len+len)<b_size){
                //the run fits the block size
                assert(len<=b_size);
                bk_len+=len;
                active_blocks[bk_id].emplace_back(sym, len);
                len=0;
            } else {// we complete a new block

                //last run of the active block
                size_t split_run_len = b_size-bk_len;
                if(split_run_len!=0){
                    active_blocks[bk_id].emplace_back(sym, split_run_len);
                }
                acc_runs+=active_blocks[bk_id++].size();

                //we compete a new block sequence
                if(acc_runs>b_runs){
                    if(bk_id==1){//only one block in the sequence and it exceeds the limit of runs
                        create_subtree(active_blocks[0], b_size, s_factor, b_runs, sym_pointers, bwt_rep);//recursive partitioning
                        active_blocks[0].clear();
                        bk_id=0;
                        acc_runs=0;
                    } else {//many blocks in the block sequence

                        create_leaf(active_blocks, bk_id-1, sym_pointers, bwt_rep);//we do not consider the last block in the sequence

                        active_blocks[0].swap(active_blocks[bk_id-1]);//move the last block to the beginning of the sequence
                        for(size_t i=1;i<bk_id;i++){//clear the other blocks in the sequence
                            active_blocks[i].clear();
                        }

                        //new block sequence
                        acc_runs = active_blocks[0].size();
                        bk_id=1;

                        if(active_blocks[0].size()>b_runs){//process the remaining block it already exceeds the limit
                            create_subtree(active_blocks[0], b_size, s_factor, b_runs, sym_pointers, bwt_rep);
                            active_blocks[0].clear();
                            bk_id=0;
                            acc_runs=0;
                        }
                    }
                }

                len -=split_run_len;
                bk_len=0;
            }
        }
    }

    ~c_state(){
        //handle the last sequence of blocks
        assert(bk_len<=b_size);
        acc_runs+=active_blocks[bk_id].size();
        bk_id+=!active_blocks[bk_id].empty();

        if(bk_len==0 && acc_runs==0) return;

        if(acc_runs<=b_runs){//the last block sequence has less than the maximum number of allowed runs
            create_leaf(active_blocks, bk_id, sym_pointers, bwt_rep);
        }else{//the last block sequence exceeds the maximum number of allowed runs

            if(bk_id==1){//there is only one block in the sequence, and it exceeds the limit of runs. Break it recursively
                create_subtree(active_blocks[0], b_size, s_factor, b_runs, sym_pointers, bwt_rep);
            } else {//multiple blocks in the sequence, and the sum of their runs exceed the maximum number of runs

                //we do not consider the last block in the sequence
                create_leaf(active_blocks, bk_id-1, sym_pointers, bwt_rep);

                //handle the last block
                if(active_blocks[bk_id-1].size()>b_runs){
                    create_subtree(active_blocks[bk_id-1], b_size, s_factor, b_runs, sym_pointers, bwt_rep);
                }else{
                    active_blocks[0].swap(active_blocks[bk_id-1]);
                    create_leaf(active_blocks, 1, sym_pointers, bwt_rep);
                }
            }
        }
    }
};

template<class bwt_dt_type>
void create_subtree(block_type& block, size_t b_size, size_t s_factor, size_t b_runs,
                    std::vector<uint64_t>& sym_pointers, bwt_dt_type& bwt_rep) {

    b_size = b_size/s_factor;
    assert(b_size>=b_runs);

    size_t sym, len;
    c_state b_state(b_size, s_factor, b_runs, sym_pointers, bwt_rep);
    for(auto & run : block){
        sym = run.first;
        len = run.second;
        b_state.process_run(sym, len);
    }
}

template<class bwt_dt_type>
void create_leaf(std::vector<block_type>& block_sequence, size_t n_blocks, std::vector<uint64_t>& sym_ptr, bwt_dt_type& bwt_rep){

    size_t sym, len;
    assert(n_blocks>0);

    for(size_t j=0;j<(block_sequence[0].size()-1);j++){
        sym = block_sequence[0][j].first;
        len = block_sequence[0][j].second;

        //TODO insert the run
        sym_ptr[sym]+=len;
    }
    bwt_rep.eff_runs+=block_sequence[0].size()-1;

    for(size_t i=1;i<n_blocks;i++){

        sym = block_sequence[i-1].back().first;
        len = block_sequence[i-1].back().second;

        if(sym==block_sequence[i][0].first){
            //collapse the runs as they are the same
            block_sequence[i][0].second +=len;
        } else {
            //process the run of the previous block
            //TODO insert the run
            sym_ptr[sym]+=len;
            bwt_rep.eff_runs++;
        }

        for(size_t j=0;j<(block_sequence[i].size()-1);j++){
            sym = block_sequence[i][j].first;
            len = block_sequence[i][j].second;
            //TODO insert the run
            sym_ptr[sym]+=len;
        }
        bwt_rep.eff_runs+=block_sequence[i].size()-1;
    }

    //process the last run
    sym = block_sequence[n_blocks-1].back().first;
    len = block_sequence[n_blocks-1].back().second;
    //TODO insert the run
    sym_ptr[sym]+=len;
    bwt_rep.eff_runs++;
}

template<class bwt_dt_type>
void build_from_grl_bwt(bwt_dt_type& bwt_rep, std::string& bwt_file){

    size_t b_runs = bwt_dt_type::max_block_runs;
    size_t b_size = bwt_dt_type::block_size;
    size_t s_factor = bwt_dt_type::scale_factor;

    bwt_buff_reader bwt_buff(bwt_file);
    size_t n_runs = bwt_buff.size();
    bwt_rep.orig_runs = n_runs;

    size_t sym, len;
    std::vector<uint64_t> sym_freqs(256, 0);
    std::vector<uint8_t> sym_map(256, 0);
    for(size_t i=0;i<n_runs;i++){
        bwt_buff.read_run(i, sym, len);
        sym_freqs[sym]+=len;
    }

    size_t sigma=0;
    for(size_t i=0;i<sym_freqs.size();i++){
        if(sym_freqs[i]!=0){
            sym_map[i] = sigma;
            sym_freqs[sigma++] = sym_freqs[i];
        }
    }
    sym_freqs.resize(sigma+1);
    size_t acc=0, tmp;
    for(size_t i=0;i<sigma;i++){
        tmp = sym_freqs[i];
        sym_freqs[i]=acc;
        acc+=tmp;
    }
    sym_freqs[sigma] = tmp;

    for(size_t s=0;s<sigma;s++){
        std::cout<<s<<" "<<sym_freqs[s]<<std::endl;
    }

    {
        c_state b_state(b_size, s_factor, b_runs, sym_freqs, bwt_rep);
        for(size_t i=0;i<n_runs;i++){
            bwt_buff.read_run(i, sym, len);
            b_state.process_run(sym_map[sym], len);
        }
    }
    std::cout<<bwt_rep.orig_runs<<" "<<bwt_rep.eff_runs<<" -> "<<float(bwt_rep.orig_runs)/float(bwt_rep.eff_runs)<<std::endl;
    for(size_t s=0;s<sigma;s++){
        std::cout<<s<<" "<<sym_freqs[s]<<std::endl;
    }
}

template<class bwt_dt_type>
void build_from_rl_plain(bwt_dt_type& bwt_rep, std::string& file){
}

template<class bwt_dt_type>
void build_from_plain(bwt_dt_type& bwt_rep, std::string& file){
}

template<class bwt_dt_type, class run_type>
void build_from_vector(std::vector<run_type>& rl_bwt){
}

template<class bwt_dt_type>
void build_dybl(bwt_dt_type& bwt, std::string& bwt_file, IFORMAT f){
    switch (f) {
        case IFORMAT::GRL_BWT:
            build_from_grl_bwt<bwt_dt_type>(bwt, bwt_file);
            break;
        case IFORMAT::RL_PLAIN:
            build_from_rl_plain<bwt_dt_type>(bwt, bwt_file);
            break;
        case IFORMAT::PLAIN:
            build_from_plain<bwt_dt_type>(bwt, bwt_file);
            break;
        default:
            std::cout<<"Unknown format"<<std::endl;
            exit(1);
    }
}


template<size_t b_size, size_t s_factor, size_t b_runs>
struct rlbwt_dybl{

    static constexpr size_t block_size = b_size;
    static constexpr size_t scale_factor = s_factor;
    static constexpr size_t max_block_runs = b_runs;

    size_t orig_runs=0;
    size_t eff_runs=0;

    //TODO static asserts in block_size, scale_factor, and b_runs
};

#endif //BWT_DTS_BENCHMARKS_DY_BL_RLBWT_H
