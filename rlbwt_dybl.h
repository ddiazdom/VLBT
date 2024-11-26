//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef BWT_DTS_BENCHMARKS_DY_BL_RLBWT_H
#define BWT_DTS_BENCHMARKS_DY_BL_RLBWT_H

#include <cmath>

using block_type = std::vector<std::pair<uint32_t, uint64_t>>;
using buffer_type = std::vector<std::vector<block_type>>;

enum INPUT_FORMAT{
   GRL_BWT=0,
   RL_PLAIN=1,
   PLAIN=2
};

template<class bwt_dt_type>
struct rl_tree{//state of the compression

    size_t bk_len=0;//length of the active block
    size_t bk_id=0;//id of the active block
    size_t acc_runs=0;//sum of the runs in the active blocks
    size_t n_children=0;//number of children of this node
    size_t written_bits=0;//number of bits required for the subtree

    const size_t lvl;//level of the subtree
    const size_t b_size;//size of each block
    const size_t s_factor;//shrinking factor for further subdivision
    const size_t b_runs;//maximum number of runs in a sequence of blocks
    buffer_type& block_buff;
    std::vector<block_type>& active_blocks;//the run-length compressed c_blocks
    bwt_dt_type& bwt_rep;//data structure encoding the representation
    std::vector<uint64_t>& sym_pointers;//information that we store in the headers of every run

    explicit rl_tree(size_t _lvl, size_t _b_size, size_t _s_factor, size_t _b_runs,
                     std::vector<uint64_t>& _sym_pointers, bwt_dt_type& _bwt_rep, buffer_type& _block_buff): lvl(_lvl),
                                                                                                             b_size(_b_size),
                                                                                                             s_factor(_s_factor),
                                                                                                             b_runs(_b_runs),
                                                                                                             block_buff(_block_buff),
                                                                                                             active_blocks(block_buff[lvl]),
                                                                                                             bwt_rep(_bwt_rep),
                                                                                                             sym_pointers(_sym_pointers){
        assert(b_size>=b_runs);
    }

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
                active_blocks[bk_id].emplace_back(sym, split_run_len);
                bk_len+=split_run_len;
                assert(split_run_len>0 && bk_len==b_size);

                acc_runs+=active_blocks[bk_id++].size();

                //we compete a new block sequence
                if(acc_runs>b_runs){

                    if(bk_id==1){//only one block in the sequence and it exceeds the limit of runs
                        create_subtree(active_blocks[0]);//recursive partitioning
                        active_blocks[0].clear();
                        bk_id=0;
                        acc_runs=0;
                    } else {//many blocks in the block sequence
                        create_leaf(bk_id-1);//we do not consider the last block in the sequence
                        active_blocks[0].swap(active_blocks[bk_id-1]);//move the last block to the beginning of the sequence
                        for(size_t i=1;i<bk_id;i++){//clear the other blocks in the sequence
                            active_blocks[i].clear();
                        }

                        //new block sequence
                        acc_runs = active_blocks[0].size();
                        bk_id=1;

                        if(active_blocks[0].size()>b_runs){//process the remaining block it already exceeds the limit
                            create_subtree(active_blocks[0]);
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

    void finish_scan(){
        //handle the last sequence of blocks
        assert(bk_len<=b_size);
        acc_runs+=active_blocks[bk_id].size();
        bk_id+=!active_blocks[bk_id].empty();

        if(bk_len==0 && acc_runs==0) return;

        if(acc_runs<=b_runs){//the last block sequence has less than the maximum number of allowed runs
            create_leaf(bk_id);
        }else{//the last block sequence exceeds the maximum number of allowed runs

            if(bk_id==1){//there is only one block in the sequence, and it exceeds the limit of runs. Break it recursively
                create_subtree(active_blocks[0]);
            } else {//multiple blocks in the sequence, and the sum of their runs exceed the maximum number of runs

                //we do not consider the last block in the sequence
                create_leaf(bk_id-1);

                //handle the last block
                if(active_blocks[bk_id-1].size()>b_runs){
                    create_subtree(active_blocks[bk_id-1]);
                }else{
                    active_blocks[0].swap(active_blocks[bk_id-1]);
                    create_leaf(1);
                }
            }
        }
        for(size_t i=0;i<bk_id;i++){
            active_blocks[i].clear();
        }

        size_t node_bits = s_factor + (n_children*bwt_rep.sigma*sym_width(b_size));//the rank information
        node_bits+= sym_width(written_bits)*n_children;//the pointers for the children of this node
        written_bits+=node_bits;
    }

    inline void create_leaf(size_t n_blocks){

        n_children++;
        size_t sym, len, n_runs=0;
        assert(n_blocks>0);
        std::vector<bool> block_sigma(256, false);
        size_t longest_run = 0;

        for(size_t j=0;j<(active_blocks[0].size()-1);j++){
            sym = active_blocks[0][j].first;
            len = active_blocks[0][j].second;

            block_sigma[sym]=true;
            if(len>longest_run) longest_run = len;
            //TODO insert the run
            sym_pointers[sym]+=len;
        }
        n_runs+=active_blocks[0].size()-1;

        for(size_t i=1;i<n_blocks;i++){

            sym = active_blocks[i-1].back().first;
            len = active_blocks[i-1].back().second;

            if(sym==active_blocks[i][0].first){
                //collapse the runs as they are the same
                active_blocks[i][0].second +=len;
            } else {
                //process the run of the previous block
                //TODO insert the run
                block_sigma[sym] = true;
                if(len>longest_run) longest_run = len;
                sym_pointers[sym]+=len;
                n_runs++;
            }

            for(size_t j=0;j<(active_blocks[i].size()-1);j++){
                sym = active_blocks[i][j].first;
                len = active_blocks[i][j].second;
                //TODO insert the run
                block_sigma[sym]=true;
                if(len>longest_run) longest_run = len;
                sym_pointers[sym]+=len;
            }
            n_runs+=active_blocks[i].size()-1;
        }

        //process the last run
        sym = active_blocks[n_blocks-1].back().first;
        len = active_blocks[n_blocks-1].back().second;

        //TODO insert the run
        block_sigma[sym]=true;
        if(len>longest_run) longest_run = len;
        sym_pointers[sym]+=len;
        n_runs++;

        size_t b_sigma=0;
        for(auto && bit : block_sigma){
            b_sigma+=bit;
        }

        size_t bytes_per_run = INT_CEIL((sym_width(b_sigma)+sym_width(longest_run)), 8);

        bwt_rep.eff_runs += n_runs;
        assert(n_runs<=bwt_dt_type::max_block_runs);

        bwt_rep.r_freq[n_runs]++;
        bwt_rep.lvl_freq[lvl]++;
        bwt_rep.sigma_freq[bytes_per_run]++;

        size_t node_bits = (n_runs*bytes_per_run)*8;
        written_bits+=node_bits;
    }

    inline void create_subtree(block_type& block) {

        n_children++;

        //each child block covers up to b_runs runs, so all the runs in the subtree can use more than
        //s_factor*b_runs*(log(sigma)+log(b_size)) bits of space

        //each child block has recursive headers, where each header uses, at most,
        //sigma*log(b_size) + 1 bits, where the +1 indicates leaf or internal node.
        //the number of nodes in the substree is s_factor + s_factor^{2} + s_factor^{r}, with r being s_factor^{r}=b_runs
        //the total number internal nodes then becomes (s_factor^{r}-1)/(s_factor-1)
        //Thus, the final space is:
        // (s_factor^{r} = b_run) * (sigma*log(b_size)+1) = headers of the leaves
        // (s_factor^{r-1})/(s_factor-1) * (sigma*log(b_size)+1) = headers of internal nodes

        //new internal node in the three
        //bwt_rep.lvl_freq[lvl]++;
        rl_tree int_node(lvl+1, b_size/s_factor, s_factor, b_runs, sym_pointers, bwt_rep, block_buff);

        size_t sym, len;
        for(auto & run : block){
            sym = run.first;
            len = run.second;
            int_node.process_run(sym, len);
        }
        int_node.finish_scan();

        written_bits+=int_node.written_bits;
        bwt_rep.n_children[int_node.n_children]++;
    }
};

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

    size_t sigma=0, sigma_bits=0;
    for(size_t i=0;i<sym_freqs.size();i++){
        if(sym_freqs[i]!=0){
            sym_map[i] = sigma;
            sym_freqs[sigma++] = sym_freqs[i];
            sigma_bits+= sym_width(sym_freqs[i]);
        }
    }
    sym_freqs.resize(sigma+1);
    bwt_rep.sigma = sigma;

    size_t acc=0, tmp;
    for(size_t i=0;i<sigma;i++){
        tmp = sym_freqs[i];
        sym_freqs[i]=acc;
        acc+=tmp;
    }
    sym_freqs[sigma] = acc;
    bwt_rep.tot_syms = acc;

    for(size_t s=0;s<sigma;s++){
        std::cout<<s<<" "<<sym_freqs[s]<<std::endl;
    }

    buffer_type block_buff(bwt_rep.levels, std::vector<block_type>(b_runs));
    rl_tree root(0, b_size, s_factor, b_runs, sym_freqs, bwt_rep, block_buff);
    for(size_t i=0;i<n_runs;i++){
        bwt_buff.read_run(i, sym, len);
        root.process_run(sym_map[sym], len);
    }
    root.finish_scan();

    std::cout<<bwt_rep.orig_runs<<" "<<bwt_rep.eff_runs<<" -> "<<float(bwt_rep.orig_runs)/float(bwt_rep.eff_runs)<<std::endl;
    for(size_t s=0;s<sigma;s++){
        std::cout<<s<<" "<<sym_freqs[s]<<std::endl;
    }

    for(size_t r=0;r<=bwt_dt_type::max_block_runs;r++){
        std::cout<<r<<" : "<<bwt_rep.r_freq[r]<<std::endl;
    }

    std::cout<<"Number of tree levels: "<<bwt_rep.levels<<std::endl;
    size_t tot_leaves=0;
    for(size_t i=0;i<20;i++){
        tot_leaves+=bwt_rep.lvl_freq[i];
    }

    for(size_t i=0;i<20;i++){
        if(bwt_rep.lvl_freq[i]>0){
            std::cout<<i<<" -> "<<double(bwt_rep.lvl_freq[i])/double(tot_leaves)<<std::endl;
        }
    }

    std::cout<<"Bytes per run in a block"<<std::endl;
    for(size_t i=0;i<20;i++){
        if(bwt_rep.sigma_freq[i]!=0){
            std::cout<<"bytes: "<<i<<" "<<double(bwt_rep.sigma_freq[i])/double(tot_leaves)<<std::endl;
        }
    }

    std::cout<<"Children freq"<<std::endl;
    for(size_t i=0;i<20;i++){
        if(bwt_rep.n_children[i]!=0){
            std::cout<<"n_children: "<<i<<" "<<bwt_rep.n_children[i]<<std::endl;
        }
    }
    std::cout<<"Written bits in the tree "<<INT_CEIL(root.written_bits, 8)<<std::endl;
}

template<class bwt_dt_type>
void build_from_rl_plain(bwt_dt_type& bwt_rep, std::string& file){
}

template<class bwt_dt_type>
void build_from_plain(bwt_dt_type& bwt_rep, std::string& file){
}

template<class bwt_dt_type>
void build_dybl(bwt_dt_type& bwt, std::string& bwt_file, INPUT_FORMAT f){
    switch (f) {
        case INPUT_FORMAT::GRL_BWT:
            build_from_grl_bwt<bwt_dt_type>(bwt, bwt_file);
            break;
        case INPUT_FORMAT::RL_PLAIN:
            build_from_rl_plain<bwt_dt_type>(bwt, bwt_file);
            break;
        case INPUT_FORMAT::PLAIN:
            build_from_plain<bwt_dt_type>(bwt, bwt_file);
            break;
        default:
            std::cout<<"Unknown format"<<std::endl;
            exit(1);
    }
}


template<size_t b_size, size_t b_runs, size_t s_factor>
struct rlbwt_dybl{

    static constexpr size_t block_size = b_size;
    static constexpr size_t scale_factor = s_factor;
    static constexpr size_t max_block_runs = b_runs;

    size_t orig_runs=0;
    size_t eff_runs=0;
    size_t tot_syms=0;
    size_t sigma=0;

    const size_t levels=0;
    std::vector<size_t> stream;

    //TODO I will not use this information in the future
    uint64_t r_freq[b_runs+1]={0};//number of runs in a leaf
    uint64_t lvl_freq[20]={0};//the height of each leaf
    uint64_t sigma_freq[20]={0};//the alphabet of each leaf
    uint64_t n_children[100]={0};

    rlbwt_dybl():levels(size_t(ceil(log(b_size)/log(s_factor)) - ceil(log(b_runs)/log(s_factor)))+1){
        // logarithm function to calculate value
        float lg = log(b_size) / log(s_factor);
        float lg2 = log(b_runs) / log(s_factor);
        assert(lg==floor(lg));
        assert(lg2==floor(lg2));
    }

    //TODO static asserts in block_size, scale_factor, and b_runs
};

#endif //BWT_DTS_BENCHMARKS_DY_BL_RLBWT_H
