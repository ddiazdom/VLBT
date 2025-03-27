//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef RLBWT_DYBL_H
#define RLBWT_DYBL_H

#include <cmath>

using block_type = std::vector<std::pair<uint32_t, size_t>>;

enum INPUT_FORMAT{
   GRL_BWT=0,
   RL_PLAIN=1,
   PLAIN=2
};

template<class bwt_dt_type>
struct rl_node {//state of the compression

    size_t bk_len=0;//length of the active block
    size_t bk_id=0;//id of the active block
    size_t acc_runs=0;//sum of the runs in the active blocks
    size_t n_children=0;//number of children of this node
    size_t node_n_bits=0;//number of bits required for the subtree

    const size_t lvl;//level of the subtree
    const size_t b_size;
    const size_t s_factor = bwt_dt_type::scale_factor;//shrinking factor for further subdivision
    const size_t b_runs = bwt_dt_type::max_block_runs;//maximum number of runs in a sequence of blocks

    rl_node *next_node = nullptr;

    bwt_dt_type& bwt_rep;//data structure encoding the representation
    std::vector<block_type> active_blocks;//run-length compressed blocks conforming a tree node
    std::vector<bool> sigma_buff;//buffer to compute the leaf's effective alphabet
    std::vector<std::vector<uint64_t>> block_ranks;//rank information we store in the header of every internal node
    std::vector<uint8_t> packed_alphabet;//leaf's packed alphabet
    std::vector<uint64_t> block_ptr;//pointers to the node's children

    explicit rl_node(size_t _lvl, size_t _b_size, bwt_dt_type& _bwt_rep):
                     lvl(_lvl),
                     b_size(_b_size),
                     bwt_rep(_bwt_rep),
                     active_blocks(b_runs),
                     sigma_buff(bwt_rep.sigma, false),
                     packed_alphabet(bwt_rep.sigma, 0) {
        if(lvl==0){
            size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, b_size);
            block_ptr.resize(n_blocks);
            block_ranks = std::vector<std::vector<uint64_t>>(n_blocks+1, std::vector<uint64_t>(bwt_rep.sigma, 0));
        }else{
            block_ptr.resize(s_factor);
            block_ranks = std::vector<std::vector<uint64_t>>(s_factor+1, std::vector<uint64_t>(bwt_rep.sigma, 0));
        }
    }

    inline void process_run(const size_t& sym, size_t& len) {

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

                acc_runs+=active_blocks[bk_id].size();
                bk_id++;

                //we compete a new block sequence
                if(acc_runs>b_runs){

                    if(bk_id==1){//only one block in the sequence and it exceeds the limit of runs
                        create_subtree(active_blocks[0]);//recursive partitioning
                        active_blocks[0].clear();
                        bk_id=0;
                        acc_runs=0;
                    } else {//multiple blocks in the block sequence
                        create_leaf(bk_id-1);//we do not consider the last block in the sequence
                        active_blocks[0].swap(active_blocks[bk_id-1]);//move the last block to the beginning of the sequence
                        for(size_t i=1;i<bk_id;i++){//clear the other blocks in the sequence
                            active_blocks[i].clear();
                        }

                        //new block sequence
                        acc_runs = active_blocks[0].size();
                        bk_id=1;

                        if(active_blocks[0].size()>b_runs){//process the remaining block if it exceeds the limit of runs
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

    inline void finish_int_node(){

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

        //children information
        for(size_t i=0;i<n_children;i++){
            //TODO add the pointers
        }

        //get the effective alphabet of the block
        std::vector<bool> eff_b_alph(bwt_rep.sigma, false);
        std::vector<std::vector<bool>> child_alpha(n_children, std::vector<bool>(bwt_rep.sigma, false));
        bool symbol_exists;
        for(size_t c=0;c<n_children;c++){
            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(c==0){
                    symbol_exists = (block_ranks[c][s]>0);
                }else{
                    symbol_exists = (block_ranks[c-1][s]!=block_ranks[c][s]);
                }
                eff_b_alph[s]=eff_b_alph[s] | symbol_exists;
                child_alpha[c][s]=symbol_exists;
            }
        }
        //
        size_t p_sigma=0;
        for(auto const& bit : eff_b_alph){
            p_sigma+=bit;
        }

        //get the effective alphabet of the child
        size_t alt_header=p_sigma+p_sigma*n_children;
        for(size_t c=0;c<n_children;c++) {
            //TODO add the rank information
            //size_t n_empty=0;
            size_t c_sigma=0;
            //std::cout<<c<<" -> ";
            for(size_t s=0;s<bwt_rep.sigma;s++){
                /*if(c==0){
                    n_empty+=block_ranks[c][s]==0;
                }else{
                    n_empty+=block_ranks[c-1][s]==block_ranks[c][s];
                }*/
                //std::cout<<""<<s<<"="<<block_ranks[c][s]<<" ";
                c_sigma+=child_alpha[c][s];
            }
            alt_header+=c_sigma*sym_width(b_size*s_factor);
            //std::cout<<n_empty<<", "<<e<<", "<<l<<" | ";
        }
        //std::cout<<""<<std::endl;

        for(size_t i=0;i<n_children;i++){
            //TODO add the stream of each children
        }

        //bwt_ep.sigma*sym_width(b_size*s_factor) bits to store the rank info of every child
        //(n_children*sigma_bits) bits for the children's rank information
        //size_t sigma_bits = bwt_rep.sigma*sym_width(b_size*s_factor);
        //std::cout<<sigma_bits<<" "<<alt_header<<std::endl;
        //size_t header_bits = (n_children*sigma_bits);
        size_t header_bits = alt_header;

        //s_factor bits to indicate which children are collapsed
        header_bits += s_factor;

        //pt_bits indicates how many bits we use to encode pointers;
        size_t pt_bits = sym_width(node_n_bits/8);

        //pt_bits*n_children are the pointers
        header_bits += pt_bits + (pt_bits*n_children);
        node_n_bits+=header_bits;

        //byte align *this internal node
        node_n_bits = INT_CEIL(node_n_bits, 8)*8;
    }

    inline void finish_root() {

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
                // because it is the one exceeding the number of runs
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

        assert(node_n_bits==(INT_CEIL(node_n_bits, 8)*8));

        //pointer information
        size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, b_size);//original number of blocks in the first level of the tree
        std::cout<<"Total n_blocks "<<n_blocks<<" versus "<<n_children<<std::endl;

        //headers with the rank information
        size_t header_bits = n_children*bwt_rep.sigma_bits;
        //pt_bits indicates how many bits we use to encode pointers:
        //bits_node/8 is the pointer and b_runs indicate collapsed blocks
        size_t pt_bits = sym_width(node_n_bits/8) + sym_width(b_runs-1);
        header_bits+= pt_bits + (pt_bits*n_blocks);
        node_n_bits+=header_bits;

        //byte align *this internal node
        node_n_bits = INT_CEIL(node_n_bits, 8)*8;
    }

    inline void create_leaf(size_t n_blocks) {

        assert(n_blocks>0);
        assert(node_n_bits==(INT_CEIL(node_n_bits, 8)*8));

        size_t sym, len, n_runs=0, longest_run=0;
        for(size_t j=0;j<(active_blocks[0].size()-1);j++){
            sym = active_blocks[0][j].first;
            len = active_blocks[0][j].second;

            sigma_buff[sym]=true;
            assert(sym<bwt_rep.sigma);
            if(len>longest_run) longest_run = len;
            block_ranks[n_children][sym]+=len;
        }
        n_runs+=active_blocks[0].size()-1;

        for(size_t i=1;i<n_blocks;i++){

            //read the rightmost run of the previous block
            sym = active_blocks[i-1].back().first;
            len = active_blocks[i-1].back().second;

            //collapse the run with the first run of the current block if they have the same symbol
            if(sym==active_blocks[i][0].first){
                //collapse the runs as they are the same
                active_blocks[i][0].second +=len;
                active_blocks[i-1].pop_back();
            } else {
                //otherwise process the last run of the previous block as an independent run
                sigma_buff[sym] = true;
                assert(sym<bwt_rep.sigma);
                if(len>longest_run) longest_run = len;
                block_ranks[n_children][sym]+=len;
                n_runs++;
            }

            for(size_t j=0;j<(active_blocks[i].size()-1);j++){
                sym = active_blocks[i][j].first;
                len = active_blocks[i][j].second;
                sigma_buff[sym]=true;
                if(len>longest_run) longest_run = len;
                block_ranks[n_children][sym]+=len;
            }
            n_runs+=active_blocks[i].size()-1;
        }

        //process the last run
        sym = active_blocks[n_blocks-1].back().first;
        len = active_blocks[n_blocks-1].back().second;

        sigma_buff[sym]=true;
        assert(sym<bwt_rep.sigma);
        if(len>longest_run) longest_run = len;
        block_ranks[n_children][sym]+=len;
        n_runs++;

        //compute the block's alphabet size
        size_t b_sigma=0, tmp_s=0;
        for(auto && bit : sigma_buff){
            packed_alphabet[tmp_s++]=b_sigma;
            b_sigma+=bit;
            bit=false;
        }

        for(size_t i=0;i<n_blocks;i++){
            for(auto & run : active_blocks[i]){
                //pack the run symbol
                run.first = packed_alphabet[run.first];
                //TODO store the run
            }
        }

        //add the ranks within the block to the parent node
        for(size_t s=0;s<bwt_rep.sigma;s++){
            block_ranks[n_children+1][s]=block_ranks[n_children][s];
        }

        size_t bytes_per_run = INT_CEIL((sym_width(b_sigma)+sym_width(longest_run)), 8);
        bwt_rep.eff_runs += n_runs;
        assert(n_runs<=bwt_dt_type::max_block_runs);

        //TODO testing
        size_t vbyte_total=0;
        if(bytes_per_run==2){
            size_t tmp[4]={0};
            for(size_t i=0;i<n_blocks;i++){
                for(auto & run : active_blocks[i]){
                    bytes_per_run = INT_CEIL((sym_width(run.first)+sym_width(run.second)), 8);
                    tmp[bytes_per_run]++;
                }
            }
            vbyte_total = tmp[1] + tmp[2]*2 + INT_CEIL(n_runs, 8);
        }
        //

        bwt_rep.r_freq[n_runs]++;
        bwt_rep.lvl_freq[lvl]++;
        bwt_rep.sigma_freq[bytes_per_run]++;

        // leaf header
        // sigma bits encode the block's effective alphabet,
        // 3 bits indicate the number of bytes that the block uses to encode its runs
        // 1 bit indicates that this is a leaf
        node_n_bits += bwt_rep.sigma+4;
        //the runs are byte-aligned
        node_n_bits = INT_CEIL(node_n_bits, 8)*8;

        //bits encoding the actual runs
        if(vbyte_total<(bytes_per_run*n_runs)){
            node_n_bits += vbyte_total*8;
        }else{
            node_n_bits += bytes_per_run*n_runs*8;
        }

        //the leaf pointer should be byte-aligned
        //notice the pointer points to the end of the block.
        block_ptr[n_children++] = node_n_bits/8;
    }

    inline void create_subtree(block_type& block) {

        assert(node_n_bits==(INT_CEIL(node_n_bits, 8)*8));//check it is byte-aligned
        for(auto & run : block){
            next_node->process_run(run.first, run.second);
        }
        next_node->finish_int_node();

        //add the rank information of the active child node (next_node) to the parent's rank information block_ranks[n_children]
        for(size_t s=0;s<bwt_rep.sigma;s++){
            block_ranks[n_children][s] += next_node->block_ranks[next_node->n_children][s];
            block_ranks[n_children+1][s] = block_ranks[n_children][s];
        }

        //the pointer to the active child node (next_node) should be aligned
        node_n_bits+=next_node->node_n_bits;
        block_ptr[n_children++]=(node_n_bits/8);
        //effective alphabet of the area under this node
        for(size_t i=0;i<next_node->sigma_buff.size();i++){
            sigma_buff[i] = next_node->sigma_buff[i];
        }

        bwt_rep.children_freq[next_node->n_children]++;
        next_node->reset();
    }

    inline void reset(){
        for(size_t i=0;i<bk_id;i++){
            active_blocks[i].clear();
        }

        for(unsigned long long & rank : block_ranks[0]){
            rank = 0;
        }

        bk_len = 0;
        bk_id = 0;
        acc_runs = 0;
        n_children = 0;
        node_n_bits = 0;
    }
};

template<class bwt_dt_type>
void build_from_grlbwt(bwt_dt_type& bwt_rep, std::string& bwt_file){

    size_t b_runs = bwt_dt_type::max_block_runs;
    size_t b_size = bwt_dt_type::block_size;

    bwt_buff_reader bwt_buff(bwt_file);
    size_t n_runs = bwt_buff.size();
    bwt_rep.orig_runs = n_runs;

    //compute symbol frequencies
    size_t sym, len;
    std::vector<uint64_t> sym_freqs(256, 0);
    std::vector<uint8_t> sym_map(256, 0);
    for(size_t i=0;i<n_runs;i++){
        bwt_buff.read_run(i, sym, len);
        sym_freqs[sym]+=len;
    }
    //

    //compute the effective alphabet and the number of bits we require to encode it
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
    bwt_rep.sigma_bits = sigma_bits;
    //

    //compute the array C[1..\sigma]
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
    //

    std::vector<rl_node<bwt_dt_type>> tmp_nodes;
    tmp_nodes.reserve(bwt_rep.levels+1);
    for(size_t i=0;i<=bwt_rep.levels;i++){
        tmp_nodes.push_back(rl_node(i, b_size, bwt_rep));
        b_size/=bwt_dt_type::scale_factor;
    }
    for(size_t i=0;i<bwt_rep.levels;i++){
        tmp_nodes[i].next_node = &tmp_nodes[i+1];
    }
    tmp_nodes[bwt_rep.levels].packed_alphabet.resize(bwt_rep.sigma);

    for(size_t i=0;i<n_runs;i++){
        bwt_buff.read_run(i, sym, len);
        tmp_nodes[0].process_run(sym_map[sym], len);
    }
    tmp_nodes[0].finish_root();

    std::cout<<bwt_rep.orig_runs<<" "<<bwt_rep.eff_runs<<" -> "<<float(bwt_rep.orig_runs)/float(bwt_rep.eff_runs)<<std::endl;
    for(size_t s=0;s<sigma;s++){
        std::cout<<s<<" "<<tmp_nodes[0].block_ranks[tmp_nodes[0].n_children][s]<<std::endl;
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
        if(bwt_rep.children_freq[i]!=0){
            std::cout<<"n_children: "<<i<<" "<<bwt_rep.children_freq[i]<<std::endl;
        }
    }
    std::cout<<"Written bytes in the data structure "<<INT_CEIL(tmp_nodes[0].node_n_bits, 8)<<std::endl;
}

template<class bwt_dt_type>
void build_from_rl_plain(bwt_dt_type& bwt_rep, std::string& file){

}

template<class bwt_dt_type>
void build_from_plain(bwt_dt_type& bwt_rep, std::string& file){

}

template<class bwt_dt_type>
void build_rlbwt_vlb(bwt_dt_type& bwt, std::string& bwt_file, INPUT_FORMAT f){

    switch (f) {
        case INPUT_FORMAT::GRL_BWT:
            build_from_grlbwt<bwt_dt_type>(bwt, bwt_file);
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
struct rlbwt_vlb {

    static constexpr size_t block_size = b_size;
    static constexpr size_t scale_factor = s_factor;
    static constexpr size_t max_block_runs = b_runs;

    size_t orig_runs=0;
    size_t eff_runs=0;
    size_t tot_syms=0;
    size_t sigma=0;
    size_t sigma_bits=0;

    const size_t levels=0;
    std::vector<size_t> stream;

    //TODO I will not use this information in the future
    uint64_t r_freq[b_runs+1]={0};//number of runs in a leaf
    uint64_t lvl_freq[20]={0};//the height of each leaf
    uint64_t sigma_freq[20]={0};//the alphabet of each leaf
    uint64_t children_freq[100]={0};//children frequency = how many nodes with 1,2,...,x children

    rlbwt_vlb():levels(size_t(ceil(log(b_size)/log(s_factor)) - ceil(log(b_runs)/log(s_factor)))+1){
        // logarithm function to calculate value
        float lg = log(b_size) / log(s_factor);
        float lg2 = log(b_runs) / log(s_factor);
        assert(lg==floor(lg));
        assert(lg2==floor(lg2));
    }

    //TODO static asserts in block_size, scale_factor, and b_runs
};

#endif //RLBWT_DYBL_H