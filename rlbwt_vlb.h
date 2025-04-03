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
    size_t node_sigma=0;//number of symbols under the parent node
    size_t node_n_bits=0;//number of bits required for the subtree rooted under this node
    size_t consumed_syms=0;//number of symbols scanned

    bool lm_tree_branch=true;//is this node in the leftmost branch of its tree
    bool rm_tree_branch=true;//is this node in the rightmost branch of its tree

    const size_t lvl;//level of the subtree
    const size_t b_size;//block size for the level
    const size_t s_factor = bwt_dt_type::scale_factor;//shrinking factor for further subdivision
    const size_t b_runs = bwt_dt_type::max_block_runs;//maximum number of runs in a sequence of blocks

    rl_node *tmp_node = nullptr;

    bwt_dt_type& bwt_rep;//data structure encoding the representation
    std::vector<block_type> active_blocks;//run-length compressed blocks conforming a tree node
    std::vector<bool> node_sigma_bv;//buffer to compute the leaf's effective alphabet
    std::vector<std::vector<bool>> succ_pred_info;// bits indicating successor/predecessor info for each symbol in the alphabet
    std::vector<uint64_t> block_ranks;//rank information we store in the header of every internal node
    std::vector<uint8_t> packed_alphabet;//leaf's packed alphabet
    std::vector<uint64_t> block_ptr;//pointers to the node's children

    //TODO remove later
    std::vector<std::vector<uint64_t>> tree_sigma_dist;
    //

    explicit rl_node(size_t _lvl, size_t _b_size, bwt_dt_type& _bwt_rep):
                     lvl(_lvl),
                     b_size(_b_size),
                     bwt_rep(_bwt_rep),
                     active_blocks(b_runs),
                     node_sigma_bv(bwt_rep.sigma, false),
                     succ_pred_info(bwt_rep.sigma, std::vector<bool>(s_factor, false)),
                     packed_alphabet(bwt_rep.sigma, 0),
                     block_ranks(bwt_rep.sigma, 0){
        if(lvl==0){
            //number of trees in the forst
            block_ptr.resize(INT_CEIL(bwt_rep.tot_syms, b_size));
            node_sigma = bwt_rep.sigma;
            //TODO remove later
            tree_sigma_dist.resize(bwt_rep.sigma);
            //
        } else{
            block_ptr.resize(s_factor);
        }
    }

    inline void process_block_seq() {

        if(bk_id==1){//only one block in the sequence and it exceeds the limit of runs
            create_int_node(active_blocks[0]);//recursive partitioning
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
                create_int_node(active_blocks[0]);
                active_blocks[0].clear();
                bk_id=0;
                acc_runs=0;
            }
        }
    }

    inline void finish_run_scan(){

        //handle the last sequence of blocks
        assert(bk_len<=b_size);
        acc_runs+=active_blocks[bk_id].size();
        bk_id+=!active_blocks[bk_id].empty();

        if(acc_runs>=b_runs){//the last block sequence has more than the maximum number of allowed runs
            process_block_seq();
        }

        assert(acc_runs<b_runs);
        if(acc_runs>0){//the block sequence still has some runs left
            create_leaf(bk_id);
            for(size_t i=0;i<bk_id;i++){//clear the other blocks in the sequence
                active_blocks[i].clear();
            }
            bk_id=0;
            acc_runs=0;
        }
        bk_len=0;
        assert(aligned<8>(node_n_bits));
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
                if(acc_runs>=b_runs){
                    process_block_seq();
                }
                len -=split_run_len;
                bk_len=0;
            }
        }
    }

    inline void finish_int_node(size_t parent_sigma, std::vector<uint64_t>& parent_ranks){

        //TODO
        //1 bit to indicate it is a internal node
        //s bits to indicate which children were collapsed
        //parent_sigma bits to indicate the effective alphabet of the node with respect to the alphabet of its parent
        //the rank information for the node
        //n_children*node_sigma to indicate successor/predecessor sibling for each symbol
        //n_children pointers to the children

        //children information
        for(size_t i=0;i<n_children;i++){
            //TODO add the pointers
        }

        for(auto && s : node_sigma_bv){
            if(s){
                //TODO add the rank information for the symbols within block
            }
        }

        for(size_t i=0;i<n_children;i++){
            //TODO add the stream of each children
        }

        size_t header_bits=1;//to indicate it is an internal node

        //parent_sigma bits denote which symbols of the parent are in the node
        header_bits += parent_sigma;

        assert(lvl>0);
        if(lvl==1) {
            //TODO fix this
            //global rank information
            header_bits += node_sigma*sym_width(bwt_rep.tot_syms);
            //global predecessor/successor information
            header_bits += (parent_sigma-node_sigma) * sym_width(INT_CEIL(bwt_rep.tot_syms, b_size))*2;

            //TODO remove later (just testing)
            bwt_rep.tree_rank_header_overhead+=node_sigma*sym_width(bwt_rep.tot_syms);
            bwt_rep.tree_su_pr_header_overhead+=(parent_sigma-node_sigma) * sym_width(INT_CEIL(bwt_rep.tot_syms, b_size))*2;
            //
        } else {
            //bsize*s_factor is the block size of the parent
            //(node_sigma * sym_widths(b_size*s_factor)) for the ranks of the symbols under the node
            header_bits += node_sigma*sym_width(b_size*s_factor);
        }

        //these bits store the successor/predecessor information for each symbol in each child of the node
        header_bits+= node_sigma*n_children;

        //s_factor bits to indicate which children are collapsed
        header_bits += s_factor;

        //pt_bits indicates how many bits we use to encode pointers;
        size_t pt_bits = sym_width(node_n_bits/8);

        //pt_bits*n_children are the pointers
        header_bits += pt_bits + (pt_bits*n_children);
        //byte align *this internal node
        header_bits = INT_CEIL(header_bits, 8)*8;

        node_n_bits+=header_bits;

        bwt_rep.header_overhead+=header_bits;
    }

    inline void finish_forest() {

        //pointer information
        //TODO
        size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, b_size);//original number of blocks in the first level of the tree
        //

        //headers with the rank information
        //pt_bits indicates how many bits we use to encode pointers:
        //bits_node/8 is the pointer and b_runs indicate collapsed blocks
        size_t pt_bits = sym_width(node_n_bits/8) + sym_width(b_runs-1);
        size_t header_bits= pt_bits + (pt_bits*n_blocks);
        //byte align *this internal node
        header_bits = INT_CEIL(header_bits, 8)*8;

        std::vector<uint64_t> counts(10001);
        size_t tot=0;
        for(size_t s=0;s<bwt_rep.sigma;s++){
            std::cout<<"symbol "<<s<<" appearing in "<<tree_sigma_dist[s].size()<<" trees width diffs: ";
            for(size_t j=1;j<std::min<size_t>(tree_sigma_dist[s].size(), 100);j++){
                size_t diff = tree_sigma_dist[s][j]-tree_sigma_dist[s][j-1];
                if(j<100){
                    std::cout<<tree_sigma_dist[s][j]-tree_sigma_dist[s][j-1]<<" ";
                }
                counts[std::min<size_t>(diff, 10000)]++;
                tot++;
            }
            std::cout<<""<<std::endl;
        }

        double acc=0;
        for(size_t i=0;i<10000;i++){
            double fr = double(counts[i])/double(tot);
            acc+=fr;
            if(counts[i]!=0){
                std::cout<<"Dist "<<i<<" "<<fr<<" "<<acc<<std::endl;
            }
        }

        double fr = double(counts[10000])/double(tot);
        acc+=fr;
        std::cout<<"Dist +999 "<<fr<<" "<<acc<<std::endl;

        node_n_bits+=header_bits;
        bwt_rep.header_overhead+=header_bits;
    }


    inline void create_leaf_int(std::vector<block_type>& blocks, size_t n_blocks, size_t parent_sigma){

        assert(node_n_bits==0);
        size_t sym, len, n_runs=0, longest_run=0;
        for(size_t j=0;j<(blocks[0].size()-1);j++){
            sym = blocks[0][j].first;
            len = blocks[0][j].second;

            node_sigma_bv[sym]=true;
            assert(sym<bwt_rep.sigma);
            if(len>longest_run) longest_run = len;
            block_ranks[sym]+=len;
        }

        n_runs+=blocks[0].size()-1;
        for(size_t i=1;i<n_blocks;i++){

            //read the rightmost run of the previous block
            sym = blocks[i-1].back().first;
            len = blocks[i-1].back().second;

            //collapse the run with the first run of the current block if they have the same symbol
            if(sym==blocks[i][0].first){
                //collapse the runs as they are the same
                blocks[i][0].second +=len;
                blocks[i-1].pop_back();
            } else {
                //otherwise process the last run of the previous block as an independent run
                node_sigma_bv[sym] = true;
                assert(sym<bwt_rep.sigma);
                if(len>longest_run) longest_run = len;
                block_ranks[sym]+=len;
                n_runs++;
            }

            for(size_t j=0;j<(blocks[i].size()-1);j++){
                sym = blocks[i][j].first;
                len = blocks[i][j].second;
                node_sigma_bv[sym]=true;
                if(len>longest_run) longest_run = len;
                block_ranks[sym]+=len;
            }
            n_runs+=blocks[i].size()-1;
        }

        //process the last run
        sym = blocks[n_blocks-1].back().first;
        len = blocks[n_blocks-1].back().second;

        node_sigma_bv[sym]=true;
        assert(sym<bwt_rep.sigma);
        if(len>longest_run) longest_run = len;
        block_ranks[sym]+=len;
        n_runs++;

        //compute the block's alphabet size
        size_t tmp_s=0;
        node_sigma=0;
        for(auto && s : node_sigma_bv){
            packed_alphabet[tmp_s++]=node_sigma;
            node_sigma+=s;
        }

        for(size_t i=0;i<n_blocks;i++){
            for(auto & run : blocks[i]){
                //pack the run symbol
                run.first = packed_alphabet[run.first];
                //TODO store the run
            }
        }

        size_t header_bits=1;//to indicate this node is a leaf
        header_bits+=3;//to indicate the encoding of the run
        header_bits+=parent_sigma;//to indicate the leaf's effective alphabet

        if(lvl==1){
            //TODO fix
            //global rank information (i.e., previous trees)
            header_bits+= sym_width(bwt_rep.tot_syms)*node_sigma;
            //global predecessor/successor information
            header_bits+= (parent_sigma-node_sigma) * sym_width(INT_CEIL(bwt_rep.tot_syms, b_size))*2;

            //TODO remove later (just testing)
            bwt_rep.tree_rank_header_overhead+=sym_width(bwt_rep.tot_syms)*node_sigma;
            bwt_rep.tree_su_pr_header_overhead+=(parent_sigma-node_sigma) * sym_width(INT_CEIL(bwt_rep.tot_syms, b_size))*2;
            //
        } else {
            //local rank information (i.e., previous siblings)
            header_bits += sym_width(b_size*s_factor)*node_sigma;
        }

        //the runs are byte-aligned
        header_bits = INT_CEIL(header_bits, 8)*8;
        node_n_bits = header_bits;

        assert(n_runs<=bwt_dt_type::max_block_runs);
        size_t max_bytes_per_run = INT_CEIL((sym_width(node_sigma)+sym_width(longest_run)), 8);
        //TODO testing
        size_t leaf_enc=0;
        if(max_bytes_per_run==2){
            size_t bytes_per_run;
            size_t tmp[3]={0};
            for(size_t i=0;i<n_blocks;i++){
                for(auto & run : blocks[i]){
                    bytes_per_run = INT_CEIL((sym_width(run.first)+sym_width(run.second)), 8);
                    tmp[bytes_per_run]++;
                }
            }
            assert((tmp[1]+tmp[2])==n_runs);
            size_t vbyte_total = tmp[1] + tmp[2]*2 + INT_CEIL(n_runs, 8);
            //byte encoding for the runs of this leaf
            if(vbyte_total<(max_bytes_per_run*n_runs)){
                node_n_bits += vbyte_total*8;
                bwt_rep.runs_overhead+= vbyte_total*8;
            }else{
                node_n_bits += max_bytes_per_run*n_runs*8;
                bwt_rep.runs_overhead+= max_bytes_per_run*n_runs*8;
                leaf_enc=max_bytes_per_run;
            }
        } else {
            node_n_bits += max_bytes_per_run*n_runs*8;
            bwt_rep.runs_overhead+= max_bytes_per_run*n_runs*8;
            leaf_enc=max_bytes_per_run;
        }
        //

        //gather some statistics
        bwt_rep.header_overhead+=header_bits;
        bwt_rep.eff_runs += n_runs;
        bwt_rep.r_freq[n_runs]++;
        bwt_rep.lvl_freq[lvl-1]++;//lvl=0 is the forest, so it doesn't count. lvl=1 is a root of a tree
        bwt_rep.leaf_enc_freq[leaf_enc]++;
    }

    inline void create_leaf(size_t n_blocks) {

        assert(n_blocks>0);
        assert(aligned<8>(node_n_bits));

        bool lm_child = consumed_syms==0;
        bool rm_child = ((consumed_syms+(b_size*n_blocks))==(b_size*s_factor));

        //std::cout<<lm_child<<"/"<<rm_child<<" / "<<b_size<<" / "<<b_size*s_factor<<" / "<<n_children<<std::endl;

        tmp_node->lm_tree_branch = lm_tree_branch && lm_child;
        tmp_node->rm_tree_branch = rm_tree_branch && rm_child;

        tmp_node->create_leaf_int(active_blocks, n_blocks, node_sigma);

        //add the rank information of the active child node (next_node) to the
        // parent's rank information
        for(size_t s=0;s<bwt_rep.sigma;s++){
            block_ranks[s]+=tmp_node->block_ranks[s];
        }

        if(lm_tree_branch){
            //
        }

        if(rm_tree_branch){
            //
        }


        //add successor/predecessor information
        if(lvl>0){//lvl=0 is the forest, so it doesn't include this information
            //compute successor/predecessor info for the parent node
            size_t s_comp=0;
            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(node_sigma_bv[s]){
                    //TODO fix this because it is not correct
                    succ_pred_info[s_comp++][n_children]=tmp_node->block_ranks[s]>0;
                    assert((tmp_node->block_ranks[s]>0) == tmp_node->node_sigma_bv[s]);
                }
            }
            assert(s_comp==node_sigma);

            if(lvl==1){
                //TODO add successor/predecessor pointer to other trees
            }
        }else{
            //TODO remove later, just testing
            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(tmp_node->node_sigma_bv[s]){
                    tree_sigma_dist[s].push_back(n_children);
                }
            }
            //
        }

        //the leaf pointer should be byte-aligned
        //notice the pointer points to the end of the block.
        node_n_bits+=tmp_node->node_n_bits;
        block_ptr[n_children++] = node_n_bits/8;

        tmp_node->reset();
        consumed_syms+=b_size*n_blocks;
    }

    inline void get_node_alphabet(const block_type& block){
        //compute the alphabet of the block first.
        // We need it beforehand
        for(auto & run : block){
           node_sigma_bv[run.first] = true;
        }
        node_sigma = 0;
        for(auto const& bit : node_sigma_bv){
            node_sigma+=bit;
        }
    }

    inline void create_int_node(block_type& block) {

        assert(aligned<8>(node_n_bits));//check it is byte-aligned

        bool lm_child = consumed_syms == 0;
        bool rm_child = ((consumed_syms+b_size)==(b_size*s_factor));

        tmp_node->lm_tree_branch = lm_tree_branch && lm_child;
        tmp_node->rm_tree_branch = rm_tree_branch && rm_child;

        //std::cout<<lm_child<<"/"<<rm_child<<" -> "<<b_size<<" / "<<b_size*s_factor<<" / "<<n_children<<std::endl;

        tmp_node->get_node_alphabet(block);
        for(auto & run : block){
            tmp_node->process_run(run.first, run.second);
        }
        tmp_node->finish_run_scan();
        tmp_node->finish_int_node(node_sigma, block_ranks);

        //add the rank information of the active child node (next_node) to the
        // parent's rank information
        for(size_t s=0;s<bwt_rep.sigma;s++){
            block_ranks[s] += tmp_node->block_ranks[s];
        }

        if(lm_tree_branch){
            //
        }

        if(rm_tree_branch){
            //
        }

        //add successor/predecessor information (lvl=0 is the forest, so it does not count)
        if(lvl>0) {
            size_t s_comp=0;
            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(node_sigma_bv[s]){
                    succ_pred_info[s_comp++][n_children]=tmp_node->block_ranks[s]>0;
                    assert((tmp_node->block_ranks[s]>0) == tmp_node->node_sigma_bv[s]);
                }
            }
            assert(s_comp==node_sigma);
            if(lvl==1){
                //TODO add successor/predecessor pointer to other trees
            }
        } else {
            //TODO remove later, just testing
            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(tmp_node->node_sigma_bv[s]){
                    tree_sigma_dist[s].push_back(n_children);
                }
            }
            //
        }

        //the pointer to the active child node (next_node) should be aligned
        node_n_bits+=tmp_node->node_n_bits;
        block_ptr[n_children++]=(node_n_bits/8);

        bwt_rep.children_freq[tmp_node->n_children]++;

        tmp_node->reset();
        consumed_syms+=b_size;
    }

    inline void reset(){
        //TODO replace these loops with a memset
        for(unsigned long long & rank : block_ranks) rank = 0;
        for(auto && s : node_sigma_bv) s=false;
        for(size_t s=0;s<bwt_rep.sigma;s++){
            for(size_t c=0;c<s_factor;c++){
                succ_pred_info[s][c]= false;
            }
        }
        //
        n_children = 0;
        node_n_bits = 0;
        node_sigma = 0;
        consumed_syms = 0;
        lm_tree_branch = lvl==1;//lvl=1 means the root of the tree
        lm_tree_branch = lvl==1;
    }
};

template<class bwt_dt_type>
void forest_stats(bwt_dt_type& bwt_rep,  std::vector<rl_node<bwt_dt_type>>& tmp_nodes){

    for(size_t s=0;s<bwt_rep.sigma;s++){
        std::cout<<"\tsymbol "<<s<<", rank: "<<tmp_nodes[0].block_ranks[s]<<std::endl;
    }

    std::cout<<"Number_of_runs_in_a_leaf dist: "<<std::endl;
    assert(bwt_rep.r_freq[0]==0);
    for(size_t r=1;r<=bwt_dt_type::max_block_runs;r++){
        std::cout<<"\t"<<r<<" : "<<bwt_rep.r_freq[r]<<std::endl;
    }
    std::cout<<"Total number of runs versus original number of runs: "<<bwt_rep.eff_runs<<" / "<<bwt_rep.orig_runs<<std::endl;
    std::cout<<"Increase in the number of runs "<<(double(bwt_rep.eff_runs)/double(bwt_rep.orig_runs)-1)*100<<"%"<<std::endl;

    std::cout<<"Tree_depth dist: "<<bwt_rep.levels<<std::endl;
    size_t tot_leaves=0;
    for(size_t i=0;i<20;i++){
        tot_leaves+=bwt_rep.lvl_freq[i];
    }
    for(size_t i=0;i<20;i++){
        if(bwt_rep.lvl_freq[i]>0){
            std::cout<<"\t"<<i<<": "<<double(bwt_rep.lvl_freq[i])/double(tot_leaves)<<std::endl;
        }
    }

    std::cout<<"Leaf_encoding dist: "<<std::endl;
    for(size_t i=0;i<20;i++){
        if(bwt_rep.leaf_enc_freq[i]!=0){
            if(i>0){
                std::cout<<"\tFixed "<<i<<" bytes:\t\t\t\t"<<double(bwt_rep.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
            }else{
                std::cout<<"\tVariable-length 1/2 bytes: "<<double(bwt_rep.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
            }
        }
    }

    std::cout<<"Number_of_children dist: "<<std::endl;
    size_t del_nodes=0, tot_nodes=0;
    for(size_t i=0;i<20;i++){
        if(bwt_rep.children_freq[i]!=0){
            std::cout<<"\t"<<i<<": "<<bwt_rep.children_freq[i]<<std::endl;
            del_nodes+=(bwt_rep.scale_factor-i)*bwt_rep.children_freq[i];
            tot_nodes+=bwt_rep.children_freq[i];
        }
    }

    size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, bwt_rep.block_size);//original number of blocks in the first level of the tree
    std::cout<<"Effective number of trees versus full number of trees (n/b): "<<tmp_nodes[0].n_children<<" / "<<n_blocks<<std::endl;
    std::cout<<"Percentage of removed trees: "<<(1-double(tmp_nodes[0].n_children)/double(n_blocks))*100<<"% "<<std::endl;

    tot_nodes*=bwt_rep.scale_factor;
    tot_nodes+=n_blocks;
    del_nodes=n_blocks-tmp_nodes[0].n_children;

    std::cout<<"Percentage of removed nodes: "<<(double(del_nodes)/double(tot_nodes))*100<<"% "<<std::endl;
    std::cout<<"Written bytes in the data structure: "<<INT_CEIL(tmp_nodes[0].node_n_bits, 8)<<std::endl;
    std::cout<<"Written bytes without the tree's succ/pred info: "<<INT_CEIL((tmp_nodes[0].node_n_bits-bwt_rep.tree_su_pr_header_overhead), 8)<<std::endl;
    std::cout<<"Headers' contribution to the final space: "<<(double(bwt_rep.header_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"% "<<std::endl;
    std::cout<<"Runs' contribution to the final space: "<<(double(bwt_rep.runs_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"% "<<std::endl;
    std::cout<<"Tree rank headers' contribution to the final space: "<<(double(bwt_rep.tree_rank_header_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"% "<<std::endl;
    std::cout<<"Tree succ/pred headers' contribution to the final space: "<<(double(bwt_rep.tree_su_pr_header_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"% "<<std::endl;
    assert(bwt_rep.header_overhead+bwt_rep.runs_overhead==tmp_nodes[0].node_n_bits);
    std::cout<<"space_usage:"<<float(tmp_nodes[0].node_n_bits)/float(bwt_rep.tot_syms)<<" bps"<<std::endl;\

}

template<class bwt_dt_type>
void build_from_grlbwt(bwt_dt_type& bwt_rep, std::string& bwt_file){

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
    size_t b_size = bwt_dt_type::block_size;
    for(size_t i=0;i<=bwt_rep.levels;i++){
        tmp_nodes.push_back(rl_node(i, b_size, bwt_rep));
        b_size/=bwt_dt_type::scale_factor;
    }
    for(size_t i=0;i<bwt_rep.levels;i++){
        tmp_nodes[i].tmp_node = &tmp_nodes[i+1];
    }

    //compute the forest
    for(size_t i=0;i<n_runs;i++){
        bwt_buff.read_run(i, sym, len);
        tmp_nodes[0].process_run(sym_map[sym], len);
    }
    tmp_nodes[0].finish_run_scan();
    tmp_nodes[0].finish_forest();
    //

    forest_stats(bwt_rep, tmp_nodes);
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
    uint64_t leaf_enc_freq[20]={0};//encoding of each leaf
    uint64_t children_freq[100]={0};//children frequency = how many nodes with 1,2,...,x children
    uint64_t header_overhead=0;//number of bits used by the headers of the nodes
    uint64_t runs_overhead=0;
    uint64_t tree_rank_header_overhead=0;
    uint64_t tree_su_pr_header_overhead=0;

    rlbwt_vlb():levels(size_t(ceil(log(b_size)/log(s_factor)) - ceil(log(b_runs)/log(s_factor)))+1){
        // logarithm function to calculate value
        float lg = log(b_size) / log(s_factor);
        assert(lg==floor(lg));
        float lg2 = log(b_runs) / log(s_factor);
        assert(lg2==floor(lg2));
    }

    //TODO static asserts in block_size, scale_factor, and b_runs
};

#endif //RLBWT_DYBL_H