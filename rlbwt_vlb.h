//
// Created by Diaz, Diego on 20.11.2024.
//

#ifndef RLBWT_DYBL_H
#define RLBWT_DYBL_H

#include <cmath>
#ifdef __linux__
#include <malloc.h>
#endif

using run_type = std::pair<uint32_t, size_t>;
using block_type = std::vector<run_type>;

enum INPUT_FORMAT{
   GRL_BWT=0,
   RL_PLAIN=1,
   PLAIN=2
};

enum node_type {
    INTERNAL,
    LEAF
};

typedef bitstream<size_t> stream_type;

//statistics about the data structure
template<class bwt_dt_type>
struct stat_collector{
    uint64_t rpl_freq[bwt_dt_type::max_block_runs+1]={0};//number of runs in a leaf
    uint64_t leaf_depth_freq[20]={0};//the depth of each leaf
    uint64_t leaf_enc_freq[20]={0};//encoding of each leaf
    uint64_t children_freq[100]={0};//children frequency = how many nodes with 1,2,...,x children
    uint64_t header_overhead=0;//number of bits used by the headers of the nodes
    uint64_t runs_overhead=0;
    uint64_t rank_overhead=0;
    uint64_t ext_su_pr_overhead=0;
    uint64_t int_su_pr_overhead=0;
    uint64_t tree_pointers_overhead=0;
    uint64_t trees_overhead=0;
    uint64_t eff_runs=0;
    uint64_t max_n_blocks=0;
    uint64_t ext_pred_freq[257]={0};
    uint64_t ext_succ_freq[257]={0};
};

template<class bwt_dt_type>
struct rl_node {//state of the compression

    size_t bk_len=0;//length of the active block
    size_t bk_id=0;//id of the active block
    size_t acc_runs=0;//sum of the runs in the active blocks
    size_t n_children=0;//number of children of this node
    size_t node_sigma=0;//number of symbols under the parent node
    size_t node_n_bits=0;//number of bits required for the subtree rooted under this node
    size_t consumed_syms=0;//number of symbols under this node that have been scanned so far
    size_t child_rank=0;//this node is the child_rank of its parent
    size_t cov_symbols=0;//how many symbols of the input BWT does this node cover

    bool lm_tree_branch=true;//true if this node in the leftmost branch of its tree
    bool rm_tree_branch=true;//true if this node in the rightmost branch of its tree
    bool lm_child=false;//true if this node is the leftmost child of its parent
    bool rm_child=false;//true if this node is the rightmost child of its parent
    bool leaf=false;//true if the node is a leaf
    std::ostream& ofs;

    const size_t lvl;//level of the subtree
    const size_t b_size;//block size for the level
    const size_t s_factor = bwt_dt_type::scale_factor;//shrinking factor for further subdivision
    const size_t b_runs = bwt_dt_type::max_block_runs;//maximum number of runs in a sequence of blocks

    rl_node *tmp_node = nullptr;

    stream_type buffer;
    stream_type children_buffer;

    bwt_dt_type& bwt_rep; //data structure encoding the representation
    std::vector<block_type> active_blocks; //run-length compressed blocks conforming a tree node
    std::vector<bool> node_sigma_bv; //buffer to compute the leaf's effective alphabet
    std::vector<std::vector<bool>> int_succ_pred_info; //bits indicating successor/predecessor info for each symbol in the alphabet

    //next_ext_pred[s]=false (resp., next_ext_succ[s]=false) means the symbol s needs predecessor (resp, successor) information
    //these vectors only consider symbols that are *in* the alphabet of the node
    std::vector<bool> need_ext_pred;//symbols in the alphabet of the need that need external predecessor info (symbols not in the left branch)
    std::vector<bool> need_ext_succ;//symbols in the alphabet of the need that need external successor info (symbols not in the right branch)

    std::vector<uint64_t> block_ranks;//rank information we store in the header of every internal node
    std::vector<uint64_t> block_ptr;//pointers (byte offsets) to the node's children
    std::vector<uint64_t> tree_offset;//number of symbols in the text before each tree
    std::vector<uint8_t> packed_alphabet;//leaf's packed alphabet

    // the vector sigma_trees[s], with s \in \Sigma, is a strictly increasing
    // sequence encoding the trees in the forest containing the symbol s
    std::vector<std::vector<uint64_t>> sigma_trees;

    //list of symbols of each tree (as a bitvector) that require external predecessor information
    //That is, each symbol need_ext_pred[s] \cup the symbols not appearing in the tree
    std::vector<bool> ext_pred_info;

    //list of the symbols of each tree (as a bitvector) that require external successor information
    //That is, each symbol need_ext_succ[s] \cup the symbols not appearing in the tree
    //ext_pred_info and ext_succ_info are information for the nodes with level 1 (i.e., tree roots)
    std::vector<bool> ext_succ_info;

    //a struct to collect statistics about the data structure
    stat_collector<bwt_dt_type>& stats;
    std::vector<uint64_t>& C;//the standard C[1..\sigma] array of the FM-index

    explicit rl_node(size_t _lvl, size_t _b_size, bwt_dt_type& _bwt_rep, stat_collector<bwt_dt_type>& st,
                     std::vector<uint64_t>& C_, std::ofstream& _ofs):
                     lvl(_lvl),
                     b_size(_b_size),
                     bwt_rep(_bwt_rep),
                     active_blocks(b_runs),
                     node_sigma_bv(bwt_rep.sigma, false),
                     int_succ_pred_info(bwt_rep.sigma, std::vector<bool>(s_factor, false)),
                     packed_alphabet(bwt_rep.sigma, 0),
                     block_ranks(bwt_rep.sigma, 0),
                     need_ext_pred(bwt_rep.sigma, true),
                     need_ext_succ(bwt_rep.sigma, true),
                     stats(st),
                     C(C_),
                     ofs(_ofs){
        if(lvl==0){
            //number of trees in the forest
            size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, b_size);
            block_ptr.resize(n_blocks+1);
            block_ptr[0] = 0;
            node_sigma = bwt_rep.sigma;
            node_sigma_bv = std::vector<bool>(bwt_rep.sigma, true);

            sigma_trees.resize(bwt_rep.sigma);
            for(size_t s=0;s<bwt_rep.sigma;s++){
                sigma_trees[s].reserve(n_blocks);
            }
            ext_pred_info = std::vector<bool>(n_blocks*bwt_rep.sigma, false);
            ext_succ_info = std::vector<bool>(n_blocks*bwt_rep.sigma, false);
            tree_offset = std::vector<uint64_t>(n_blocks+1, 0);
        } else{
            block_ptr.resize(s_factor+1);
        }
    }

    inline void process_block_seq() {
        if(bk_id==1){//only one block in the sequence and it exceeds the limit of runs
            create_node<INTERNAL>(1);//recursive partitioning
            active_blocks[0].clear();
            bk_id=0;
            acc_runs=0;
        } else {//multiple blocks in the block sequence

            create_node<LEAF>(bk_id-1);//we do not consider the last block in the sequence
            active_blocks[0].swap(active_blocks[bk_id-1]);//move the last block to the beginning of the sequence
            for(size_t i=1;i<bk_id;i++){//clear the other blocks in the sequence
                active_blocks[i].clear();
            }

            //new block sequence
            acc_runs = active_blocks[0].size();
            bk_id=1;
            if(active_blocks[0].size()>b_runs){//process the remaining block if it exceeds the limit of runs
                create_node<INTERNAL>(1);
                active_blocks[0].clear();
                bk_id=0;
                acc_runs=0;
            }
        }
    }

    inline void destroy(){
        buffer.destroy();
        children_buffer.destroy();
        for(auto bk : active_blocks){
            destroy_vector(bk);
        }
        destroy_vector(node_sigma_bv);
        for(auto vec : int_succ_pred_info){
            destroy_vector(vec);
        }
        destroy_vector(need_ext_pred);
        destroy_vector(need_ext_succ);
        destroy_vector(block_ranks);
        destroy_vector(block_ptr);
        destroy_vector(tree_offset);
        destroy_vector(packed_alphabet);
#ifdef __linux__
        malloc_trim(0);
#endif
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
            create_node<LEAF>(bk_id);
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

    inline void finish_int_node(size_t parent_sigma,
                                const std::vector<bool>& parent_sigma_bv,
                                const std::vector<uint64_t>& parent_rank_info){

        assert(n_children<=s_factor);
        assert(lvl>0);

        //HEADER DESCRIPTION:
        //1 bit to indicate it is an internal node
        //parent_sigma bits to indicate the effective alphabet of the node with respect to the alphabet of its parent
        //r_width*node_sigma bits store the rank information
        //s bits to indicate which children were collapsed
        //n_children*node_sigma to indicate successor/predecessor sibling for each symbol
        //n_children pointers to the children

        //===some preliminary information
        //width in bits for the rank information
        size_t r_width;
        if(lvl==1) {
            //the root of the tree
            r_width = sym_width(bwt_rep.max_freq);
        } else {
            //internal node that is not the root
            //bsize*s_factor is the block size of the parent
            r_width = sym_width(b_size*s_factor);
        }
        size_t rank_bits = r_width*node_sigma;
        //amount of bits for the succ/pred info
        size_t su_pr_bv_bits = node_sigma*n_children;

        //number of bits we use to encode pointers to the children
        size_t pt_width=20;
        size_t pt_bits = sym_width(node_n_bits/8);
        assert(pt_bits<pt_width);

        //amount of information (in bits) for the header
        size_t header_bits = 1+parent_sigma+rank_bits+s_factor+su_pr_bv_bits+pt_width+(pt_bits*n_children);
        header_bits = INT_CEIL(header_bits, 8)*8;//byte-aligned

        buffer.reserve_in_bits(header_bits+node_n_bits);

        //THE ENCODING STARTS HERE
        //==1 bit (false) to indicate this is an internal node
        size_t bit_pos = 0;
        buffer.write(bit_pos, bit_pos, 0);
        bit_pos++;
        //==

        //==parent_sigma bits encode the internal node's effective alphabet
        for(size_t i=0;i<bwt_rep.sigma;i++){
            if(parent_sigma_bv[i]){
                buffer.write(bit_pos, bit_pos, node_sigma_bv[i]);
                bit_pos++;
            }
        }
        assert(bit_pos==(1+parent_sigma));
        //==

        //==rank information
        //rank_bits=r_width*node_sigma bits store the rank information
        for(size_t i=0;i<bwt_rep.sigma;i++){
            if(node_sigma_bv[i]){
                buffer.write(bit_pos, bit_pos+r_width-1, parent_rank_info[i]);
                bit_pos+=r_width;
            }
        }
        assert(bit_pos==(1+parent_sigma+rank_bits));
        //==

        //==
        //s_factor bits to indicate which children were collapsed
        //these bits give us the real number of children
        //TODO compute this information
        for(size_t i=0;i<s_factor;i++){
            buffer.write(bit_pos, bit_pos, true);
            bit_pos++;
        }
        assert(bit_pos==(1+parent_sigma+rank_bits+s_factor));
        //==

        //==int succ/prec info
        //these bits store the successor/predecessor information for each symbol in each child of the node
        for(size_t s_comp=0;s_comp<node_sigma;s_comp++){
            for(size_t c=0;c<n_children;c++){
                buffer.write(bit_pos, bit_pos, int_succ_pred_info[s_comp][c]);
                bit_pos++;
            }
        }
        assert(bit_pos==(1+parent_sigma+rank_bits+s_factor+su_pr_bv_bits));
        //==

        //==pt_width*n_children bits store pointers (byte offsets) to the children of this node
        buffer.write(bit_pos, bit_pos+pt_width-1, pt_bits);
        bit_pos+=pt_width;
        for(size_t i=0;i<n_children;i++){
            //NOTE consider the humber of bytes in the header when computing the byte position of a child.
            //i.e., child_byte_pos = h_bits/8 + block_ptr[c], where c is the child we need to find
            buffer.write(bit_pos, bit_pos+pt_bits-1, block_ptr[i]);
            bit_pos+=pt_bits;
        }
        assert(bit_pos==(1+parent_sigma+rank_bits+s_factor+su_pr_bv_bits+20+(pt_bits*n_children)));
        //==

        //move to the next byte-aligned position
        size_t header_bytes = INT_CEIL(bit_pos, 8);
        assert((header_bytes*8)==header_bits);

        //==add the stream of the children
        buffer.concatenate(header_bytes, children_buffer, node_n_bits/8);
        //==

        node_n_bits+=header_bits;

        //gather statistics
        stats.header_overhead+=header_bits;
        stats.rank_overhead+=rank_bits;
        stats.int_su_pr_overhead+=su_pr_bv_bits;
    }

    inline void compute_ext_succ_pred_info(){

        assert(lvl==0);

        //compute symbols that are in few trees
        std::vector<bool> low_freq_syms(bwt_rep.sigma, false);
        std::vector<std::pair<uint64_t, int64_t>> active_pred(bwt_rep.sigma, {0, -1});
        std::vector<std::pair<uint64_t, int64_t>> active_succ(bwt_rep.sigma, {0, 0});

        size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, b_size);//original number of blocks in the first level of the tree

        size_t acc_bits=0;
        for(size_t s=0;s<bwt_rep.sigma;s++){
            double per =  double(sigma_trees[s].size())/double(n_children);
            //symbol present in <=1% of the trees
            if(per<=0.01){
                acc_bits+= sigma_trees[s].size()*(sym_width(n_blocks)+ sym_width(bwt_rep.tot_syms));
                for(unsigned long long tree : sigma_trees[s]){
                    //std::cout<<"symbol:"<<s<<" tree:"<<tree<<" offset:"<<tree_offset[tree]<<std::endl;
                }
                //std::cout<<""<<std::endl;
            }
            low_freq_syms[s]= per<=0.01;
            
            sigma_trees[s].push_back(n_children);
            active_succ[s] = {0, sigma_trees[s][0]};
        }

        std::cout<<"Computing ext. succ/pred info"<<std::endl;

        size_t l_sym=0, r_sym=0;
        size_t l_sa_bound, r_sa_bound = C[r_sym+1]-1;

        size_t l_tree_bound;
        size_t r_tree_bound;
        for(int64_t b=0;b<n_children;b++){

            l_tree_bound = tree_offset[b];
            r_tree_bound = tree_offset[b+1]-1;

            while(C[l_sym+1]<l_tree_bound){
                ++l_sym;
            }
            l_sa_bound = C[l_sym];

            while(r_sa_bound<r_tree_bound){
                r_sa_bound = C[++r_sym+1]-1;
            }

            //std::cout<<l_tree_bound<<", "<<r_tree_bound<<" -> "<<l_sa_bound<<", "<<r_sa_bound<<" "<<l_sym<<"/"<<r_sym<<std::endl;

            size_t max_dist=0, pred_samp=0, succ_samp=0;
            for(size_t s=0;s<bwt_rep.sigma;s++){
                size_t sym_pos =(b*bwt_rep.sigma)+s;
                ext_pred_info[sym_pos]=ext_pred_info[sym_pos] && !low_freq_syms[s];
                if(ext_pred_info[sym_pos]){
                    int64_t dist = b-active_pred[s].second;
                    assert(dist>0 && dist<n_children);
                    bool out_of_range = (tree_offset[active_pred[s].second+1]-1)<l_sa_bound;

                    ext_pred_info[sym_pos] = dist>5 && !out_of_range;
                    if(ext_pred_info[sym_pos]){
                        if(dist>max_dist) max_dist = dist;
                        //std::cout<<"block:"<<b<<", symbol:"<<s<<", pred_tree:"<<p_tree<<" "<<dist<<std::endl;
                        pred_samp++;
                    }
                }

                if(sigma_trees[s][active_pred[s].first]==b){
                    active_pred[s].first++;
                    active_pred[s].second = b;
                }
            }

            //std::cout<<" ----- "<<std::endl;
            for(size_t s=0;s<bwt_rep.sigma;s++){

                if(active_succ[s].second==b){
                    active_succ[s].first++;
                    assert(active_succ[s].first<sigma_trees[s].size());
                    active_succ[s].second = sigma_trees[s][active_succ[s].first];
                }

                size_t sym_pos =(b*bwt_rep.sigma)+s;
                ext_succ_info[sym_pos] = ext_succ_info[sym_pos] && !low_freq_syms[s];
                if(ext_succ_info[sym_pos]){
                    int64_t dist = active_succ[s].second-b;
                    assert(dist>0 && dist<n_children);
                    bool out_of_range = tree_offset[active_succ[s].second]>r_sa_bound;
                    ext_succ_info[sym_pos] = dist>5 && !out_of_range;
                    if(ext_succ_info[sym_pos]){
                        //std::cout<<"block:"<<b<<", symbol:"<<s<<", succ_tree:"<<s_tree<<" "<<dist<<std::endl;
                        succ_samp++;
                    }
                }
            }

            //std::cout<<"\n max_dist:"<<max_dist<<std::endl;
            acc_bits+= (pred_samp+succ_samp) * sym_width(max_dist);
            acc_bits+= 2*bwt_rep.sigma;
            acc_bits+= sym_width(max_dist);

            stats.ext_pred_freq[pred_samp]++;
            stats.ext_succ_freq[succ_samp]++;
        }
        stats.header_overhead+=acc_bits;
        stats.ext_su_pr_overhead=acc_bits;
        node_n_bits+=acc_bits;
    }

    inline void finish_forest() {

        tree_offset[n_children]=bwt_rep.tot_syms;

        //free unused memory
        auto *node = tmp_node;
        while(node!= nullptr){
            node->destroy();
            node = node->tmp_node;
        }
        //

        compute_ext_succ_pred_info();

        //pointer information
        //TODO
        size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, b_size);//original number of blocks in the first level of the tree
        //
        //headers with the rank information
        //pt_bits indicates how many bits we use to encode pointers:
        //bits_node/8 is the pointer and b_runs indicate collapsed blocks
        size_t pt_bits = sym_width(node_n_bits/8) + sym_width(b_runs-1);
        size_t header_bits = pt_bits + (pt_bits*n_blocks);
        //byte align *this internal node
        header_bits = INT_CEIL(header_bits, 8)*8;
        node_n_bits+=header_bits;

        //gather statistics
        stats.tree_pointers_overhead+=header_bits;
        stats.header_overhead+=header_bits;
    }

    inline void create_leaf(std::vector<block_type>& blocks,
                            size_t n_blocks,
                            const size_t parent_sigma,
                            const std::vector<bool>& parent_sigma_bv,
                            const std::vector<uint64_t>& parent_rank_info
                            ){

        assert(node_n_bits==0);
        assert(lvl>0);

        size_t sym, len, n_runs=0;//, longest_run=0;
        for(size_t j=0;j<(blocks[0].size()-1);j++){
            sym = blocks[0][j].first;
            len = blocks[0][j].second;
            assert(sym<bwt_rep.sigma);
            node_sigma_bv[sym]=true;
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
                assert(sym<bwt_rep.sigma);
                node_sigma_bv[sym] = true;
                block_ranks[sym]+=len;
                n_runs++;
            }

            for(size_t j=0;j<(blocks[i].size()-1);j++){
                sym = blocks[i][j].first;
                len = blocks[i][j].second;
                assert(sym<bwt_rep.sigma);
                node_sigma_bv[sym]=true;
                block_ranks[sym]+=len;
            }
            n_runs+=blocks[i].size()-1;
        }

        //process the last run
        sym = blocks[n_blocks-1].back().first;
        len = blocks[n_blocks-1].back().second;
        assert(sym<bwt_rep.sigma);
        node_sigma_bv[sym]=true;
        block_ranks[sym]+=len;
        n_runs++;

        assert(n_runs<=bwt_dt_type::max_block_runs);
        //compute the block's alphabet size
        size_t tmp_s=0;
        node_sigma=0;
        for(auto && s : node_sigma_bv){
            packed_alphabet[tmp_s++]=node_sigma;
            node_sigma+=s;
        }

        size_t bfr_dist[9]={0}, bytes;
        for(size_t i=0;i<n_blocks;i++){
            for(auto & run : blocks[i]){
                run.first = packed_alphabet[run.first];
                //I need to use a fixed number of bits for the symbols (i.e., sym_width(node_sigma) bits)
                bytes = INT_CEIL((sym_width(node_sigma)+sym_width(run.second)), 8);
                bfr_dist[bytes]++;
            }
        }
        assert(bfr_dist[0]==0);

        size_t total_vbytes=0, max_vbytes=0;
        for(size_t b=1;b<9;b++){
            total_vbytes +=bfr_dist[b]*b;
            if(bfr_dist[b]>0) max_vbytes = b;
        }
        assert(max_vbytes>0 && max_vbytes<6);

        if(max_vbytes>1){
            //number of control masks of 1 byte for fast vbyte decoding;
            total_vbytes += INT_CEIL(n_runs, (8/sym_width(max_vbytes-1)));
        }

        //alternative encoding using a fixed number of bytes per run
        size_t total_fbytes = max_vbytes*n_runs;

        //byte encoding for the runs of this leaf
        size_t leaf_enc, run_bits;
        bool fix_len_enc=false;
        if(total_vbytes<total_fbytes){
            assert(max_vbytes>1);
            run_bits = total_vbytes*8;
            stats.runs_overhead += run_bits;
            leaf_enc=5+max_vbytes;//+5 is to difference them from fix-length blocks
        }else{
            run_bits = total_fbytes*8;
            stats.runs_overhead += run_bits;
            leaf_enc=max_vbytes;
            fix_len_enc = true;
        }

        //THE ENCODING STARTS HERE
        size_t r_width;
        //width in bits to store the range values
        if(lvl==1){
            //the leaf is the root of the tree
            //rank information: previous trees
            r_width = sym_width(bwt_rep.max_freq);
        } else {
            //the leaf is the child of an internal node
            //rank information: previous siblings
            r_width = sym_width(b_size*s_factor);
        }

        size_t leaf_enc_width=4;//we use 4 bits to encode the encoding type for the sequence of runs in this leaf
        size_t rank_bits = r_width*node_sigma;
        size_t header_bits = 1+leaf_enc_width+parent_sigma+rank_bits;
        header_bits = INT_CEIL(header_bits, 8)*8;//byte-aligned

        //allocate bytes for the information of this leaf
        buffer.reserve_in_bits(header_bits + run_bits);

        //start writing the in the buffer
        size_t bit_pos = 0;
        //1 bit (true) to indicate this node is a leaf
        buffer.write(bit_pos, bit_pos, 1);
        bit_pos++;

        //parent_sigma bits to encode the leaf's effective alphabet
        buffer.write(bit_pos, bit_pos+leaf_enc_width-1, leaf_enc);
        bit_pos+=leaf_enc_width;
        for(size_t i=0;i<bwt_rep.sigma;i++){
            if(parent_sigma_bv[i]){
                buffer.write(bit_pos, bit_pos, node_sigma_bv[i]);
                bit_pos++;
            }
        }
        //assert(bit_pos==(1+4+parent_sigma));

        //write the rank information
        for(size_t i=0;i<bwt_rep.sigma;i++){
            if(node_sigma_bv[i]){
                buffer.write(bit_pos, bit_pos+r_width-1, parent_rank_info[i]);
                bit_pos+=r_width;
            }
        }
        //assert(bit_pos==(1+4+parent_sigma+rank_bits));

        //the runs are byte-aligned
        size_t byte_pos = INT_CEIL(bit_pos, 8);
        assert((byte_pos*8)==header_bits);

        auto *byte_stream = (uint8_t *) buffer.stream;
        size_t written_bytes = insert_runs(blocks, n_blocks, &byte_stream[byte_pos], sym_width(node_sigma), max_vbytes , fix_len_enc);
        assert((written_bytes*8)==run_bits);

        node_n_bits = header_bits+run_bits;

        if(lm_tree_branch){
            need_ext_pred = node_sigma_bv;
        }

        if(rm_tree_branch){
            need_ext_succ = node_sigma_bv;
        }

        //gather statistics
        stats.header_overhead+=header_bits;
        stats.rank_overhead+=rank_bits;
        stats.eff_runs += n_runs;
        stats.rpl_freq[n_runs]++;
        stats.leaf_depth_freq[lvl-1]++;//lvl=0 is the forest, so it doesn't count. lvl=1 is a root of a tree
        stats.leaf_enc_freq[leaf_enc]++;//the encoding type for a leaf
    }

    size_t insert_runs(std::vector<block_type>& blocks, size_t n_blocks, uint8_t *stream,
                       size_t sigma_bytes, size_t max_bytes, bool fix_len_enc){

        size_t written_bytes=0;
        if(fix_len_enc){
            size_t enc_run;
            for(size_t i=0;i<n_blocks;i++){
                for(auto & run : blocks[i]){
                    enc_run =  run.second<<node_sigma | run.first;
                    memcpy(stream, &enc_run, max_bytes);
                    stream+=max_bytes;
                    written_bytes+=max_bytes;
                }
            }
        }else{
            uint8_t control_bits = sym_width(max_bytes-1);
            size_t vb_lens[8]={0};

            //we pack the stream in groups of (at most) 8 elements,
            // Each code uses (at most) 8 bytes, thus we need 8*8=64 tmp_bytes
            uint64_t code;
            uint8_t tmp_stream[64];
            uint8_t control = 0, acc_width=0, p=0, byte_pos=0;

            for(size_t i=0;i<n_blocks;i++){
                for(auto & run : blocks[i]){

                    vb_lens[p] = INT_CEIL((sigma_bytes+sym_width(run.second)), 8);
                    code = run.second<<node_sigma | run.first;
                    memcpy(&tmp_stream[byte_pos], &code, vb_lens[p]);
                    control |= (vb_lens[p] << acc_width);

                    byte_pos+=vb_lens[p];
                    p++;
                    acc_width+=control_bits;

                    if(acc_width+control_bits>8){
                        *stream=control;
                        stream++;
                        memcpy(stream, &tmp_stream[0], byte_pos);
                        stream+=byte_pos;
                        written_bytes+=byte_pos+1;
                        p=0;
                        control = 0;
                        acc_width = 0;
                        byte_pos = 0;
                    }
                }
            }

            if(control!=0){
                *stream=control;
                stream++;
                memcpy(stream, &tmp_stream[0], byte_pos);
                written_bytes+=byte_pos+1;
            }
        }
        return written_bytes;
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

    void print_node_info(){

        std::string pad = std::string(lvl+1, '\t');

        std::cout<<pad<<(leaf? "leaf" :"internal_node")<<std::endl;
        std::cout<<pad<<"lm_child:"<<lm_child<<std::endl;
        std::cout<<pad<<"rm_child:"<<rm_child<<std::endl;
        std::cout<<pad<<"covered_symbols:"<<cov_symbols<<std::endl;
        std::cout<<pad<<"child_rank:"<<child_rank<<std::endl;
        std::cout<<pad<<"node_sigma:"<<int(node_sigma)<<std::endl;
        std::cout<<pad<<"level:"<<int(lvl)<<std::endl;
        std::cout<<pad<<"alphabet:(";
        size_t p=0;
        for(size_t s=0;s<bwt_rep.sigma;s++){
            if(node_sigma_bv[s]){
                std::cout<<(p++>0 ? ", ":"")<<s;
            }
        }
        std::cout<<")"<<std::endl;
        std::cout<<pad<<"lm_branch:"<<lm_tree_branch<<std::endl;
        std::cout<<pad<<"rm_branch:"<<rm_tree_branch<<std::endl;

        if(!leaf){
            std::cout<<pad<<"inner succ/prec info:";
            size_t s_comp=0;
            for(size_t s=0;s<node_sigma_bv.size();s++){
                if(node_sigma_bv[s]){
                    std::cout<<s<<"=(";
                    for(size_t c=0;c<n_children;c++){
                        std::cout<<(c>0 ? ", ":"")<<int_succ_pred_info[s_comp][c];
                    }
                    std::cout<<") ";
                    s_comp++;
                }
            }
            std::cout<<""<<std::endl;
        }

        if(lvl==1){
            std::cout<<pad<<"ext pred:(";
            p=0;
            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(node_sigma_bv[s] && !need_ext_pred[s]){
                    std::cout<<(p++>0 ? ", ":"")<<s;
                }
            }
            std::cout<<")"<<std::endl;

            std::cout<<pad<<"ext succ:(";
            p=0;
            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(node_sigma_bv[s] && !need_ext_succ[s]) {
                    std::cout<<(p++>0 ?", ":"")<<s;
                }
            }
            std::cout<<")"<<std::endl;
        }
        std::cout<<""<<std::endl;
    }

    template<node_type type>//internal or leaf
    inline void create_node(size_t n_blocks) {

        assert(aligned<8>(node_n_bits));//check it is byte-aligned

        //lvl=0 means the forest, and tmp_node is then the root v of a tree.
        // Therefore, v is the leftmost and rightmost branches of the tree
        tmp_node->cov_symbols = b_size*n_blocks;
        tmp_node->lm_child = lvl==0 || consumed_syms == 0;
        tmp_node->rm_child = lvl==0 || (consumed_syms+ tmp_node->cov_symbols)==(b_size*s_factor);
        tmp_node->lm_tree_branch = lm_tree_branch && tmp_node->lm_child;
        tmp_node->rm_tree_branch = rm_tree_branch && tmp_node->rm_child;
        tmp_node->child_rank = n_children;
        tmp_node->leaf = type==LEAF;

        if constexpr (type==INTERNAL){
            assert(n_blocks==1);
            tmp_node->get_node_alphabet(active_blocks[0]);
            for(auto & run : active_blocks[0]){
                tmp_node->process_run(run.first, run.second);
            }
            tmp_node->finish_run_scan();
            tmp_node->finish_int_node(node_sigma, node_sigma_bv, block_ranks);
            stats.children_freq[tmp_node->n_children]++;
        }else{
            assert(n_blocks>=1);
            tmp_node->create_leaf(active_blocks, n_blocks, node_sigma, node_sigma_bv, block_ranks);
        }

        //add the rank information of the active child node (next_node) to the
        // parent's rank information
        for(size_t s=0;s<bwt_rep.sigma;s++){
            block_ranks[s] += tmp_node->block_ranks[s];
        }

        //get the symbols that need predecessor/successor information to other trees
        if(tmp_node->lm_tree_branch){
            for(size_t s=0;s<bwt_rep.sigma;s++){
                need_ext_pred[s] = need_ext_pred[s] & tmp_node->need_ext_pred[s];
            }
        }

        if(tmp_node->rm_tree_branch){
            for(size_t s=0;s<bwt_rep.sigma;s++){
                need_ext_succ[s] = need_ext_succ[s] & tmp_node->need_ext_succ[s];
            }
        }
        //
        assert(aligned<8>(node_n_bits+tmp_node->node_n_bits));

        if(lvl>0){
            //add internal successor/predecessor information for the current node
            //lvl=0 is the forest, so it does not count
            assert(n_children<s_factor);
            size_t s_comp=0;
            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(node_sigma_bv[s]){
                    int_succ_pred_info[s_comp++][n_children]=tmp_node->block_ranks[s]>0;
                    assert((tmp_node->block_ranks[s]>0) == tmp_node->node_sigma_bv[s]);
                }
            }
            assert(s_comp==node_sigma);

            //append the stream of tmp_node to the stream of children for this node
            children_buffer.reserve_in_bits(node_n_bits+tmp_node->node_n_bits);
            children_buffer.concatenate(node_n_bits/8, tmp_node->buffer, tmp_node->node_n_bits/8);
        } else {

            if(n_blocks>stats.max_n_blocks){
                stats.max_n_blocks = n_blocks;
            }

            tree_offset[n_children] = consumed_syms;

            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(tmp_node->node_sigma_bv[s]){
                    sigma_trees[s].push_back(n_children);
                }
                //symbols within this tree requiring predecessor information to trees on the left side
                ext_pred_info[(n_children*bwt_rep.sigma)+s] = !tmp_node->node_sigma_bv[s] || !tmp_node->need_ext_pred[s];

                //symbols within this tree requiring successor information to trees on the right side
                ext_succ_info[(n_children*bwt_rep.sigma)+s] = !tmp_node->node_sigma_bv[s] || !tmp_node->need_ext_succ[s];
            }

            //store to disk
            assert(ofs.tellp()==block_ptr[n_children]);
            assert(aligned<8>(tmp_node->node_n_bits));
            ofs.write((char *)tmp_node->buffer.stream, tmp_node->node_n_bits/8);
            stats.trees_overhead+=tmp_node->node_n_bits;
        }

        //print the node information for debugging purposes
        //tmp_node->print_node_info();
        //

        //the pointer to the active child node (next_node) should be aligned
        node_n_bits+=tmp_node->node_n_bits;
        block_ptr[++n_children]=(node_n_bits/8);
        consumed_syms+=tmp_node->cov_symbols;
        tmp_node->reset();
    }

    inline void reset(){
        memset(block_ranks.data(), 0, block_ranks.size()*sizeof(uint64_t));
        std::fill(node_sigma_bv.begin(), node_sigma_bv.end(), false);
        std::fill(need_ext_pred.begin(), need_ext_pred.end(), true);
        std::fill(need_ext_succ.begin(), need_ext_succ.end(), true);
        for(size_t s=0;s<bwt_rep.sigma;s++){
            std::fill(int_succ_pred_info[s].begin(), int_succ_pred_info[s].end(), false);
        }
        n_children = 0;
        node_n_bits = 0;
        node_sigma = 0;
        consumed_syms = 0;
        lm_tree_branch = lvl==0;//lvl=1 means the root of the tree
        lm_tree_branch = lvl==0;
        block_ptr[0] = 0;
    }
};

template<class bwt_dt_type, class node_type>
void forest_stats(bwt_dt_type& bwt_rep,  std::vector<node_type>& tmp_nodes, stat_collector<bwt_dt_type>& stats){

    for(size_t s=0;s<bwt_rep.sigma;s++){
        std::cout<<"\tsymbol "<<s<<", rank: "<<tmp_nodes[0].block_ranks[s]<<std::endl;
    }

    std::cout<<"Number_of_runs_in_a_leaf dist:"<<std::endl;
    assert(stats.rpl_freq[0]==0);
    for(size_t r=1;r<=bwt_dt_type::max_block_runs;r++){
        std::cout<<"\t"<<r<<" : "<<stats.rpl_freq[r]<<std::endl;
    }
    std::cout<<"Total number of runs versus original number of runs: "<<stats.eff_runs<<" / "<<bwt_rep.orig_runs<<std::endl;
    std::cout<<"Increase in the number of runs "<<(double(stats.eff_runs)/double(bwt_rep.orig_runs)-1)*100<<"%"<<std::endl;

    std::cout<<"Leaf_depth dist:"<<std::endl;
    size_t tot_leaves=0;
    for(size_t i=0;i<20;i++){
        tot_leaves+=stats.leaf_depth_freq[i];
    }
    for(size_t i=0;i<20;i++){
        if(stats.leaf_depth_freq[i]>0){
            std::cout<<"\t"<<i<<": "<<double(stats.leaf_depth_freq[i])/double(tot_leaves)<<std::endl;
        }
    }

    std::cout<<"Leaf_encoding dist:"<<std::endl;
    for(size_t i=0;i<20;i++){
        if(stats.leaf_enc_freq[i]!=0){
            if(i<=5){
                std::cout<<"\tFixed "<<i<<" bytes:\t\t\t\t"<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
            }else{
                std::cout<<"\tVariable-length 1-"<<i-5<<" bytes: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
            }
        }
    }

    std::cout<<"Number_of_children dist:"<<std::endl;
    size_t del_nodes=0, tot_nodes=0;
    for(size_t i=0;i<20;i++){
        if(stats.children_freq[i]!=0){
            std::cout<<"\t"<<i<<": "<<stats.children_freq[i]<<std::endl;
            del_nodes+=(bwt_rep.scale_factor-i)*stats.children_freq[i];
            tot_nodes+=stats.children_freq[i];
        }
    }

    std::cout<<"Ext_pred_info_dist"<<std::endl;
    for(size_t i=0;i<256;i++){
        if(stats.ext_pred_freq[i]){
            std::cout<<"  sym_samp:"<<i<<" freq:"<<stats.ext_pred_freq[i]<<std::endl;
        }
    }

    std::cout<<"Ext_succ_info_dist"<<std::endl;
    for(size_t i=0;i<256;i++){
        if(stats.ext_succ_freq[i]){
            std::cout<<"  sym_samp:"<<i<<" freq:"<<stats.ext_succ_freq[i]<<std::endl;
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
    std::cout<<"Space breakdown"<<std::endl;
    //std::cout<<"Written bytes without the tree's succ/pred info: "<<INT_CEIL((tmp_nodes[0].node_n_bits-bwt_rep.tree_su_pr_header_overhead), 8)<<std::endl;
    std::cout<<"\tRuns: "<<INT_CEIL(stats.runs_overhead, 8)<<" ("<<(double(stats.runs_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"%)"<<std::endl;
    std::cout<<"\tHeaders: "<<INT_CEIL(stats.header_overhead, 8)<<" ("<<(double(stats.header_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"%)"<<std::endl;
    std::cout<<"\t\tTree pointers: "<<INT_CEIL(stats.tree_pointers_overhead, 8)<<" ("<<(double(stats.tree_pointers_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"%)"<<std::endl;
    size_t ptr_bv_ov = stats.header_overhead - (stats.rank_overhead + stats.int_su_pr_overhead + stats.ext_su_pr_overhead+ stats.tree_pointers_overhead);
    std::cout<<"\t\tInt. Pointers and bitvectors: "<<INT_CEIL(ptr_bv_ov, 8)<<" ("<<(double(ptr_bv_ov)/double(tmp_nodes[0].node_n_bits))*100<<"%)"<<std::endl;
    std::cout<<"\t\tRank: "<<INT_CEIL(stats.rank_overhead, 8)<<" ("<<(double(stats.rank_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"%)"<<std::endl;
    std::cout<<"\t\tInt. succ/pred: "<<INT_CEIL(stats.int_su_pr_overhead, 8)<<" ("<<(double(stats.int_su_pr_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"%)"<<std::endl;
    std::cout<<"\t\tExt. succ/pred: "<<INT_CEIL(stats.ext_su_pr_overhead, 8)<<" ("<<(double(stats.ext_su_pr_overhead)/double(tmp_nodes[0].node_n_bits))*100<<"%)"<<std::endl;
    assert(stats.header_overhead+stats.runs_overhead==tmp_nodes[0].node_n_bits);
    std::cout<<"\t\tOverhead of the trees without ext. succ/pred info nor tree pointers: "<<INT_CEIL(stats.trees_overhead, 8)<<std::endl;
    std::cout<<"space_usage:"<<float(tmp_nodes[0].node_n_bits)/float(bwt_rep.tot_syms)<<" bps"<<std::endl;\
}

template<class bwt_dt_type>
void build_from_grlbwt(bwt_dt_type& bwt_rep, std::string& bwt_file){

    stat_collector<bwt_dt_type> dt_sts;
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
    size_t sigma=0, sigma_bits=0, max_freq=0;
    for(size_t i=0;i<sym_freqs.size();i++){
        if(sym_freqs[i]!=0){
            if(sym_freqs[i]>max_freq) max_freq = sym_freqs[i];
            sym_map[i] = sigma;
            sym_freqs[sigma++] = sym_freqs[i];
            sigma_bits+= sym_width(sym_freqs[i]);
        }
    }
    sym_freqs.resize(sigma+1);
    bwt_rep.sigma = sigma;
    bwt_rep.max_freq = max_freq;
    //bwt_rep.sigma_bits = sigma_bits;
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

    tmp_workspace tws("./", false);
    std::ofstream ofs(tws.get_file("trees"), std::ios::binary);

    std::vector<rl_node<bwt_dt_type>> tmp_nodes;
    tmp_nodes.reserve(bwt_rep.levels+1);
    size_t b_size = bwt_dt_type::block_size;
    for(size_t i=0;i<=bwt_rep.levels;i++){
        tmp_nodes.push_back(rl_node(i, b_size, bwt_rep, dt_sts, sym_freqs, ofs));
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
    std::cout<<"trees stored in "<<tws.get_file("trees")<<std::endl;
    bwt_rep.eff_runs = dt_sts.eff_runs;

    forest_stats(bwt_rep, tmp_nodes, dt_sts);
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
    size_t tot_syms=0;
    size_t sigma=0;
    size_t max_freq=0;
    size_t eff_runs=0;

    size_t levels=0;
    size_t *stream=nullptr;

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