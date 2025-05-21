//
// Created by Diaz, Diego on 17.4.2025.
//

#ifndef RLBWT_VLB_CONSTRUCT_RLBWT_VLB_H
#define RLBWT_VLB_CONSTRUCT_RLBWT_VLB_H

#include "rlbwt_vlb.h"
#include "bwt_io.h"
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
    uint64_t lfs_offset=0;//symbols with low frequency for which we encode the blocks where they occur explicitly
    uint64_t ext_suc_overhead=0;
    uint64_t int_su_pr_overhead=0;
    uint64_t tree_pointers_overhead=0;
    uint64_t trees_overhead=0;
    uint64_t max_n_blocks=0;
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
    size_t syms_before=0;//number of symbols in the text preceeding this node
    size_t child_rank=0;//this node is the child_rank of its parent
    size_t cov_symbols=0;//how many symbols of the input BWT does this node cover
    size_t child_mark_acc=0;
    size_t leaf_enc=0;

    bool lm_tree_branch=true;//true if this node in the leftmost branch of its tree
    bool rm_tree_branch=true;//true if this node in the rightmost branch of its tree
    bool lm_child=false;//true if this node is the leftmost child of its parent
    bool rm_child=false;//true if this node is the rightmost child of its parent
    bool leaf=false;//true if the node is a leaf

    const size_t lvl;//level of the subtree
    const size_t b_size;//block size for the level
    const size_t s_factor = bwt_dt_type::scale_factor;//shrinking factor for further subdivision
    const size_t b_runs = bwt_dt_type::max_block_runs;//maximum number of runs in a sequence of blocks
    const uint8_t run_width = bwt_dt_type::run_width;//number of bits we require to encode symbols in the range [0..b_runs]

    rl_node *tmp_node = nullptr;

    stream_type buffer;
    stream_type children_buffer;

    bwt_dt_type& bwt_rep; //data structure encoding the representation
    std::vector<block_type> active_blocks; //run-length compressed blocks conforming a tree node
    std::vector<bool> node_sigma_bv; //buffer to compute the leaf's effective alphabet
    std::vector<std::vector<bool>> int_succ_pred_info; //bits indicating internal successor/predecessor info for each symbol in the alphabet

    //next_ext_pred[s]=false (resp., next_ext_succ[s]=false) means the symbol s needs predecessor (resp, successor) information
    //these vectors only consider symbols that are *in* the alphabet of the node
    std::vector<bool> need_ext_succ;//symbols in the alphabet of the need that need external successor info (symbols not in the right branch)
    std::vector<bool> child_marks;//int a block of size b_size, it is a bit vector B[1..s_factor] that marks the original blocks of size b_size/s_factor in the collapsed blocks

    std::vector<uint64_t> block_ranks;//rank information we store in the header of every internal node
    std::vector<uint64_t> block_ptr;//pointers (byte offsets) to the node's children
    std::vector<uint8_t> packed_alphabet;//leaf's packed alphabet

    //this information is for the tree
    std::ostream *ofs= nullptr;
    std::ifstream *ifs= nullptr;
    std::vector<uint64_t> tree_offset;//number of symbols in the text before each tree

    // the vector sigma_trees[s], with s \in \Sigma, is a strictly increasing
    // sequence encoding the blocks in the tree containing the symbol s
    std::vector<std::vector<uint64_t>> sigma_trees;

    //list of the symbols of each tree (as a bitvector) that require external successor information
    //That is, each symbol need_ext_succ[s] \cup the symbols not appearing in the tree
    //ext_pred_info and ext_succ_info are information for the nodes with level 1 (i.e., tree roots)
    std::vector<bool> ext_succ_info;
    //

    //a struct to collect statistics about the data structure
    stat_collector<bwt_dt_type>& stats;

    explicit rl_node(size_t _lvl, size_t _b_size, bwt_dt_type& _bwt_rep, stat_collector<bwt_dt_type>& st):
            lvl(_lvl),
            b_size(_b_size),
            bwt_rep(_bwt_rep),
            active_blocks(b_runs),
            node_sigma_bv(bwt_rep.sigma, false),
            int_succ_pred_info(bwt_rep.sigma, std::vector<bool>(s_factor, false)),
            packed_alphabet(bwt_rep.sigma, 0),
            block_ranks(bwt_rep.sigma, 0),
            need_ext_succ(bwt_rep.sigma, true),
            stats(st){
        if(lvl==0){
            //number of blocks in the tree representation
            size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, b_size);
            block_ptr.resize(n_blocks+2);
            block_ptr[0] = 0;
            node_sigma = bwt_rep.sigma;
            node_sigma_bv = std::vector<bool>(bwt_rep.sigma, true);

            sigma_trees.resize(bwt_rep.sigma);
            for(size_t s=0;s<bwt_rep.sigma;s++){
                sigma_trees[s].reserve(n_blocks);
            }
            ext_succ_info = std::vector<bool>(n_blocks*bwt_rep.sigma, false);
            tree_offset = std::vector<uint64_t>(n_blocks+2, 0);
        } else{
            block_ptr.resize(s_factor+1);
            child_marks.resize(s_factor+1);
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
        destroy_vector(need_ext_succ);
        destroy_vector(block_ranks);
        destroy_vector(block_ptr);
        destroy_vector(tree_offset);
        destroy_vector(packed_alphabet);
        destroy_vector(child_marks);
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
            r_width = sym_width(b_size*s_factor*s_factor);
        }
        size_t rank_bits = r_width*node_sigma;
        //amount of bits for the int succ/pred info
        size_t su_pr_bv_bits = node_sigma*n_children;

        //number of bits we use to encode pointers to the children
        size_t pt_bits = sym_width(node_n_bits/8);
        assert(pt_bits<(1<<bwt_rep.int_pt_width));

        //amount of information (in bits) for the header
        size_t header_bits = 1+parent_sigma+rank_bits+s_factor+su_pr_bv_bits+bwt_rep.int_pt_width+(pt_bits*n_children);
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
        for(size_t i=0;i<s_factor;i++){
            buffer.write(bit_pos, bit_pos, child_marks[i]);
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
        buffer.write(bit_pos, bit_pos+bwt_rep.int_pt_width-1, pt_bits);
        bit_pos+=bwt_rep.int_pt_width;
        for(size_t i=0;i<n_children;i++){
            //NOTE consider the humber of bytes in the header when computing the byte position of a child.
            //i.e., child_byte_pos = h_bits/8 + block_ptr[c], where c is the child we need to find
            buffer.write(bit_pos, bit_pos+pt_bits-1, block_ptr[i]);
            bit_pos+=pt_bits;
        }
        assert(bit_pos==(1+parent_sigma+rank_bits+s_factor+su_pr_bv_bits+bwt_rep.int_pt_width+(pt_bits*n_children)));
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

    inline void record_low_freq_symbols(size_t& bit_pos, std::vector<bool>& low_freq_syms){

        assert(lvl==0);
        //compute how many (and which) symbols have low frequency
        //(i.e., symbol present in <=1% of the trees)
        //we will include the extra value INT_CEIL(tot_syms, block_size)*block_size for consistency
        size_t n_syms=0, n_elms=0, sym_bit_pos=bit_pos;
        for(size_t s=0;s<bwt_rep.sigma;s++){
            double per =  double(sigma_trees[s].size())/double(n_children);
            low_freq_syms[s]= per<=0.01;
            n_elms+=(sigma_trees[s].size()+1)*low_freq_syms[s];
            n_syms+=per<=0.01;
            bit_pos++;
        }
        //

        //a dummy limit>=tot_syms to access the dummy block (the one that is fake)
        size_t limit = INT_CEIL(bwt_rep.tot_syms, b_size)*b_size;
        size_t w=sym_width(limit);
        size_t elm_bits = n_elms*w;
        size_t c_bit_pos=bit_pos;//control bits

        //40 bits to encode the bit offset(s) in the stream where the x elements of s lie
        //NOTE: (offset(s+1)-offset(s))/(w) is equal to x
        bit_pos+= (n_syms+1)*40;
        buffer.reserve_in_bits(bit_pos+elm_bits);

        size_t data_start = bit_pos;
        for(size_t s=0;s<bwt_rep.sigma;s++){
            if(low_freq_syms[s]){

                std::cout<<"symbol:"<<s<<" c_bit_pos:"<<c_bit_pos<<" b_pos:"<<bit_pos<<std::endl;
                buffer.write(c_bit_pos, c_bit_pos+39, bit_pos);
                c_bit_pos+=40;
                for(unsigned long tree_id : sigma_trees[s]){
                    //encode the tree where s occurs
                    buffer.write(bit_pos, bit_pos+w-1, tree_offset[tree_id]);
                    std::cout<<"\t b_pos:"<<bit_pos<<" "<<tree_offset[tree_id]<<std::endl;
                    bit_pos+=w;
                    //std::cout<<tree_offset[tree_id]<<" ";
                }
                buffer.write(bit_pos, bit_pos+w-1, limit);
                std::cout<<"\t b_pos:"<<bit_pos<<" "<<limit<<"\n"<<std::endl;
                bit_pos+=w;
            }
            buffer.write(sym_bit_pos, sym_bit_pos, low_freq_syms[s]);
            sym_bit_pos++;
        }
        buffer.write(c_bit_pos, c_bit_pos+39, bit_pos);
        c_bit_pos+=40;
        assert(c_bit_pos==data_start);
        assert(bit_pos==(c_bit_pos+elm_bits));
    }

    inline void compute_ext_succ_pred_info(std::vector<uint64_t>& concat_ext_suc_pred_info,
                                           std::vector<bool>& low_freq_syms,
                                           std::vector<uint64_t>& bk_boundaries){
        assert(lvl==0);
        std::vector<uint64_t> trees_extra_bits(n_children, 0);

        //compute symbols that are in few trees
        std::vector<std::pair<uint64_t, int64_t>> active_succ(bwt_rep.sigma, {0, 0});
        for(size_t s=0;s<bwt_rep.sigma;s++){
            sigma_trees[s].push_back(n_children);
            active_succ[s] = {0, sigma_trees[s][0]};
        }

        std::cout<<"Computing ext. succ/pred info"<<std::endl;
        size_t l_sym=0, r_sym=0;
        size_t r_sa_bound = bk_boundaries[r_sym+1]-1;
        size_t l_tree_bound, r_tree_bound;
        int64_t max_dist=0;

        for(int64_t b=0;b<n_children;b++){

            l_tree_bound = tree_offset[b];
            r_tree_bound = tree_offset[b+1]-1;

            while(bk_boundaries[l_sym+1]<l_tree_bound){
                ++l_sym;
            }

            while(r_sa_bound<r_tree_bound){
                r_sa_bound = bk_boundaries[++r_sym+1]-1;
            }

            //std::cout<<l_tree_bound<<", "<<r_tree_bound<<" -> "<<l_sa_bound<<", "<<r_sa_bound<<" "<<l_sym<<"/"<<r_sym<<std::endl;
            size_t succ_samp=0;
            int64_t max_tree_dist=0, real_dist;

            for(size_t s=0;s<bwt_rep.sigma;s++){

                if(active_succ[s].second==b){
                    active_succ[s].first++;
                    assert(active_succ[s].first<sigma_trees[s].size());
                    active_succ[s].second = sigma_trees[s][active_succ[s].first];
                }

                size_t sym_pos =(b*bwt_rep.sigma)+s;
                ext_succ_info[sym_pos] = ext_succ_info[sym_pos] && !low_freq_syms[s];
                if(ext_succ_info[sym_pos]){
                    int64_t tree_dist = active_succ[s].second-b;
                    assert(tree_dist>0 && tree_dist<n_children);
                    bool out_of_range = tree_offset[active_succ[s].second]>r_sa_bound;
                    ext_succ_info[sym_pos] = tree_dist>5 && !out_of_range;
                    if(ext_succ_info[sym_pos]){
                        real_dist = (tree_offset[b+tree_dist]-tree_offset[b])/b_size;
                        //if(tree_dist>max_tree_dist) max_tree_dist = tree_dist;
                        if(real_dist>max_tree_dist) max_tree_dist = real_dist;
                        //if(tree_offset[b]==285605888){
                        //    std::cout<<"what wea? block:"<<b<<", symbol:"<<s<<", succ_tree:"<<real_dist<<" lo guarde en:"<<concat_ext_suc_pred_info.size()<<std::endl;
                        //}
                        concat_ext_suc_pred_info.push_back(real_dist);
                        succ_samp++;
                    }
                }
            }

            concat_ext_suc_pred_info.push_back(max_tree_dist);
            if(max_tree_dist>max_dist) max_dist = max_tree_dist;
            size_t t_ext_bits = succ_samp*sym_width(max_tree_dist) + bwt_rep.sigma;
            trees_extra_bits[b] = t_ext_bits;
            stats.ext_succ_freq[succ_samp]++;
        }

        bwt_rep.mtd_bits = std::max<uint8_t>(1, sym_width(sym_width((size_t)max_dist)));
        size_t acc_bits=0;
        for(size_t b=0;b<n_children;b++){
            trees_extra_bits[b]+=bwt_rep.mtd_bits;
            acc_bits += INT_CEIL(trees_extra_bits[b], 8)*8;
        }

        stats.header_overhead+=acc_bits;
        stats.ext_suc_overhead+=acc_bits;
        node_n_bits+=acc_bits;
    }

    inline size_t compute_dummy_tree_bits(){
        //the trailing bits include:
        //bwt_rep.sigma bits indicating ext succ info (all set to false). We need these bits for consistency.
        //bwt_rep.mtd_bits indicate the width of the ext. succ info (fake)
        //1 bit indicates this is a (dummy) leaf
        //bwt_rep.sigma indicate which symbols have rank info (all set to true)
        //bwt_rep.sigma*sym_width(bwt_rep.max_freq) to encode the ranks answers
        size_t r_width = sym_width(bwt_rep.max_freq);
        size_t rank_bits = r_width*bwt_rep.sigma;
        size_t ext_succ_bits = bwt_rep.sigma + bwt_rep.mtd_bits;
        ext_succ_bits = INT_CEIL(ext_succ_bits, 8)*8;//trees are byte aligned by construction in the stream
        size_t trailing_bits = ext_succ_bits + 1 + bwt_rep.sigma + rank_bits;
        return INT_CEIL(trailing_bits, 8)*8;//byte-align the trailing bits
    }

    inline void append_dummy_tree(size_t& bit_pos, std::vector<uint64_t>& parent_rank_info){

        size_t trailing_bits = compute_dummy_tree_bits();
        size_t r_width = sym_width(bwt_rep.max_freq);
        size_t rank_bits = r_width*bwt_rep.sigma;

        //fake succ/pred info
        for(size_t i=0;i<(bwt_rep.sigma);i++){
            buffer.write(bit_pos, bit_pos, 0);
            bit_pos++;
        }

        buffer.write(bit_pos, bit_pos+bwt_rep.mtd_bits-1, 0);
        bit_pos+=bwt_rep.mtd_bits;
        bit_pos = INT_CEIL(bit_pos, 8)*8;

        //start writing rank info
        //1 bit (true) to indicate this node is a leaf
        buffer.write(bit_pos, bit_pos, 1);
        bit_pos++;

        //parent_sigma bits to encode the alphabet
        for(size_t i=0;i<bwt_rep.sigma;i++){
            buffer.write(bit_pos, bit_pos, 1);
            bit_pos++;
        }

        //write the rank information
        for(size_t i=0;i<bwt_rep.sigma;i++){
            //std::cout<<"symbol:"<<i<<" rank_bit_pos:"<<bit_pos<<" r_width:"<<int(r_width)<<" rank:"<<parent_rank_info[i]<<std::endl;
            buffer.write(bit_pos, bit_pos+r_width-1, parent_rank_info[i]);
            bit_pos+=r_width;
        }
        bit_pos= INT_CEIL(bit_pos, 8)*8;

        stats.header_overhead+=trailing_bits;
        stats.rank_overhead+=rank_bits;
    }

    inline void finish_tree(std::vector<uint64_t>& bk_boundaries) {

        assert(lvl==0);
        //round the last offset for consistency in the succ. info computation
        tree_offset[n_children]= INT_CEIL(bwt_rep.tot_syms, b_size)*b_size;
        tree_offset[n_children+1]= tree_offset[n_children];

        /*for(size_t i=0;i<bwt_rep.sigma;i++){
            std::cout<<"fake leaf symbol:"<<i<<" rank:"<<block_ranks[i]<<std::endl;
        }*/

        //compute symbols with low frequency and their tree positions explicitly
        size_t bit_pos=0;
        std::vector<bool> low_freq_syms(bwt_rep.sigma, false);
        record_low_freq_symbols(bit_pos, low_freq_syms);
        bwt_rep.lfs_bits = bit_pos;
        //

        //compute ext succ/pred information
        std::vector<uint64_t> concat_exp_suc_pred_info;
        compute_ext_succ_pred_info(concat_exp_suc_pred_info, low_freq_syms, bk_boundaries);
        //

        size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, b_size);//original number of blocks in the first level of the tree

        //pt_bits indicates how many bits we use to encode pointers to the trees:
        //node_n_bits/8 is the pointer and b_runs indicate collapsed blocks
        //so far, node_n_bits considers:
        // * the sum of the tree sizes in bits (excluding the external su/pred information)
        // * the sum of the ext. succ/pred information for the trees
        //bwt_rep.ext_pt_width = sym_width(node_n_bits/8) + run_width;
        bwt_rep.ext_pt_width = std::max<uint16_t>(sym_width(node_n_bits/8), 2*run_width)+1;

        //(pt_bits*n_blocks) for the pointers to the trees
        //+1 because we add a dummy tree at end for consistency
        size_t tree_ptr_bits = (bwt_rep.ext_pt_width*(n_blocks+1));
        size_t header_bits = bwt_rep.lfs_bits+tree_ptr_bits;
        bwt_rep.header_bytes = INT_CEIL(header_bits, 8);
        header_bits = bwt_rep.header_bytes*8;
        bit_pos=header_bits;

        size_t trailing_bits = compute_dummy_tree_bits();

        size_t w2;
        buffer.reserve_in_bits(header_bits+node_n_bits+trailing_bits);

        size_t s_sym_pos=0, s_trees, pos=0, max_tree_dist, stream_byte_pos, tree_new_byte_pos;
        std::streamsize tree_bytes;
        auto *stream = (char *)buffer.stream;
        for(size_t b=0;b<n_children;b++){

            assert(aligned<8>(bit_pos));
            tree_new_byte_pos = (bit_pos-header_bits)/8;

            /*if(tree_offset[b]==285605888){
                std::cout<<"This is the bit pos where I start to write: "<<bit_pos<<std::endl;
            }*/

            //add the ext successor information
            s_trees=0;
            for(size_t s=0;s<bwt_rep.sigma;s++){
                buffer.write(bit_pos, bit_pos, ext_succ_info[s_sym_pos]);
                s_trees+=ext_succ_info[s_sym_pos];
                bit_pos++;
                s_sym_pos++;
            }

            max_tree_dist = concat_exp_suc_pred_info[pos+s_trees];
            w2 = sym_width(max_tree_dist);
            buffer.write(bit_pos, bit_pos+bwt_rep.mtd_bits-1, w2);
            bit_pos+=bwt_rep.mtd_bits;

            for(size_t t=0;t<s_trees;t++){
                buffer.write(bit_pos, bit_pos+w2-1, concat_exp_suc_pred_info[pos++]);
                bit_pos+=w2;
            }

            pos++;
            bit_pos = INT_CEIL(bit_pos, 8)*8;
            //

            //add the tree
            stream_byte_pos = bit_pos/8;
            tree_bytes = block_ptr[b+1]-block_ptr[b];
            assert(INT_CEIL(stream_byte_pos+tree_bytes, 8) <= buffer.stream_cap);
            ifs->read(&stream[stream_byte_pos], tree_bytes);
            //
            bit_pos+=tree_bytes*8;
            block_ptr[b]=tree_new_byte_pos;
        }

        assert((bit_pos-header_bits)==node_n_bits);
        tree_new_byte_pos = (bit_pos-header_bits)/8;
        block_ptr[n_children]=tree_new_byte_pos;
        append_dummy_tree(bit_pos, block_ranks);//add a dummy block

        assert(pos==concat_exp_suc_pred_info.size());
        assert(bit_pos==(header_bits+node_n_bits+trailing_bits));
        block_ptr[n_children+1]=(bit_pos-header_bits)/8;

        //write the pointers to the trees
        bit_pos=bwt_rep.lfs_bits;
        size_t c=0, l=0, n_syms, r, offsets;
        for(size_t b=0;b<=n_children;b++){

            buffer.write(bit_pos, bit_pos+bwt_rep.ext_pt_width-1, (block_ptr[b]<<1));
            //std::cout<<"block:"<<b<<" real_block:"<<c<<" b_pos:"<<bit_pos<<" ptr:"<<block_ptr[b]<<" tree_offset:"<<tree_offset[b]<<" "<<tree_offset[b+1]<<std::endl;
            bit_pos+=bwt_rep.ext_pt_width;
            r = (tree_offset[b+1]-tree_offset[b])/b_size;
            c++;
            l++;
            r--;
            n_syms =tree_offset[b];
            while((n_syms+b_size)<tree_offset[b+1]){
                offsets = (l<<run_width) | r;
                offsets = (offsets<<1) | 1;
                //std::cout<<"block:"<<b<<" real_block:"<<c<<" b_pos:"<<bit_pos<<" offsets:"<<l<<" "<<r<<std::endl;
                buffer.write(bit_pos, bit_pos+bwt_rep.ext_pt_width-1, offsets);
                bit_pos+=bwt_rep.ext_pt_width;
                n_syms +=b_size;
                c++;
                l++;
                r--;
            }
            l=0;
        }
        assert(bit_pos==(bwt_rep.lfs_bits+tree_ptr_bits));
        assert(c==(n_blocks+1));

        node_n_bits+=header_bits;
        node_n_bits+=trailing_bits;
        //move the information to the bwt
        bwt_rep.stream.swap(buffer);

        //gather statistics
        stats.ext_suc_overhead+=bwt_rep.lfs_bits;
        stats.header_overhead+=header_bits;
        stats.tree_pointers_overhead+=tree_ptr_bits;
        stats.lfs_offset+=bwt_rep.lfs_bits;
    }

    void inline get_max_psum(uint64_t* tmp_psum, uint64_t* max_psum) const {
        //max prefix sum within 4 blocks with 8 runs each
        uint64_t tmp_max = std::max(std::max(tmp_psum[0], tmp_psum[1]),
                                    std::max(tmp_psum[2], tmp_psum[3]));
        if(tmp_max>max_psum[0]){
            max_psum[0] = tmp_max;
        }

        //max prefix sum within 2 blocks with 16 runs each
        tmp_psum[0]+=tmp_psum[1];
        tmp_psum[2]+=tmp_psum[3];
        tmp_max = std::max(tmp_psum[0], tmp_psum[2]);
        if(tmp_max>max_psum[1]){
            max_psum[1] = tmp_max;
        }

        //max prefix sum within one blocks withn 32 runs
        tmp_psum[0]+=tmp_psum[2];
        if(tmp_max>max_psum[2]){
            max_psum[2] = tmp_psum[0];
        }
    }

    inline uint8_t compute_leaf_enc_code(uint8_t max_bytes, bool vbyte_enc, const uint64_t *max_psum) const {

        uint8_t code = 0;

        switch (max_bytes) {
            case 1:
                code +=max_psum[1]>0xFF;
                code +=max_psum[2]>0XFF;
                return code;
            case 2:
                code =3;
                code +=max_psum[0]>0xFFFF;
                code +=max_psum[1]>0xFFFF;
                code <<=vbyte_enc;
                return code;
            case 3:
                code =9;
                code +=vbyte_enc;
                return code;
            case 4:
                code =11;
                code +=vbyte_enc;
                return code;
            case 5:
                return 13;
            case 6:
                return 14;
            case 7:
                return 15;
            default:
                std::cout<<"Unknown code for leaf encoding"<<std::endl;
                exit(1);
        }
    }

    inline void create_leaf(std::vector<block_type>& blocks,
                            size_t n_blocks,
                            const size_t parent_sigma,
                            const std::vector<bool>& parent_sigma_bv,
                            const std::vector<uint64_t>& parent_rank_info){

        //for 1 byte: check that the sum of 16 (or 32 for AVX) consecutive run lens is <=256
        //for 2 bytes: check that the sum of 8 (or 16 for AVX) consecutive run lens is <=2^16-1
        //for 3-4 bytes: check that the sum of 4 (or 8 for AVX) consecutive run lens is <=2^32-1
        //for 4-8 bytes: check that the sum of 2 (or 4 for AVX) consecutive run lens is <=2^64-1
        //we need to do this check to work with SIMD instructions. It it does not fit, use the next encoding that can fits them

        assert(node_n_bits==0);
        assert(lvl>0);

        size_t sym, len, n_runs=0;
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
        uint64_t tmp_psum[4]={0};//8,16,24,32
        uint64_t max_psum[3]={0};//8,16,32
        uint64_t bk = 0;

        for(size_t i=0;i<n_blocks;i++){
            for(auto & run : blocks[i]){
                run.first = packed_alphabet[run.first];
                //I need to use a fixed number of bits for the symbols (i.e., sym_width(node_sigma) bits)
                bytes = INT_CEIL((sym_width(node_sigma)+sym_width(run.second)), 8);
                bfr_dist[bytes]++;
                tmp_psum[bk>>3] += run.second;
                bk++;

                if(bk==32){
                    get_max_psum(tmp_psum, max_psum);
                    memset(tmp_psum, 0, 32);
                    bk=0;
                }
            }
        }
        get_max_psum(tmp_psum, max_psum);
        assert(bfr_dist[0]==0);

        size_t total_vbytes=0, max_bytes=0;
        for(size_t b=1;b<9;b++){
            total_vbytes +=bfr_dist[b]*b;
            if(bfr_dist[b]>0) max_bytes = b;
        }
        assert(max_bytes>0 && max_bytes<6);

        if(max_bytes>1){
            //number of control masks of 1 byte for fast vbyte decoding;
            total_vbytes += INT_CEIL(n_runs, (8/sym_width(max_bytes-1)));
        }
        //alternative encoding using a fixed number of bytes per run
        size_t total_fbytes = max_bytes*n_runs;

        //byte encoding for the runs of this leaf
        size_t run_bits;
        bool fix_len_enc=false;
        if(total_vbytes<total_fbytes){
            assert(max_bytes>1);
            run_bits = total_vbytes*8;
        }else{
            run_bits = total_fbytes*8;
            fix_len_enc = true;
        }

        stats.runs_overhead += run_bits;
        leaf_enc = compute_leaf_enc_code(max_bytes, !fix_len_enc, max_psum);

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
            r_width = sym_width(b_size*s_factor*s_factor);
        }

        size_t rank_bits = r_width*node_sigma;
        size_t header_bits = 1+bwt_rep.leaf_enc_width+parent_sigma+rank_bits;
        header_bits = INT_CEIL(header_bits, 8)*8;//byte-aligned

        //allocate bytes for the information of this leaf
        buffer.reserve_in_bits(header_bits + run_bits);

        //start writing the in the buffer
        size_t bit_pos = 0;
        //1 bit (true) to indicate this node is a leaf
        buffer.write(bit_pos, bit_pos, 1);
        bit_pos++;

        //parent_sigma bits to encode the leaf's effective alphabet
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

        //store the leaf encoding
        buffer.write(bit_pos, bit_pos+bwt_rep.leaf_enc_width-1, leaf_enc);
        bit_pos+=bwt_rep.leaf_enc_width;
        //

        //the runs are byte-aligned
        size_t byte_pos = INT_CEIL(bit_pos, 8);
        assert((byte_pos*8)==header_bits);

        auto *byte_stream = (uint8_t *) buffer.stream;
        size_t written_bytes = insert_runs(blocks, n_blocks, &byte_stream[byte_pos], sym_width(node_sigma), max_bytes , fix_len_enc);
        assert((written_bytes*8)==run_bits);

        node_n_bits = header_bits+run_bits;

        if(rm_tree_branch){
            need_ext_succ = node_sigma_bv;
        }
        bwt_rep.eff_runs += n_runs;

        //gather statistics
        if(n_blocks>stats.max_n_blocks) stats.max_n_blocks = n_blocks;
        stats.header_overhead+=header_bits;
        stats.rank_overhead+=rank_bits;
        stats.rpl_freq[n_runs]++;
        stats.leaf_depth_freq[lvl-1]++;//lvl=0 is the tree, so it doesn't count. lvl=1 is a root of a block
        stats.leaf_enc_freq[leaf_enc]++;//the encoding type for a leaf
    }

    size_t insert_runs(std::vector<block_type>& blocks, size_t n_blocks, uint8_t *stream, size_t sigma_bits,
                       size_t max_bytes, bool fix_len_enc){

        assert(sigma_bits== sym_width(node_sigma));

        size_t written_bytes=0;
        if(fix_len_enc){
            size_t enc_run;
            for(size_t i=0;i<n_blocks;i++){
                for(auto & run : blocks[i]){
                    enc_run =  run.second<<sigma_bits | run.first;
                    memcpy(stream, &enc_run, max_bytes);
                    stream+=max_bytes;
                    written_bytes+=max_bytes;
                }
            }
        }else{
            uint8_t n_ctrl_bits = sym_width(max_bytes-1);
            size_t vb_lens[8]={0};

            //we pack the stream in groups of (at most) 8 elements,
            // Each code uses (at most) 8 bytes, thus we need 8*8=64 tmp_bytes
            uint64_t code;
            uint8_t tmp_stream[64];
            uint8_t ctrl_bits = 0, acc_width=0, p=0, byte_pos=0;

            for(size_t i=0;i<n_blocks;i++){
                for(auto & run : blocks[i]){

                    vb_lens[p] = INT_CEIL((sigma_bits+sym_width(run.second)), 8);
                    code = run.second<<sigma_bits | run.first;
                    memcpy(&tmp_stream[byte_pos], &code, vb_lens[p]);
                    ctrl_bits |= ((vb_lens[p]-1U) << acc_width);
                    byte_pos+=vb_lens[p];
                    p++;
                    acc_width+=n_ctrl_bits;

                    if(acc_width+n_ctrl_bits>8){
                        *stream=ctrl_bits;
                        stream++;

                        memcpy(stream, &tmp_stream[0], byte_pos);
                        stream+=byte_pos;
                        written_bytes+=byte_pos+1;
                        p=0;
                        ctrl_bits = 0;
                        acc_width = 0;
                        byte_pos = 0;
                    }
                }
            }

            if(p!=0){
                *stream=ctrl_bits;
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

    void print_node_info(std::vector<block_type> bkl, size_t n_blocks, std::vector<uint64_t>& parent_rank_info){

        std::string pad = std::string(lvl+1, '\t');

        std::cout<<pad<<(leaf? "leaf" :"internal_node")<<std::endl;
        std::cout<<pad<<"lm_child:"<<lm_child<<std::endl;
        std::cout<<pad<<"rm_child:"<<rm_child<<std::endl;
        std::cout<<pad<<"sym_before:"<<syms_before<<std::endl;
        std::cout<<pad<<"covered_symbols:"<<cov_symbols<<std::endl;
        std::cout<<pad<<"child_rank:"<<child_rank<<std::endl;
        std::cout<<pad<<"node_sigma:"<<int(node_sigma)<<std::endl;
        std::cout<<pad<<"level:"<<int(lvl)<<std::endl;
        std::cout<<pad<<"size in bytes:"<<INT_CEIL(node_n_bits, 8)<<std::endl;
        std::cout<<pad<<"alphabet:(";
        size_t p=0;
        for(size_t s=0;s<bwt_rep.sigma;s++){
            if(node_sigma_bv[s]){
                std::cout<<(p++>0 ? ", ":"")<<"(sym:"<<s<<" rank:"<<parent_rank_info[s]<<")";
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
            std::cout<<pad<<"n_children:"<<n_children<<std::endl;
            std::cout<<pad<<"child info:";
            for(size_t k=0;k<s_factor;k++){
                std::cout<<(k>0 ? ", ":"")<<int(child_marks[k]);
            }
            std::cout<<""<<std::endl;
            std::cout<<pad<<"pointer to children:";
            for(size_t k=0;k<n_children;k++){
                std::cout<<(k>0 ? ", ":"")<<block_ptr[k];
            }
            std::cout<<""<<std::endl;
        }else{
            std::cout<<pad<<"leaf encoding:"<<int(leaf_enc)<<std::endl;
            std::cout<<pad<<"runs: ";
            for(size_t k=0;k<n_blocks;k++){
                for(auto & l : bkl[k]){
                    std::cout<<"(packed_sym:"<<int(l.first)<<",len:"<<l.second<<") ";
                }
            }
            std::cout<<""<<std::endl;
        }

        if(lvl==1){
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
        //lvl=0 means the tree root, and tmp_node is then the root v of a block in the tree.
        // Therefore, v is the leftmost and rightmost branches of the tree
        tmp_node->cov_symbols = b_size*n_blocks;
        tmp_node->lm_child = lvl==0 || consumed_syms == 0;
        tmp_node->rm_child = lvl==0 || (consumed_syms+ tmp_node->cov_symbols)==(b_size*s_factor);
        tmp_node->lm_tree_branch = lm_tree_branch && tmp_node->lm_child;
        tmp_node->rm_tree_branch = rm_tree_branch && tmp_node->rm_child;
        tmp_node->child_rank = n_children;
        tmp_node->leaf = type==LEAF;
        tmp_node->syms_before=syms_before+consumed_syms;

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

        //print the node information for debugging purposes
        //tmp_node->print_node_info(active_blocks, n_blocks, block_ranks);
        //

        //add the rank information of the active child node (tmp_node) to the
        // parent's rank information
        for(size_t s=0;s<bwt_rep.sigma;s++){
            block_ranks[s] += tmp_node->block_ranks[s];
        }

        //get the symbols that need successor information to other trees
        if(tmp_node->rm_tree_branch){
            for(size_t s=0;s<bwt_rep.sigma;s++){
                need_ext_succ[s] = need_ext_succ[s] & tmp_node->need_ext_succ[s];
            }
        }
        //
        assert(aligned<8>(node_n_bits+tmp_node->node_n_bits));

        if(lvl>0){

            child_marks[child_mark_acc]=true;
            child_mark_acc+=n_blocks;
            assert(child_mark_acc<=s_factor);

            //add internal successor/predecessor information for the current node
            //lvl=0 is the tree root, so it does not count
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

            tree_offset[n_children] = consumed_syms;

            for(size_t s=0;s<bwt_rep.sigma;s++){
                if(tmp_node->node_sigma_bv[s]){
                    sigma_trees[s].push_back(n_children);
                }
                //symbols within this tree requiring successor information to trees on the right side
                ext_succ_info[(n_children*bwt_rep.sigma)+s] = !tmp_node->node_sigma_bv[s] || !tmp_node->need_ext_succ[s];
            }

            //store to disk
            assert(ofs->tellp()==block_ptr[n_children]);
            assert(aligned<8>(tmp_node->node_n_bits));
            ofs->write((char *)tmp_node->buffer.stream, tmp_node->node_n_bits/8);
            stats.trees_overhead+=tmp_node->node_n_bits;
        }

        //the pointer to the active child node (next_node) should be aligned
        node_n_bits+=tmp_node->node_n_bits;
        block_ptr[++n_children]=(node_n_bits/8);
        consumed_syms+=tmp_node->cov_symbols;
        tmp_node->reset();
    }

    inline void reset(){
        memset(block_ranks.data(), 0, block_ranks.size()*sizeof(uint64_t));
        std::fill(child_marks.begin(), child_marks.end(), false);
        std::fill(node_sigma_bv.begin(), node_sigma_bv.end(), false);
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
        child_mark_acc=0;
        leaf_enc=0;
        syms_before=0;
    }
};

template<class bwt_dt_type, class node_type>
struct tree_dt{

    stat_collector<bwt_dt_type> stats;
    tmp_workspace twd;
    bwt_dt_type& bwt_rep;
    std::vector<uint64_t> C;
    node_type *root= nullptr;

    explicit tree_dt(const std::string& tmp_dir, bwt_dt_type& _bwt_rep): twd(tmp_dir),
                                                                         bwt_rep(_bwt_rep),
                                                                         C(256, 0){}

    void build_from_grlbwt(std::string& bwt_file){

        bwt_buff_reader bwt_buff(bwt_file);
        size_t n_runs = bwt_buff.size();
        bwt_rep.orig_runs = n_runs;
        bwt_rep.packed_alpha.resize(256);

        //compute symbol frequencies
        size_t sym, len;
        for(size_t i=0;i<n_runs;i++){
            bwt_buff.read_run(i, sym, len);
            C[sym]+=len;
        }
        //

        //compute the effective alphabet and the number of bits we require to encode it
        size_t sigma=0, sigma_bits=0, max_freq=0;
        bwt_rep.unpacked_alpha.reserve(256);
        for(size_t i=0;i<C.size();i++){
            if(C[i]!=0){
                if(C[i]>max_freq) max_freq = C[i];
                bwt_rep.packed_alpha[i] = sigma;
                bwt_rep.unpacked_alpha.push_back(i);
                C[sigma++] = C[i];
                sigma_bits+= sym_width(C[i]);
            }
        }
        bwt_rep.unpacked_alpha.shrink_to_fit();
        C.resize(sigma+1);
        //
        //compute the array C[1..\sigma]
        size_t acc=0, tmp;
        for(size_t i=0;i<sigma;i++){
            tmp = C[i];
            C[i]=acc;
            acc+=tmp;
        }
        C[sigma] = acc;
        bwt_rep.tot_syms = acc;
        bwt_rep.sigma = sigma;
        bwt_rep.max_freq = max_freq;
        /*for(size_t s=0;s<sigma;s++){
            std::cout<<s<<" "<<C[s]<<std::endl;
        }*/
        //

        root = new node_type(0, bwt_dt_type::block_size, bwt_rep, stats);
        std::ofstream ofs(twd.get_file("trees"), std::ios::binary);
        root->ofs = &ofs;

        node_type *current = root;
        size_t b_size=bwt_dt_type::block_size/bwt_dt_type::scale_factor;
        for(size_t i=1;i<=bwt_rep.levels;i++){
            current->tmp_node = new node_type(i, b_size, bwt_rep, stats);
            current = current->tmp_node;
            b_size/=bwt_dt_type::scale_factor;
        }

        //compute the tree
        for(size_t i=0;i<n_runs;i++){
            bwt_buff.read_run(i, sym, len);
            root->process_run(bwt_rep.packed_alpha[sym], len);
        }
        root->finish_run_scan();
        ofs.close();

        std::ifstream ifs(twd.get_file("trees"), std::ios::binary);
        root->ifs = &ifs;
        root->finish_tree(C);
        ifs.close();
        //
    }

    void build_from_rl_plain(std::string& bwt_file){

    }

    void build_from_plain(std::string& bwt_file){

    }

    void report_stats(){
        for(size_t s=0;s<bwt_rep.sigma;s++){
            std::cout<<"\tsymbol "<<s<<", rank: "<<root->block_ranks[s]<<std::endl;
        }

        std::cout<<"Number_of_runs_in_a_leaf dist:"<<std::endl;
        assert(stats.rpl_freq[0]==0);
        for(size_t r=1;r<=bwt_dt_type::max_block_runs;r++){
            std::cout<<"\t"<<r<<" : "<<stats.rpl_freq[r]<<std::endl;
        }
        std::cout<<"Total number of runs versus original number of runs: "<<bwt_rep.eff_runs<<" / "<<bwt_rep.orig_runs<<std::endl;
        std::cout<<"Increase in the number of runs "<<(double(bwt_rep.eff_runs)/double(bwt_rep.orig_runs)-1)*100<<"%"<<std::endl;

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
        for(size_t i=0;i<16;i++){
            switch(i) {
                case 0:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t1 byte: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 1:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t1 byte, overflow 16: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 2:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t1 byte, overflow 32: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 3:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 4:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, overflow 8: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 5:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, overflow 16: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 6:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 7:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, overflow 8, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 8:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, overflow 16, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 9:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t3 bytes: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 10:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t3 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 11:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t4 bytes: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 12:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t4 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 13:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t5 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 14:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t6 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 15:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t7 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
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

        std::cout<<"Ext_succ_info_dist"<<std::endl;
        for(size_t i=0;i<256;i++){
            if(stats.ext_succ_freq[i]){
                std::cout<<"  sym_samp:"<<i<<" freq:"<<stats.ext_succ_freq[i]<<std::endl;
            }
        }

        size_t n_blocks = INT_CEIL(bwt_rep.tot_syms, bwt_rep.block_size);//original number of blocks in the first level of the tree
        std::cout<<"Effective number of trees versus full number of trees (n/b): "<<root->n_children<<" / "<<n_blocks<<std::endl;
        std::cout<<"Percentage of removed trees: "<<(1-double(root->n_children)/double(n_blocks))*100<<"% "<<std::endl;

        tot_nodes*=bwt_rep.scale_factor;
        tot_nodes+=n_blocks;
        del_nodes=n_blocks-root->n_children;

        std::cout<<"Percentage of removed nodes: "<<(double(del_nodes)/double(tot_nodes))*100<<"% "<<std::endl;
        std::cout<<"Written bytes in the data structure: "<<INT_CEIL(root->node_n_bits, 8)<<std::endl;
        std::cout<<"Space breakdown"<<std::endl;
        std::cout<<"\tRuns: "<<INT_CEIL(stats.runs_overhead, 8)<<" ("<<(double(stats.runs_overhead)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        std::cout<<"\tHeaders: "<<INT_CEIL(stats.header_overhead, 8)<<" ("<<(double(stats.header_overhead)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        std::cout<<"\t\tTree pointers: "<<INT_CEIL(stats.tree_pointers_overhead, 8)<<" ("<<(double(stats.tree_pointers_overhead)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        size_t ptr_bv_ov = stats.header_overhead - (stats.rank_overhead + stats.int_su_pr_overhead + stats.ext_suc_overhead+ stats.tree_pointers_overhead);
        std::cout<<"\t\tInt. Pointers and bitvectors: "<<INT_CEIL(ptr_bv_ov, 8)<<" ("<<(double(ptr_bv_ov)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        std::cout<<"\t\tRank: "<<INT_CEIL(stats.rank_overhead, 8)<<" ("<<(double(stats.rank_overhead)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        std::cout<<"\t\tInt. succ/pred: "<<INT_CEIL(stats.int_su_pr_overhead, 8)<<" ("<<(double(stats.int_su_pr_overhead)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        std::cout<<"\t\tExt. succ: "<<INT_CEIL(stats.ext_suc_overhead, 8)<<" ("<<(double(stats.ext_suc_overhead)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        std::cout<<"\t\t\tLow freq. symbols: "<<INT_CEIL(stats.lfs_offset, 8)<<" ("<<(double(stats.lfs_offset)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        assert(stats.header_overhead+stats.runs_overhead==root->node_n_bits);
        std::cout<<"\t\tTrees without ext. succ/pred info nor tree pointers: "<<INT_CEIL(stats.trees_overhead, 8)<<std::endl;
        std::cout<<"space_usage:"<<float(root->node_n_bits)/float(bwt_rep.tot_syms)<<" bps"<<std::endl;
    }
};


template<class bwt_dt_type>
void build_rlbwt_vlb(bwt_dt_type& bwt_rep, std::string& bwt_file, INPUT_FORMAT f, std::string tmp_dir="./"){
    tree_dt<bwt_dt_type, rl_node<bwt_dt_type>> tree(tmp_dir, bwt_rep);
    switch (f) {
        case INPUT_FORMAT::GRL_BWT:
            tree.build_from_grlbwt(bwt_file);
            break;
        case INPUT_FORMAT::RL_PLAIN:
            tree.build_from_rl_plain(bwt_file);
            break;
        case INPUT_FORMAT::PLAIN:
            tree.build_from_plain(bwt_file);
            break;
        default:
            std::cout<<"Unknown format"<<std::endl;
            exit(1);
    }
    tree.report_stats();
}
#endif //RLBWT_VLB_CONSTRUCT_RLBWT_VLB_H
