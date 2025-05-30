//
// Created by Diaz, Diego on 30.5.2025.
//

#ifndef VLBT_COMMON_H
#define VLBT_COMMON_H
enum BWT_FORMAT{
    GRL_BWT=0,
    RL_PLAIN=1,
    PLAIN=2
};

enum node_type {
    INTERNAL,
    LEAF
};

//statistics about the data structure
template<class bwt_type>
struct bwt_stat_collector{
    uint64_t rpl_freq[bwt_type::max_block_runs+1]={0};//number of runs in a leaf
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
    uint64_t samp_overhead=0;
};
#endif //VLBT_COMMON_H
