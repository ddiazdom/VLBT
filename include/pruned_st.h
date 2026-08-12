/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */

#ifndef VLBT_PRUNED_ST_H
#define VLBT_PRUNED_ST_H

#include "utils.h"
#include "bwt_streams.h"
#include <stack>

//suffix tree node represented as the range of the leaves it covers
struct st_node_t {
    uint64_t start=0;
    uint64_t end=0;
    uint64_t depth=0;
    st_node_t()=default;
    st_node_t(uint64_t start_, uint64_t end_, uint64_t depth_) : start(start_),
                                                                 end(end_),
                                                                 depth(depth_){}

    bool is_child(st_node_t& v) const{
        return start<=v.start && v.end<=end;
    }

    bool is_sibling(st_node_t& v) const {
        return v.depth == depth && v.start==end+1;
    }

    bool unrelated(st_node_t& v) const {
        return v.depth<depth;
    }

    //+1 because the backward search considers rank(lb, sym) and rank(rb+1, sym)
    //If this st_node (start, end) immediately precedes (a, b) (that is, end+1=a),
    //then (a, b) might need `sym` to answer rank(a=rb+1, sym)
    [[nodiscard]] bool intersect(const uint64_t start_, const uint64_t end_) const {
        return !(end+1 < start_ || end_ < start);
    }
};

//pruned suffix tree where the internal nodes with depth>threshold are pruned and the whole subtree is treated as a leaf
struct pruned_suffix_tree{
    size_t k=0;
    std::vector<st_node_t>& nodes_in_dfs;
    std::stack<st_node_t> stack;

    explicit pruned_suffix_tree(std::vector<st_node_t>& nodes_): nodes_in_dfs(nodes_){
        stack.push(nodes_in_dfs[k++]);
    }

    [[nodiscard]] bool intersect(uint64_t start, uint64_t end) const {
        if(k>=nodes_in_dfs.size()) return false;
        return stack.top().intersect(start, end);
    }

    void operator++(){
        if(k<nodes_in_dfs.size()){
            assert(!stack.empty());

            if(stack.top().unrelated(nodes_in_dfs[k])){
                stack.pop();//return to the parent
            } else {
                if(stack.top().is_sibling(nodes_in_dfs[k])){
                    stack.pop();
                }
                stack.push(nodes_in_dfs[k]);
                k++;
            }
        } else if(k==nodes_in_dfs.size()){
            stack.pop();//climbing the rightmost branch of the pruned suffix tree
            if(stack.empty()){
                stack.emplace(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), 0);
                k++;
            }
        }
    }

    st_node_t operator*() const {
        return stack.top();
    }
};

//a block in the partition of VLBT represented as a range over the BWT
struct block_range_t{
    uint64_t start;
    uint64_t end;
    uint64_t lb;
    uint64_t rb;
};

inline std::vector<st_node_t> compute_nodes_of_pruned_st(rle_vbyte_bwt_reader& bwt_buff, const std::vector<uint64_t>& sa_ranges,
                                                         const std::vector<uint8_t>& packed_alpha, const size_t n_iterations){


    std::vector<std::vector<st_node_t>> st_nodes_in_dfs(n_iterations+2);

    size_t depth=0;
    st_nodes_in_dfs[depth].emplace_back(0, sa_ranges.back()-1, depth);
    size_t tot_nodes=1;
    depth++;

    for(size_t i=0;i<(sa_ranges.size()-1);i++){
        st_nodes_in_dfs[depth].emplace_back(sa_ranges[i], sa_ranges[i+1]-1, depth);
    }
    tot_nodes += st_nodes_in_dfs[depth].size();
    depth++;

    for(size_t i=0;i<n_iterations;i++){

        std::vector<uint64_t> symbols(256, 0);
        std::vector<uint64_t> block_symbols(256, 0);

        size_t k = 0;
        uint64_t range_end = st_nodes_in_dfs[depth-1][k].end+1;
        size_t acc=0;
        uint64_t sym, len;

        bwt_buff.rewind();
        while(bwt_buff.next_run(sym, len)){

            sym = packed_alpha[sym];

            assert(acc<range_end);

            while((acc+len)>=range_end){
                size_t l = range_end - acc;
                block_symbols[sym]+=l;
                acc+=l;

                for(size_t s=0;s<256;s++){
                    if(block_symbols[s]>0){
                        uint64_t lb = sa_ranges[s]+symbols[s];
                        uint64_t rb = lb+block_symbols[s]-1;
                        st_nodes_in_dfs[depth].emplace_back(lb, rb, depth);
                        symbols[s]+=block_symbols[s];
                    }
                }

                k++;
                range_end++;
                if(k<st_nodes_in_dfs[depth-1].size()){
                    range_end = st_nodes_in_dfs[depth-1][k].end+1;
                }
                len -=l;
                memset(block_symbols.data(), 0, 256*sizeof(uint64_t));
            }
            block_symbols[sym]+=len;
            acc+=len;
        }
        assert(acc==(range_end-1));

        std::sort(st_nodes_in_dfs[depth].begin(), st_nodes_in_dfs[depth].end(), [](auto const& a, auto const& b){
            return a.start<b.start;
        });

        tot_nodes+=st_nodes_in_dfs[depth].size();
        depth++;
    }

    st_nodes_in_dfs[0].reserve(tot_nodes);
    for(size_t d=1;d<depth;d++){
        for(size_t j=0;j<st_nodes_in_dfs[d].size();j++){
            st_nodes_in_dfs[0].push_back(st_nodes_in_dfs[d][j]);
        }
        destroy_vector(st_nodes_in_dfs[d]);
    }

    std::sort(st_nodes_in_dfs[0].begin(), st_nodes_in_dfs[0].end(), [](auto const& a, auto const& b){
        if(a.start==b.start){
            return a.depth<b.depth;
        }
        return a.start<b.start;
    });

    size_t k=0;
    for(size_t j=1;j<st_nodes_in_dfs[0].size();j++){
        if(st_nodes_in_dfs[0][j].start!=st_nodes_in_dfs[0][k].start ||
           st_nodes_in_dfs[0][j].end!=st_nodes_in_dfs[0][k].end){
            st_nodes_in_dfs[0][++k] = st_nodes_in_dfs[0][j];
        }
    }
    st_nodes_in_dfs[0].resize(k);
    st_nodes_in_dfs[0].shrink_to_fit();

    return st_nodes_in_dfs[0];
}

#endif //VLBT_PRUNED_ST_H
