//
// Created by Diaz, Diego on 3.2.2022.
//

#ifndef TEST_RL_BCR_BWT_FM_INDEX_H
#define TEST_RL_BCR_BWT_FM_INDEX_H

#include "bwt_io.h"

template<class bwt_type>
struct fm_index{

    bwt_type& bwt;
    std::vector<uint64_t> C;
    std::vector<uint8_t> char2comp;
    uint8_t dummy=0;

    explicit fm_index(bwt_type& _bwt, std::vector<uint64_t> _C): bwt(_bwt), C(_C){
    }

    [[nodiscard]] size_t size() const {
        return bwt.size();
    }

    [[nodiscard]] size_t alphabet() const {
        return C.size();
    }

    [[nodiscard]] std::pair<uint8_t, size_t> lf(size_t idx) const {
        auto res = bwt.inverse_select(idx);
        size_t next = C[char2comp[res.second]] + res.first;
        return {res.second, next};
    }

    [[nodiscard]] inline uint8_t get_dummy() const {
        return dummy;
    }

    [[nodiscard]] inline size_t tot_strings() const {
        return C[1]-C[0];
    }

    inline size_t backward_search(std::string &pat){
        size_t l=0, r=bwt.size()-1, j=pat.size();
        uint8_t cc;
        while(j-->0 && l<=r){
            cc = char2comp[pat[j]];
            l = C[cc] + bwt.rank(l, pat[j]); // count c in bwt[0..l-1]
            r = C[cc] + bwt.rank(r+1, pat[j]) - 1; // count c in bwt[0..r]
        }
        return r+1-l;
    }
};
#endif //TEST_RL_BCR_BWT_FM_INDEX_H
