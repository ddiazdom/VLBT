//
// Created by Diaz, Diego on 3.2.2022.
//

#ifndef TEST_RL_BCR_BWT_FM_INDEX_H
#define TEST_RL_BCR_BWT_FM_INDEX_H

#include "../include/bwt_io.h"

template<class bwt_type>
struct fm_index{

    bwt_type& bwt;
    std::vector<uint64_t> C;
    const std::vector<uint8_t>& byte2comp;
    const std::vector<uint8_t>& comp2byte;
    uint8_t dummy=0;

    explicit fm_index(bwt_type& bwt_,
                      std::vector<uint64_t>& C_,
                      const std::vector<uint8_t>& byte2comp_,
                      const std::vector<uint8_t>& comp2byte_): bwt(bwt_),
                                                               C(C_),
                                                               byte2comp(byte2comp_),
                                                               comp2byte(comp2byte_){}

    [[nodiscard]] size_t size() const {
        return bwt.size();
    }

    [[nodiscard]] size_t alphabet() const {
        return C.size();
    }

    [[nodiscard]] std::pair<uint8_t, size_t> lf(size_t idx) const {
        auto res = bwt.inverse_select(idx);
        size_t next = C[byte2comp[res.second]] + res.first;
        return {res.second, next};
    }

    [[nodiscard]] inline uint8_t get_dummy() const {
        return dummy;
    }

    [[nodiscard]] inline size_t tot_strings() const {
        return C[1]-C[0];
    }

    [[nodiscard]] inline std::pair<uint64_t, uint64_t> backward_search(const std::string &pat) const {
        size_t l=0, r=bwt.size()-1, j=pat.size();
        uint8_t cc;
        while(j-->0 && l<=r){
            cc = byte2comp[uint8_t(pat[j])];
            l = C[cc] + bwt.rank(l, pat[j]); // count c in bwt[0..l-1]
            r = C[cc] + bwt.rank(r+1, pat[j]) - 1; // count c in bwt[0..r]
        }
        return {l, r};
    }

    [[nodiscard]] inline std::tuple<uint64_t, uint64_t, uint64_t> count_with_head(const std::string &pat) const {
        size_t l=0, r=bwt.size()-1, j=pat.size();
        uint8_t cc;
        std::pair<uint64_t, uint64_t> head[2]={{0,0}, {j-1, l}};
        while(j-->0 && l<=r){
            cc = byte2comp[uint8_t(pat[j])];
            auto res = bwt.template rank<true>(l, pat[j]);
            head[res.second] = {j, l};
            //std::cout<<head[1].first<<" "<<head[1].second<<std::endl;

            l = C[cc] + res.first; // count c in bwt[0..l-1]
            r = C[cc] + bwt.rank(r+1, pat[j]) - 1; // count c in bwt[0..r]
        }
        //std::cout<<"A:"<<pat<<" / "<<head[1].first<<" "<<head[1].second<<std::endl;
        return {l, r, head[1].second};
    }
};
#endif //TEST_RL_BCR_BWT_FM_INDEX_H
