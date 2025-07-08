//
// Created by Diaz, Diego on 3.2.2022.
//

#ifndef TEST_RL_BCR_BWT_FM_INDEX_H
#define TEST_RL_BCR_BWT_FM_INDEX_H

#include "../include/bwt_io.h"

template<class bwt_type, bool include_sa_samples=false, typename sa_samp_type=uint64_t>
struct fm_index{

    bwt_type& bwt;
    std::vector<uint64_t> C;
    const std::vector<uint8_t>& byte2comp;
    const std::vector<uint8_t>& comp2byte;
    std::vector<uint64_t> sa_samps;
    uint8_t dummy=0;

    explicit fm_index(bwt_type& bwt_,
                      std::vector<uint64_t>& C_,
                      const std::string sa_samples,
                      const std::vector<uint8_t>& byte2comp_,
                      const std::vector<uint8_t>& comp2byte_): bwt(bwt_),
                                                               C(C_),
                                                               byte2comp(byte2comp_),
                                                               comp2byte(comp2byte_){

        if constexpr (include_sa_samples){
            size_t n_samples = std::filesystem::file_size(sa_samples)/sizeof(sa_samp_type);
            std::ifstream sa_subsamp_ifs(sa_samples, std::ios::binary);
            std::vector<sa_samp_type> sa_samp_buffer(n_samples, 0);
            sa_subsamp_ifs.read((char *)sa_samp_buffer.data(), n_samples*sizeof(sa_samp_type));
            assert(n_samples/2==bwt.n_runs());
            sa_samp_buffer.swap(sa_samps);
        }
    }

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
        //std::cout<<l<<" "<<r<<std::endl;
        while(j-->0 && l<=r){
            cc = byte2comp[uint8_t(pat[j])];
            l = C[cc] + bwt.rank(l, pat[j]); // count c in bwt[0..l-1]
            r = C[cc] + bwt.rank(r+1, pat[j]) - 1; // count c in bwt[0..r]
            //std::cout<<l<<" "<<r<<std::endl;
        }
        return {l, r};
    }

    [[nodiscard]] inline int64_t get_sa_samp_of_succ_head(uint64_t i, uint8_t symbol) const {
        int64_t run = bwt.succ_run(i, symbol);
        if(run<0){
            return -1;
        }
        //std::cout<<sa_samps[2*run]<<std::endl;
        //auto lf = bwt.inverse_select(i);
        //i = C[byte2comp[lf.second]] + lf.first;
        //lf = bwt.inverse_select(i);
        //i = C[byte2comp[lf.second]] + lf.first;
        //size_t r = bwt.pos2run(i);
        //std::cout<<sa_samps[2*(r-1)]<<std::endl;
        return sa_samps[2*run];
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
        int64_t sa_samp = get_sa_samp_of_succ_head(head[1].second, pat[head[1].first]);
        //std::cout<<"ground_tr: \""<<pat<<"\" -> "<<l<<" "<<r<<" "<<sa_samp<<std::endl;
        return {l, r, sa_samp};
    }

};
#endif //TEST_RL_BCR_BWT_FM_INDEX_H
