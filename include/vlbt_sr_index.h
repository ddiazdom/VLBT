//
// Created by Diaz, Diego on 30.5.2025.
//

#ifndef VLBT_SR_INDEX_H
#define VLBT_SR_INDEX_H

#include "utils.h"

template<class bwt_dt_type, class phi_type>
struct vlbt_sr_index{

    typedef bwt_dt_type bwt_t;
    typedef phi_type phi_t;

    bwt_dt_type bwt_with_th;
    phi_type phi;
    uint32_t ssamp_val;

    [[nodiscard]] inline std::pair<uint64_t, uint64_t> backward_search(const std::string &pat) const {
        size_t l=0, r=bwt_with_th.size()-1, j=pat.size();
        /*uint8_t cc;
        while(j-->0 && l<=r){
            cc = byte2comp[pat[j]];
            l = C[cc] + bwt_with_th.rank(l, pat[j]); // count c in bwt[0..l-1]
            r = C[cc] + bwt_with_th.rank(r+1, pat[j]) - 1; // count c in bwt[0..r]
        }*/
        return {l, r};
    }

    inline uint64_t locate(std::string& pattern) const {

    }

    [[nodiscard]] inline uint64_t size() const {
        return bwt_with_th.size();
    }

    [[nodiscard]] inline int64_t rank(size_t i, uint8_t symbol) const {
        return bwt_with_th.rank(i, symbol);
    }

    [[nodiscard]] inline uint8_t eff2byte(uint8_t c) const {
        return bwt_with_th.eff2byte(c);
    }

    const std::vector<uint8_t>& get_packed_alpha(){
        return bwt_with_th.packed_alpha;
    }

    const std::vector<uint8_t>& get_unpacked_alpha(){
        return bwt_with_th.unpacked_alpha;
    }

    size_t serialize(std::ostream & ofs) const {
        size_t written_bytes = 0;
        written_bytes += serialize_elm(ofs, ssamp_val);
        written_bytes+=bwt_with_th.serialize(ofs);
        written_bytes+=phi.serialize(ofs);
        return written_bytes;
    }

    void load(std::istream & ifs){
        load_elm(ifs, ssamp_val);
        bwt_with_th.load(ifs);
        phi.load(ifs);
    }
};

#endif //VLBT_SR_INDEX_H
