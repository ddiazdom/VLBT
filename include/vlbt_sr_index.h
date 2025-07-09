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

    bwt_dt_type bwt;//bwt with toeholds
    phi_type phi;//phi function
    uint32_t ssamp_val=0;//value of subsampling

    [[nodiscard]] inline std::pair<uint64_t, uint64_t> count(const std::string &pat) const {
        return bwt.count(pat);
    }

    inline auto locate(std::string& pattern) const {
        return bwt.count_with_head(pattern);
    }

    [[nodiscard]] inline uint64_t size() const {
        return bwt.size();
    }

    [[nodiscard]] inline std::pair<uint64_t, uint8_t> inverse_select(size_t i) const {
        return bwt.inverse_select(i);
    }


    [[nodiscard]] inline int64_t rank(size_t i, uint8_t symbol) const {
        return bwt.rank(i, symbol);
    }

    [[nodiscard]] inline uint8_t eff2byte(uint8_t c) const {
        return bwt.eff2byte(c);
    }

    const std::vector<uint8_t>& get_packed_alpha(){
        return bwt.packed_alpha;
    }

    const std::vector<uint8_t>& get_unpacked_alpha(){
        return bwt.unpacked_alpha;
    }

    size_t serialize(std::ostream & ofs) const {
        size_t written_bytes = 0;
        written_bytes+= serialize_elm(ofs, ssamp_val);
        written_bytes+= bwt.serialize(ofs);
        written_bytes+= phi.serialize(ofs);
        return written_bytes;
    }

    void load(std::istream & ifs){
        load_elm(ifs, ssamp_val);
        bwt.load(ifs);
        phi.load(ifs);
    }
};
#endif //VLBT_SR_INDEX_H
