/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */

#ifndef VLBT_SR_INDEX_H
#define VLBT_SR_INDEX_H

#include "vlbt_bwt.h"
#include "vlbt_phi.h"

template<class bwt_dt_type, class phi_type>
struct vlbt_sr_index{

    static_assert(bwt_dt_type::tag==RLBWT_WITH_TOEHOLDS);
    static constexpr VLBT_TYPE tag = SRI_VALID_AREA;

    typedef bwt_dt_type bwt_t;
    typedef phi_type phi_t;

    bwt_dt_type bwt;//bwt with toeholds
    phi_type phi;//phi function

    [[nodiscard]] std::pair<uint64_t, uint64_t> count(const std::string &pat) const {
        return bwt.count(pat);
    }

    [[nodiscard]] std::tuple<uint64_t, uint64_t, uint64_t> count_with_head(const std::string &pat) const {
        return bwt.count_with_head(pat);
    }

    [[nodiscard]] uint8_t operator[](const size_t idx) const {
        return bwt[idx];
    }

    [[nodiscard]] std::vector<uint64_t> locate(const std::string& pattern) {
        int64_t l, r, sa_val;

        //obtain (l,r) for the range SA{i..r-1] and the value of sa_val=SA[i]
        std::tie(l, r, sa_val)  = bwt.count_with_head(pattern);

        //no occurrences
        if(l>r) std::vector<uint64_t>();

        //compute the occurrences in SA[i+1..r-1]
        const int64_t n_occ = r-l+1;
        std::vector<uint64_t> occ(n_occ);
        occ[0] = sa_val;
        if constexpr (phi_type::variant == NO_VALID_AREA) {
            for (int64_t i = 1; i<n_occ; i++) {
                sa_val = bwt.decode_sa(l+i);
                if (sa_val < 0) {
                    sa_val = phi(occ[i-1]);
                }
                assert(sa_val>=0);
                occ[i] = sa_val;
            }
        } else if constexpr (phi_type::variant == WITH_VALID_AREA) {
            for (int64_t i = 1; i<n_occ; i++) {
                sa_val = phi(occ[i-1]);
                if(sa_val<0) {//invalid, resort to lf mapping to find the sa value
                    sa_val = bwt.decode_sa(l+i);
                }
                assert(sa_val>=0);
                occ[i] = sa_val;
            }
        } else if constexpr (phi_type::variant==NO_SUBSAMPLING) {//standard r-index
            for (int64_t i = 1; i<n_occ; i++) {
                occ[i] = phi(occ[i-1]);
            }
        }
        return occ;
    }

    [[nodiscard]] uint64_t size() const {
        return bwt.size();
    }

    [[nodiscard]] uint64_t subsampling_value() const {
        return bwt.subsampling_value();
    }

    [[nodiscard]] uint64_t tot_runs() const {
        return bwt.orig_runs;
    }

    [[nodiscard]] std::pair<uint64_t, uint8_t> inverse_select(size_t i) const {
        return bwt.inverse_select(i);
    }

    [[nodiscard]] int64_t rank(size_t i, uint8_t symbol) const {
        return bwt.rank(i, symbol);
    }

    [[nodiscard]] uint8_t eff2byte(uint8_t c) const {
        return bwt.eff2byte(c);
    }

    const std::vector<uint8_t>& get_packed_alpha(){
        return bwt.packed_alpha;
    }

    size_t alphabet_size() const {
        return bwt.sigma;
    }

    const std::vector<uint8_t>& get_unpacked_alpha(){
        return bwt.unpacked_alpha;
    }

    size_t serialize(std::ostream & ofs) const {

        size_t written_bytes = 0;

        //this is to make sure the data structure is loaded
        //with the correct template arguments
        written_bytes+= serialize_elm(ofs, tag);
        size_t b_size_phi = phi_t::block_size;
        written_bytes+= serialize_elm(ofs, b_size_phi);
        size_t b_size_bwt = bwt_t::block_size;
        written_bytes+= serialize_elm(ofs, b_size_bwt);
        //

        written_bytes+= bwt.serialize(ofs);
        written_bytes+= phi.serialize(ofs);
        return written_bytes;
    }

    void load(std::istream & ifs){

        VLBT_TYPE tmp_tag;
        load_elm(ifs, tmp_tag);
        assert(tmp_tag==tag);
        size_t b_size_phi, b_size_bwt;
        load_elm(ifs, b_size_phi);
        assert(phi_t::block_size==b_size_phi);
        load_elm(ifs, b_size_bwt);
        assert(bwt_t::block_size==b_size_bwt);

        bwt.load(ifs);
        phi.load(ifs);
    }
};

template<size_t b_size_bwt, size_t b_size_phi>
using vlbt_sri_va = vlbt_sr_index<
    vlbt_bwt<RLBWT_WITH_TOEHOLDS, b_size_bwt>,
    vlbt_phi<WITH_VALID_AREA, b_size_phi>>;
#endif //VLBT_SR_INDEX_H
