//
// Created by Diaz, Diego on 30.5.2025.
//

#ifndef VLBT_SR_INDEX_H
#define VLBT_SR_INDEX_H

template<class bwt_dt_type, class phi_type>
struct vlbt_sr_index{

    static_assert(bwt_dt_type::variant==WITH_TOEHOLDS);

    typedef bwt_dt_type bwt_t;
    typedef phi_type phi_t;

    bwt_dt_type bwt;//bwt with toeholds
    phi_type phi;//phi function

    [[nodiscard]] inline std::pair<uint64_t, uint64_t> count(const std::string &pat) const {
        return bwt.count(pat);
    }

    [[nodiscard]] inline std::tuple<uint64_t, uint64_t, uint64_t> count_with_head(const std::string &pat) const {
        return bwt.count_with_head(pat);
    }

    [[nodiscard]] inline uint8_t operator[](const size_t idx) const {
        return bwt[idx];
    }

    [[nodiscard]] inline std::vector<uint64_t> locate(const std::string& pattern) const {
        int64_t l, r, sa_val;

        //obtain (l,r) for the range SA{i..r-1] and the value of sa_val=SA[i]
        std::tie(l, r, sa_val)  = bwt.count_with_head(pattern);

        //no occurrences
        if(l>r) std::vector<uint64_t>();

        //compute the occurrences in SA[i+1..r-1]
        const size_t len = r-l+1;
        std::vector<uint64_t> occ(len);
        occ[0] = sa_val;
        if constexpr (phi_type::variant == NO_VALID_AREA) {
            for (int64_t i = 1; i<len; i++) {
                sa_val = bwt.decode_sa(l+i);
                if (sa_val < 0) {
                    sa_val = phi(occ[i-1]);
                }
                assert(sa_val>=0);
                occ[i] = sa_val;
            }
        } else if constexpr (phi_type::variant == WITH_VALID_AREA) {
            for (int64_t i = 1; i<len; i++) {
                sa_val = phi(occ[i-1]);
                if(sa_val<0) {//invalid, resort to lf mapping to find the sa value
                    sa_val = bwt.decode_sa(l+i);
                }
                assert(sa_val>=0);
                occ[i] = sa_val;
            }
        } else if constexpr (phi_type::variant==NO_SUBSAMPLING) {//standard r-index
            for (int64_t i = 1; i<len; i++) {
                occ[i] = phi(occ[i-1]);
            }
        }
        return occ;
    }

    [[nodiscard]] inline uint64_t size() const {
        return bwt.size();
    }

    [[nodiscard]] inline uint64_t subsampling_value() const {
        return bwt.subsampling_value();
    }

    [[nodiscard]] inline uint64_t tot_runs() const {
        return bwt.orig_runs;
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
        written_bytes+= bwt.serialize(ofs);
        written_bytes+= phi.serialize(ofs);
        return written_bytes;
    }

    void load(std::istream & ifs){
        bwt.load(ifs);
        phi.load(ifs);
    }
};
#endif //VLBT_SR_INDEX_H
