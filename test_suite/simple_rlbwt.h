//
// Created by Diaz, Diego on 3.2.2022.
//

#ifndef TEST_SIMPLE_RLBWT_H
#define TEST_SIMPLE_RLBWT_H
#include <sdsl/construct.hpp>

template<class bwt_type>
struct simple_rlbwt{

    bwt_type bwt;
    std::vector<uint64_t> C;
    std::vector<uint8_t> byte2comp;
    uint8_t terminator{};
    typedef typename bwt_type::size_type size_type;

    simple_rlbwt()= default;

    explicit simple_rlbwt(bwt_type& bwt_): bwt(bwt_) {

        byte2comp.resize(256);

        for(size_t sym=0;sym<256;sym++){
            size_t sym_rank = bwt.rank(bwt.size(), sym);
            if(sym_rank>0) {
                if(C.empty()) terminator = sym;
                byte2comp[sym] = C.size();
                C.push_back(sym_rank);
            }
        }

        size_t acc=0;
        for(size_t i=0;i<C.size();i++){
            const size_t tmp = C[i];
            C[i] = acc;
            acc+=tmp;
        }
        C.push_back(acc);
    }


    template<class T>
    static uint64_t serialize_vec(const std::vector<T>& vec, std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") {
        size_t sz = vec.size();
        out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
        out.write(reinterpret_cast<const char *>(vec.data()), sizeof(T)*sz);
        return sizeof(sz) + sizeof(T)*sz;
    }

    template<class T>
    static void load_vec(std::vector<T>& vec, std::istream& in) {
        size_t sz;
        in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
        vec.resize(sz);
        in.read(reinterpret_cast<char *>(vec.data()), sizeof(T)*sz);
    }

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        std::uint64_t written_bytes = 0;
        written_bytes += sdsl::serialize(bwt, out, child, "bwt");
        written_bytes += serialize_vec(C, out, child, "C");
        written_bytes += serialize_vec(byte2comp, out, child, "byte2comp");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    // Load the data structure from a stream and set the supported vector
    void load(std::istream& in) {
        sdsl::load(bwt, in);
        load_vec(C, in);
        load_vec(byte2comp, in);
    }

    [[nodiscard]] size_t size() const {
        return bwt.size();
    }

    [[nodiscard]] size_t alphabet() const {
        return C.size()-1;
    }

    [[nodiscard]] std::pair<uint8_t, size_t> lf(size_t idx) const {
        auto res = bwt.inverse_select(idx);
        size_t next = C[byte2comp[res.second]] + res.first;
        return {res.second, next};
    }

    [[nodiscard]] uint8_t get_terminator_sym() const {
        return terminator;
    }

    [[nodiscard]] size_t tot_strings() const {
        return C[1]-C[0];
    }

    auto rank(size_t pos, uint8_t sym) const {
        return bwt.rank(pos, sym);
    }

    auto operator[](size_type idx) const {
        return bwt[idx];
    }

    auto inverse_select(size_type idx) const {
        return bwt.inverse_select(idx);
    }

    [[nodiscard]] std::pair<uint64_t, uint64_t> count(const std::string &pat) const {
        size_t l=0, r=bwt.size()-1, j=pat.size();
        //std::cout<<l<<" "<<r<<std::endl;
        while(j-->0 && l<=r){
            const uint8_t cc = byte2comp[static_cast<uint8_t>(pat[j])];
            const size_t lb = bwt.rank(l, pat[j]);
            const size_t rb = bwt.rank(r+1, pat[j]);
            l = C[cc] + lb; // count c in bwt[0..l-1]
            r = C[cc] +  rb - 1; // count c in bwt[0..r]
        }
        return {l, r};
    }

    [[nodiscard]] uint8_t eff2byte(const size_t cmp_sym) const {
        uint8_t byte_sym=0;
        while(cmp_sym!=byte2comp[byte_sym]) byte_sym++;
        return byte_sym;
    }
};
#endif //TEST_RL_BCR_BWT_FM_INDEX_H
