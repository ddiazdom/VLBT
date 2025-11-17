//
// Created by diego on 02-09-20.
//

#ifndef LPG_COMPRESSOR_BITSTREAM_H
#define LPG_COMPRESSOR_BITSTREAM_H

#include <iostream>
#include <limits>
#include <cstring>
#include <cassert>
#ifdef __linux__
#include <malloc.h>
#endif
#ifdef __BMI2__
#include <immintrin.h>
#endif


#include "utils.h"

template <class T> struct mem {

    static T * allocate(const size_t n) {
        if (n == 0) {return nullptr;}
        if (n > static_cast<size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        void * const pv = malloc(n * sizeof(T));
        if (!pv) { throw std::bad_alloc(); }
        return static_cast<T *>(pv);
    }

    static T * reallocate(T* old_ptr, const size_t n) {
        if (n == 0) {return nullptr;}
        if (n > static_cast<size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        void * const pv = realloc(old_ptr, n * sizeof(T));
        if (!pv) { throw std::bad_alloc(); }
        return static_cast<T *>(pv);
    }

    static void deallocate(T * const p) {
        free(p);
#ifdef __linux__
        malloc_trim(0);
#endif
    }
};

const uint64_t ps_overflow[] = {
        0x8080808080808080ULL,
        0x7f7f7f7f7f7f7f7fULL,
        0x7e7e7e7e7e7e7e7eULL,
        0x7d7d7d7d7d7d7d7dULL,
        0x7c7c7c7c7c7c7c7cULL,
        0x7b7b7b7b7b7b7b7bULL,
        0x7a7a7a7a7a7a7a7aULL,
        0x7979797979797979ULL,
        0x7878787878787878ULL,
        0x7777777777777777ULL,
        0x7676767676767676ULL,
        0x7575757575757575ULL,
        0x7474747474747474ULL,
        0x7373737373737373ULL,
        0x7272727272727272ULL,
        0x7171717171717171ULL,
        0x7070707070707070ULL,
        0x6f6f6f6f6f6f6f6fULL,
        0x6e6e6e6e6e6e6e6eULL,
        0x6d6d6d6d6d6d6d6dULL,
        0x6c6c6c6c6c6c6c6cULL,
        0x6b6b6b6b6b6b6b6bULL,
        0x6a6a6a6a6a6a6a6aULL,
        0x6969696969696969ULL,
        0x6868686868686868ULL,
        0x6767676767676767ULL,
        0x6666666666666666ULL,
        0x6565656565656565ULL,
        0x6464646464646464ULL,
        0x6363636363636363ULL,
        0x6262626262626262ULL,
        0x6161616161616161ULL,
        0x6060606060606060ULL,
        0x5f5f5f5f5f5f5f5fULL,
        0x5e5e5e5e5e5e5e5eULL,
        0x5d5d5d5d5d5d5d5dULL,
        0x5c5c5c5c5c5c5c5cULL,
        0x5b5b5b5b5b5b5b5bULL,
        0x5a5a5a5a5a5a5a5aULL,
        0x5959595959595959ULL,
        0x5858585858585858ULL,
        0x5757575757575757ULL,
        0x5656565656565656ULL,
        0x5555555555555555ULL,
        0x5454545454545454ULL,
        0x5353535353535353ULL,
        0x5252525252525252ULL,
        0x5151515151515151ULL,
        0x5050505050505050ULL,
        0x4f4f4f4f4f4f4f4fULL,
        0x4e4e4e4e4e4e4e4eULL,
        0x4d4d4d4d4d4d4d4dULL,
        0x4c4c4c4c4c4c4c4cULL,
        0x4b4b4b4b4b4b4b4bULL,
        0x4a4a4a4a4a4a4a4aULL,
        0x4949494949494949ULL,
        0x4848484848484848ULL,
        0x4747474747474747ULL,
        0x4646464646464646ULL,
        0x4545454545454545ULL,
        0x4444444444444444ULL,
        0x4343434343434343ULL,
        0x4242424242424242ULL,
        0x4141414141414141ULL,
        0x4040404040404040ULL
};

const uint8_t lt_sel[] = {
        0,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
        4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,

        0,0,0,1,0,2,2,1,0,3,3,1,3,2,2,1,
        0,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
        0,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
        5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
        0,6,6,1,6,2,2,1,6,3,3,1,3,2,2,1,
        6,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
        6,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
        5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
        0,7,7,1,7,2,2,1,7,3,3,1,3,2,2,1,
        7,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
        7,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
        5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
        7,6,6,1,6,2,2,1,6,3,3,1,3,2,2,1,
        6,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
        6,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
        5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,

        0,0,0,0,0,0,0,2,0,0,0,3,0,3,3,2,
        0,0,0,4,0,4,4,2,0,4,4,3,4,3,3,2,
        0,0,0,5,0,5,5,2,0,5,5,3,5,3,3,2,
        0,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
        0,0,0,6,0,6,6,2,0,6,6,3,6,3,3,2,
        0,6,6,4,6,4,4,2,6,4,4,3,4,3,3,2,
        0,6,6,5,6,5,5,2,6,5,5,3,5,3,3,2,
        6,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
        0,0,0,7,0,7,7,2,0,7,7,3,7,3,3,2,
        0,7,7,4,7,4,4,2,7,4,4,3,4,3,3,2,
        0,7,7,5,7,5,5,2,7,5,5,3,5,3,3,2,
        7,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
        0,7,7,6,7,6,6,2,7,6,6,3,6,3,3,2,
        7,6,6,4,6,4,4,2,6,4,4,3,4,3,3,2,
        7,6,6,5,6,5,5,2,6,5,5,3,5,3,3,2,
        6,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,

        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
        0,0,0,0,0,0,0,4,0,0,0,4,0,4,4,3,
        0,0,0,0,0,0,0,5,0,0,0,5,0,5,5,3,
        0,0,0,5,0,5,5,4,0,5,5,4,5,4,4,3,
        0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,3,
        0,0,0,6,0,6,6,4,0,6,6,4,6,4,4,3,
        0,0,0,6,0,6,6,5,0,6,6,5,6,5,5,3,
        0,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
        0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,3,
        0,0,0,7,0,7,7,4,0,7,7,4,7,4,4,3,
        0,0,0,7,0,7,7,5,0,7,7,5,7,5,5,3,
        0,7,7,5,7,5,5,4,7,5,5,4,5,4,4,3,
        0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,3,
        0,7,7,6,7,6,6,4,7,6,6,4,6,4,4,3,
        0,7,7,6,7,6,6,5,7,6,6,5,6,5,5,3,
        7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,

        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
        0,0,0,0,0,0,0,5,0,0,0,5,0,5,5,4,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
        0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,4,
        0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,5,
        0,0,0,6,0,6,6,5,0,6,6,5,6,5,5,4,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
        0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,4,
        0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,5,
        0,0,0,7,0,7,7,5,0,7,7,5,7,5,5,4,
        0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
        0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,4,
        0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,5,
        0,7,7,6,7,6,6,5,7,6,6,5,6,5,5,4,

        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
        0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,5,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
        0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,5,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
        0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
        0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
        0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,5,

        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
        0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,

        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7
};

template<class word_t, uint8_t max_dist=std::numeric_limits<word_t>::digits>
struct bit_stream{

    constexpr static uint8_t word_bits = std::numeric_limits<word_t>::digits;
    constexpr static uint8_t word_shift = __builtin_ctz(word_bits);
    const static size_t masks[65];
    using size_type = size_t;

    word_t *stream=nullptr;
    size_t stream_cap=0;//in words

    bit_stream(): stream(nullptr), stream_cap(0){};

    bit_stream(bit_stream&& other) noexcept {
        std::swap(stream, other.stream);
        std::swap(stream_cap, other.stream_cap);
    }

    bit_stream(bit_stream& other) noexcept {
        if(&other!=this && other.stream!=nullptr){
            reserve_in_words(other.stream_cap);
            memcpy(stream, other.stream, words2bytes(other.stream_cap));
        }
    }

    [[nodiscard]] inline size_t capacity_in_bits() const {
        return stream_cap*word_bits;
    }

    [[nodiscard]] inline size_t capacity_in_words() const {
        return stream_cap;
    }

    [[nodiscard]] inline size_t capacity_in_bytes() const {
        return stream_cap*sizeof(word_t);
    }

    inline void reserve_in_bits(size_t bit_size){
        reserve_in_words(bits2words(bit_size));
    }

    inline void reserve_in_bytes(size_t byte_size){
        reserve_in_words(bytes2words(byte_size));
    }

    inline void reserve_in_words(size_t n_words){
        if(n_words>stream_cap){
            if(stream==nullptr){
                //stream = (word_t *)malloc(words2bytes(n_words));
                stream = mem<word_t>::allocate(n_words);
            }else{
                assert(stream_cap!=0);
                //stream = (word_t *)realloc(stream, words2bytes(n_words));
                stream = mem<word_t>::reallocate(stream, n_words);
            }
            stream_cap = n_words;
        }
    }

    [[nodiscard]] static inline size_t bits2words(size_t n_bits) {
        return INT_CEIL(n_bits, (sizeof(word_t)*8));
    }

    [[nodiscard]] static inline size_t bytes2words(size_t n_bytes) {
        return INT_CEIL(n_bytes, sizeof(word_t));
    }

    [[nodiscard]] static inline size_t words2bytes(size_t n_words) {
        return n_words*sizeof(word_t);
    }

    void destroy(){
        if(stream!= nullptr){
            //free(stream);
            mem<word_t>::deallocate(stream);
            stream = nullptr;
        }
        stream_cap=0;
    }

    inline bit_stream& swap(bit_stream& other) {
        std::swap(stream, other.stream);
        std::swap(stream_cap, other.stream_cap);
        return *this;
    }

    inline bit_stream& operator=(bit_stream const& other){
        if(&other!=this){
            reserve_in_words(other.stream_cap);
            memcpy(stream, other.stream, words2bytes(other.stream_cap));
        }
        return *this;
    }

    inline void write(size_t i, size_t j, size_t value){
        size_t cell_i = i >> word_shift;
        size_t i_pos = (i & (word_bits-1UL));

        if constexpr (max_dist==1){
            stream[cell_i] &= ~(1UL << i_pos);
            stream[cell_i] |= value << i_pos;
        }else{
            size_t cell_j = j >> word_shift;
            if(cell_i==cell_j){
                stream[cell_i] &= ~(masks[j-i+1UL] << i_pos);
                stream[cell_i] |= (value & masks[j-i+1UL]) << i_pos;
            }else{
                size_t right = word_bits - i_pos;
                size_t left = 1+(j & (word_bits - 1UL));
                stream[cell_i] = (stream[cell_i] & ~(masks[right] << i_pos)) | (value << i_pos);
                stream[cell_j] = (stream[cell_j] & ~masks[left]) | (value >> right);
            }
        }
    }

    inline void write_chunk(const void* source, size_t i, size_t j){
        size_t tot_bits = j-i+1;
        size_t n_words = INT_CEIL(tot_bits, word_bits);
        size_t left = i & (word_bits - 1UL);
        size_t right = word_bits - left;
        size_t cell_i = i >> word_shift;

        auto tmp_src = reinterpret_cast<const word_t *>(source);

        //TODO use SIMD instructions for this segment
        for(size_t k=0; k < n_words - 1; k++){
            stream[cell_i] = (stream[cell_i] & ~(masks[right] << left)) | (tmp_src[k] << left);
            cell_i++;
            stream[cell_i] = (stream[cell_i] & ~masks[left]) | (tmp_src[k] >> right);
        }
        //

        size_t read_bits = ((n_words - 1) << word_shift);

        write(i + read_bits, j, (tmp_src[n_words-1] & masks[tot_bits-read_bits]));
    }

    [[nodiscard]] inline size_t read(size_t i, size_t j) const{
        if constexpr (max_dist==1){
            return (stream[i>>word_shift] >> (i & (word_bits - 1UL))) & 1UL;
        } else {

            //const __uint128_t combined = (static_cast<__uint128_t>(stream[j>>word_shift]) << word_bits) | stream[i>>word_shift];
            //return combined >> (i & (word_bits-1UL)) & ((1ULL << (j-i+1)) - 1);

            size_t cell_i = i >> word_shift;
            size_t i_pos = i & (word_bits - 1UL);
            size_t cell_j = j >> word_shift;
            if(cell_i == cell_j){
                return (stream[cell_i] >> i_pos) & masks[(j - i + 1UL)];
            }
            size_t right = word_bits-i_pos;
            return ((stream[cell_j] & masks[1+(j & (word_bits - 1UL))]) << right) | ((stream[cell_i] >> i_pos) & masks[right]);

            /*size_t cell_i = i >> word_shift;
            size_t i_pos = (i & (word_bits - 1UL));
            size_t cell_j = j >> word_shift;
            size_t j_pos = (j & (word_bits - 1UL));

            size_t diff_word = cell_i!=cell_j;
            size_t l_end = (j_pos+1)*diff_word;
            size_t r_end = (j-i+1)-(l_end*diff_word);
            return ((stream[cell_j] & masks[l_end])<<r_end) | ((stream[cell_i]>>i_pos) & masks[r_end]);*/
        }
    }

    inline void prefetch(size_t i) const {
        __builtin_prefetch(&stream[i>>word_shift], 0, 0);
    }

    [[nodiscard]] inline bool read_bit(size_t i) const{
        return (stream[i>>word_shift] >> (i & (word_bits - 1UL))) & 1UL;
    }


    [[nodiscard]] inline size_t pop_count(size_t i, size_t j) const {
        size_t cell_i = i >> word_shift;
        size_t i_pos = (i & (word_bits - 1UL));
        size_t cell_j = j >> word_shift;

        if(cell_i != cell_j){
            size_t count=0;
            for(size_t c=cell_i+1;c<cell_j;++c){
                count+=__builtin_popcountll(stream[c]);
            }
            return count + __builtin_popcountll(stream[cell_i] >> i_pos) + __builtin_popcountll(stream[cell_j] & masks[1+(j & (word_bits - 1UL))]);
        }

        return __builtin_popcountll((stream[cell_i] >> i_pos) & masks[(j - i + 1UL)]);
    }

    //from the SDSL
    [[nodiscard]] static inline uint32_t select64_scalar(uint64_t x, uint32_t i) {
        uint64_t s = x, b;  // s = sum
        s = s-((s>>1) & 0x5555555555555555ULL);
        s = (s & 0x3333333333333333ULL) + ((s >> 2) & 0x3333333333333333ULL);
        s = (s + (s >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
        s = 0x0101010101010101ULL*s;
        b = (s+ps_overflow[i]);//&0x8080808080808080ULL;// add something to the partial sums to cause overflow
        i = (i-1)<<8;
        if (b&0x0000000080000000ULL) // byte <=3
            if (b&0x0000000000008000ULL) //byte <= 1
                if (b&0x0000000000000080ULL)
                    return lt_sel[(x&0xFFULL) + i];
                else
                    return 8 +lt_sel[(((x>>8)&0xFFULL)  + i - ((s&0xFFULL)<<8))&0x7FFULL];//byte 1;
            else//byte >1
            if (b&0x0000000000800000ULL) //byte <=2
                return 16+lt_sel[(((x>>16)&0xFFULL) + i - (s&0xFF00ULL))&0x7FFULL];//byte 2;
            else
                return 24+lt_sel[(((x>>24)&0xFFULL) + i - ((s>>8)&0xFF00ULL))&0x7FFULL];//byte 3;
        else//  byte > 3
        if (b&0x0000800000000000ULL) // byte <=5
            if (b&0x0000008000000000ULL) //byte <=4
                return 32+lt_sel[(((x>>32)&0xFFULL) + i - ((s>>16)&0xFF00ULL))&0x7FFULL];//byte 4;
            else
                return 40+lt_sel[(((x>>40)&0xFFULL) + i - ((s>>24)&0xFF00ULL))&0x7FFULL];//byte 5;
        else// byte >5
        if (b&0x0080000000000000ULL) //byte<=6
            return 48+lt_sel[(((x>>48)&0xFFULL) + i - ((s>>32)&0xFF00ULL))&0x7FFULL];//byte 6;
        else
            return 56+lt_sel[(((x>>56)&0xFFULL) + i - ((s>>40)&0xFF00ULL))&0x7FFULL];//byte 7;
        return 0;
    }

    //from the SLDSL
    [[nodiscard]] static inline uint32_t select64(uint64_t x, size_t i) {
#ifdef __BMI2__
        // index i is 1-based here, (i-1) changes it to 0-based
        return __builtin_ctzll(_pdep_u64(1ull << (i-1), x));
#elif defined(__SSE4_2__)
        uint64_t s = x, b;
        s = s-((s>>1) & 0x5555555555555555ULL);
        s = (s & 0x3333333333333333ULL) + ((s >> 2) & 0x3333333333333333ULL);
        s = (s + (s >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
        s = 0x0101010101010101ULL*s;
        // now s contains 8 bytes s[7],...,s[0]; s[j] contains the cumulative sum
        // of (j+1)*8 least significant bits of s
        b = (s+ps_overflow[i]) & 0x8080808080808080ULL;
        // ps_overflow contains a bit mask x consisting of 8 bytes
        // x[7],...,x[0] and x[j] is set to 128-j
        // => a byte b[j] in b is >= 128 if cum sum >= j

        // __builtin_ctzll returns the number of trailing zeros, if b!=0
        int  byte_nr = __builtin_ctzll(b) >> 3;   // byte nr in [0..7]
        s <<= 8;
        i -= (s >> (byte_nr<<3)) & 0xFFULL;
        return (byte_nr << 3) + lt_sel[((i-1) << 8) + ((x>>(byte_nr<<3))&0xFFULL) ];
#else
        return select64_scalar(x, i);
#endif
    }

    [[nodiscard]] inline size_t select(size_t i, size_t j, size_t r) const {
        if((j-i+1)<=64){
            return select64(read(i, j), r);
        }else {
            size_t cell_i = i >> word_shift;
            size_t i_pos = (i & (word_bits - 1UL));
            size_t cell_j = j >> word_shift;

            /*size_t tmp=0;
            for (size_t k = i, m=0; k <= j; k++,m++) {
                tmp+= read_bit(k);
                std::cout<<m<<":"<<read_bit(k);
                if(tmp==r && read_bit(k)){
                    std::cout<<"* ";
                }else{
                    std::cout<<" ";
                }
            }
            std::cout << "" << std::endl;*/

            uint64_t word = stream[cell_i] >> i_pos;
            size_t rank = __builtin_popcountll(word), prev_rank = 0;
            uint32_t sel_answer = 0, consumed = word_bits - i_pos;
            while(rank < r && cell_i<cell_j) {
                sel_answer += consumed;
                prev_rank = rank;
                word = stream[++cell_i];
                rank += __builtin_popcountll(word);
                consumed = 64;
            }
            sel_answer += select64(word, r - prev_rank);
            return sel_answer;
        }
    }

    inline void read_chunk(void* dst, size_t i, size_t j) const{
        size_t tot_bits = j-i+1;
        size_t n_words = INT_CEIL(tot_bits, word_bits);
        size_t left = i & (word_bits - 1UL);
        size_t right = word_bits - left;
        size_t cell_i = i >> word_shift;

        auto tmp_dst = reinterpret_cast<word_t *>(dst);

        //TODO use SIMD instructions for this segment
        for(size_t k=0; k < n_words - 1; k++){
            tmp_dst[k] = (stream[cell_i] >> left) & masks[right];
            tmp_dst[k] |= (stream[++cell_i] & masks[left]) << right;
        }
        //

        size_t read_bits = ((n_words - 1) << word_shift);
        tmp_dst[n_words-1] &= ~masks[tot_bits-read_bits];
        tmp_dst[n_words-1] |= read(i + read_bits, j);
    }

    //compare a segment of the stream with an external source of bits
    inline bool compare_chunk(const void* input, size_t i, size_t bits) const {

        size_t n_words = INT_CEIL(bits, word_bits);
        size_t left = i & (word_bits - 1UL);
        size_t right = word_bits - left;
        size_t cell_i = i >> word_shift;

        auto tmp_in = reinterpret_cast<const word_t *>(input);

        //TODO use SIMD instructions for this segment
        size_t tmp_data;
        for(size_t k=0; k < n_words - 1; k++){
            tmp_data = (stream[cell_i] >> left) & masks[right];
            tmp_data |= (stream[++cell_i] & masks[left]) << right;
            if(tmp_data != tmp_in[k]) return false;
        }
        //
        size_t read_bits = ((n_words - 1) << word_shift);
        return (tmp_in[n_words - 1] & masks[(bits-read_bits)]) == read(i + read_bits, i+bits-1);
    }

    template<class stream_type>
    void copy(size_t n_bits, stream_type& dest_stream){
        size_t n_words = bits2words(n_bits);
        assert(n_words<=stream_cap);
        dest_stream.reserve_in_words(n_words);
        memcpy(dest_stream.stream, stream, words2bytes(n_words));
    }

    //compare a segment of the stream with an external source of bits
    /*inline bool compare_segment(const uint8_t* input, size_t i, size_t bits) const {

        size_t n_words = INT_CEIL(bits, word_bits);
        size_t left = i & (word_bits - 1UL);
        size_t right = word_bits - left;

        for(size_t k=0, cell_i=(i>>word_shift); k < n_words - 1; k+=8, cell_i++){
            const size_t tmp_data = ((stream[cell_i] >> left) & masks[right]) | ((stream[cell_i+1] & masks[left]) << right);
            if(memcmp(input+k, &tmp_data, sizeof(size_t))!=0) return false;
        }

        //size_t cell_i = i >> word_shift;
        size_t read_bits = ((n_words - 1) << word_shift);
        //return (tmp_in[n_words - 1] & masks[(bits-read_bits)]) == read(i + read_bits, i+bits-1);
    }*/

    //compare the segment ]a-bits..a] with the segment ]b-bits..b+bits]
    //return the bit_pos (0-based) of the rightmost different bit (return len if the segments are equal)
    inline size_t inv_com_segments(size_t a, size_t b, size_t& bits) const {
        size_t n_words = INT_CEIL(bits, word_bits);
        size_t rem_bits, data_a, data_b, read_bits=0;

        if(n_words>1){

            size_t rs_a = a & (word_bits - 1UL);
            size_t ls_a = word_bits - rs_a - 1;
            size_t cell_a = a >> word_shift;

            size_t rs_b = b & (word_bits - 1UL);
            size_t ls_b = word_bits - rs_b - 1;
            size_t cell_b = b >> word_shift;

            for(size_t k=n_words; k-->1;){
                //the mask helps to deal with corner cases
                // (if the shift is 64-bits long, for instance)
                data_a = stream[cell_a] << ls_a;
                data_a |= (stream[--cell_a] & (masks[ls_a]<<(rs_a+1UL))) >> (rs_a+1UL);

                data_b = stream[cell_b] << ls_b;
                data_b |= (stream[--cell_b] & (masks[ls_b]<<(rs_b+1UL))) >> (rs_b+1UL);

                if(data_a != data_b){
                    return read_bits + __builtin_clzl(data_a^data_b);
                }
                read_bits+=word_bits;
            }
        }

        rem_bits = bits-read_bits;
        data_a = read(a-bits+1, a-read_bits);
        data_b = read(b-bits+1, b-read_bits);

        data_a^=data_b;
        if(data_a==0){
            return bits;
        }else{
            return read_bits + (rem_bits- ((8*sizeof(unsigned long) - __builtin_clzl(data_a))));
        }
    }

    void concatenate(size_t bytes_a, bit_stream<word_t, max_dist>& stream_b, size_t bytes_b){
        size_t new_size_in_words = INT_CEIL((bytes_a+bytes_b), sizeof(word_t));
        reserve_in_words(new_size_in_words);
        auto * tmp_stream_a = (uint8_t *)stream;
        auto * tmp_stream_b = (uint8_t *)stream_b.stream;
        memcpy(&tmp_stream_a[bytes_a], tmp_stream_b, bytes_b);
    }

    size_type serialize(std::ostream &out) const{
        size_t written_bytes = serialize_elm(out, stream_cap);
        out.write((char *)stream, words2bytes(stream_cap));
        return written_bytes + words2bytes(stream_cap);
    }

    void load(std::istream &in){
        size_t tmp_size;
        load_elm(in, tmp_size);
        reserve_in_words(tmp_size);
        in.read((char *)stream, words2bytes(tmp_size));
    }
};

template<class word_t, uint8_t max_dist>
const size_t bit_stream<word_t, max_dist>::masks[65]={0x0,
                                                      0x1, 0x3, 0x7, 0xF,
                                                      0x1F, 0x3F, 0x7F, 0xFF,
                                                      0x1FF, 0x3FF, 0x7FF, 0xFFF,
                                                      0x1FFF, 0x3FFF, 0x7FFF, 0xFFFF,
                                                      0x1FFFF, 0x3FFFF, 0x7FFFF, 0xFFFFF,
                                                      0x1FFFFF, 0x3FFFFF, 0x7FFFFF, 0xFFFFFF,
                                                      0x1FFFFFF, 0x3FFFFFF, 0x7FFFFFF, 0xFFFFFFF,
                                                      0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF,
                                                      0x1FFFFFFFF, 0x3FFFFFFFF, 0x7FFFFFFFF, 0xFFFFFFFFF,
                                                      0x1FFFFFFFFF, 0x3FFFFFFFFF, 0x7FFFFFFFFF, 0xFFFFFFFFFF,
                                                      0x1FFFFFFFFFF, 0x3FFFFFFFFFF, 0x7FFFFFFFFFF, 0xFFFFFFFFFFF,
                                                      0x1FFFFFFFFFFF, 0x3FFFFFFFFFFF, 0x7FFFFFFFFFFF, 0xFFFFFFFFFFFF,
                                                      0x1FFFFFFFFFFFF, 0x3FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0xFFFFFFFFFFFFF,
                                                      0x1FFFFFFFFFFFFF, 0x3FFFFFFFFFFFFF, 0x7FFFFFFFFFFFFF, 0xFFFFFFFFFFFFFF,
                                                      0x1FFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFF,
                                                      0x1FFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF};

#endif //LPG_COMPRESSOR_BITSTREAM_H
