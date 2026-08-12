//
// Created by Diaz, Diego on 11.8.2026.
//

#ifndef VLBT_BWT_STREAMS_H
#define VLBT_BWT_STREAMS_H

#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include "vbyte.h"
#include "memory_handler.hpp"
#include "logger.h"

enum BWT_FORMAT{
    RL_VBYTE=0,
    PLAIN=1
};

template<typename sym_type=uint8_t>
struct plain_bwt_reader_t {

private:

    sym_type last_sym=0;
    uint64_t last_len=0;
    std::ifstream owned{};
    std::istream* in=nullptr;
    sym_type *buffer=nullptr;
    size_t read_elms = 0;
    size_t buff_pos = 0;
    size_t read_runs = 0;

    size_t buff_size= 1024*1024/sizeof(sym_type);
    bool finished = false;
    bool is_closed = false;

    bool refill() {
        in->read(reinterpret_cast<char *>(buffer), buff_size*sizeof(sym_type));
        const size_t read_bytes = in->gcount();
        read_elms = read_bytes/sizeof(sym_type);
        if (read_elms*sizeof(sym_type)!=read_bytes) {
            LOG_INFO("malformed BWT");
            exit(1);
        }
        buff_pos=0;
        finished = read_bytes==0;
        return !finished;
    }

public:
    plain_bwt_reader_t(const std::string &file, size_t bf_size=1024*1024) {
        buff_size = bf_size/sizeof(sym_type);
        if(buff_size==0){
            LOG_INFO("The read buffer must hold at least one symbol");
            exit(1);
        }
        owned.open(file, std::ifstream::binary);
        if (!owned) {
            LOG_INFO("Cannot open file "+file);
            exit(1);
        }
        in = &owned;
        buffer = mem::allocate<sym_type>(buff_size);
        refill();
    }

    plain_bwt_reader_t(std::istream& s, size_t bf_size=1024*1024): in(&s) {
        buff_size = bf_size/sizeof(sym_type);
        if(buff_size==0){
            LOG_INFO("The read buffer must hold at least one symbol");
            exit(1);
        }
        buffer = mem::allocate<sym_type>(buff_size);
        refill();
    }

    plain_bwt_reader_t(const plain_bwt_reader_t&) = delete;
    plain_bwt_reader_t& operator=(const plain_bwt_reader_t&) = delete;

    bool next_run() {
        if (finished) return false;
        last_sym=buffer[buff_pos];
        last_len=1;
        read_runs++;
        while (true){
            ++buff_pos;
            if (buff_pos==read_elms && !refill()) break;
            if (last_sym!=buffer[buff_pos]) break;
            ++last_len;
        }
        return true;
    }

    bool next_run(sym_type& sym, uint64_t& len) {
        const bool res = next_run();
        if (res) {
            sym = last_sym;
            len = last_len;
        }
        return res;
    }

    bool read_run(sym_type& sym, uint64_t&len) const {
        sym = last_sym;
        len = last_len;
        return !finished;
    }

    [[nodiscard]] size_t n_read_runs() const {
        return read_runs;
    }

    void close() {
        if(!is_closed) {
            owned.close();
            mem::deallocate(buffer);
        }
        is_closed = true;
    }

    ~plain_bwt_reader_t() {
        close();
    }
};

struct rle_vbyte_bwt_reader {

private:

    static constexpr uint64_t exhausted_sentinel = std::numeric_limits<uint64_t>::max();
    size_t buff_cap=1024*1024;
    size_t buff_len=0;
    size_t buff_pos=0;
    size_t tot_read_bytes=0;
    size_t last_vbyte_pos=0;
    size_t read_runs=0;
    uint64_t last_sym=0;
    uint64_t last_len=0;
    uint8_t *buffer=nullptr;
    bool closed=false;
    std::ifstream owned{};
    std::istream * in =nullptr;

    void update_last_vbyte_pos() {
        assert(buffer!=nullptr);
        assert(buff_len>0);
        last_vbyte_pos = buff_len-1;
        while(last_vbyte_pos>0 && buffer[last_vbyte_pos]<128) last_vbyte_pos--;
        if(buffer[last_vbyte_pos]<128) {
            LOG_ERROR("Malformed vbyte stream");
            exit(1);
        }
    }

    //a complete vbyte is available when the buffer holds one that terminates inside it.
    //the buff_pos<buff_len test also covers the empty buffer, where last_vbyte_pos is stale
    [[nodiscard]] bool vbyte_available() const {
        return buff_pos<buff_len && buff_pos<=last_vbyte_pos;
    }

    bool ensure_vbyte() {
        if(vbyte_available()) return true;
        refill();
        return vbyte_available();
    }

    //precondition: vbyte_available()
    void decode_elm(uint64_t& elm) {
        size_t read_vbyte = vbyte::read(buffer+buff_pos, elm);
        buff_pos+=read_vbyte;
        tot_read_bytes+=read_vbyte;
    }

    void refill() {
        const size_t rem_bytes = buff_len - buff_pos;
        memmove(buffer, buffer+buff_pos, rem_bytes);
        in->read(reinterpret_cast<char *>(buffer+rem_bytes), buff_cap-rem_bytes);
        const size_t read_bytes = in->gcount();
        buff_len = read_bytes+rem_bytes;
        buff_pos = 0;
        //the tail guard has to follow the data: after a short read the bytes
        //in [buff_len, buff_cap) are stale, not 0x80
        memset(buffer+buff_len, 0x80, 8);
        if (buff_len>0) {
            update_last_vbyte_pos();
        }
    }

public:
    explicit rle_vbyte_bwt_reader(const std::string& file, size_t bf_cap=1024*1024) {
        if(bf_cap<8){//a vbyte code is at most 8 bytes, so a smaller buffer can never assemble one
            LOG_INFO("The read buffer must be at least 8 bytes");
            exit(1);
        }
        buff_cap = bf_cap;
        owned.open(file, std::ios::binary);
        if(!owned) {
            LOG_INFO("Cannot open file "+file);
            exit(1);
        }
        in = &owned;
        buffer = mem::allocate<uint8_t>(buff_cap+8);//the +8 is to decode vbyte fast
        memset(buffer + buff_cap, 0x80, 8);//fill the tail with termination (0x80) vbyte symbols
        refill();
    }

    explicit rle_vbyte_bwt_reader(std::istream& s, size_t bf_cap=1024*1024): in(&s) {
        if(bf_cap<8){//a vbyte code is at most 8 bytes, so a smaller buffer can never assemble one
            LOG_INFO("The read buffer must be at least 8 bytes");
            exit(1);
        }
        buff_cap = bf_cap;
        buffer = mem::allocate<uint8_t>(buff_cap+8);//the +8 is to decode vbyte fast
        memset(buffer + buff_cap, 0x80, 8);//fill the tail with termination (0x80) vbyte symbols
        refill();
    }

    void rewind() {
        in->clear();
        in->seekg(0);
        if(!*in){//pipes are not seekable: the caller must materialise the data first
            LOG_INFO("This BWT stream cannot be rewound");
            exit(1);
        }
        buff_pos=0;
        buff_len=0;
        last_vbyte_pos=0;
        read_runs=0;
        last_sym=0;
        last_len=0;
        tot_read_bytes=0;
        refill();
    }

    bool next_run() {
        //the end of the input is only legal at a run boundary
        if (!ensure_vbyte()) {
            last_sym = exhausted_sentinel;
            last_len = exhausted_sentinel;
            return false;
        }
        decode_elm(last_sym);
        if (!ensure_vbyte()) {
            LOG_INFO("Truncated vbyte BWT (run symbol without length)");
            exit(1);
        }
        decode_elm(last_len);
        read_runs++;
        return true;
    }

    bool next_run(uint64_t& sym, uint64_t&len) {
        next_run();
        sym = last_sym;
        len = last_len;
        return sym!=exhausted_sentinel;
    }

    bool read_run(uint64_t& sym, uint64_t&len) const {
        sym = last_sym;
        len = last_len;
        return sym!=exhausted_sentinel;
    }

    [[nodiscard]] size_t n_read_runs() const {
        return read_runs;
    }

    rle_vbyte_bwt_reader(const rle_vbyte_bwt_reader&) = delete;
    rle_vbyte_bwt_reader& operator=(const rle_vbyte_bwt_reader&) = delete;

    void close() {
        mem::deallocate(buffer);
        owned.close();
        closed = true;
    }

    ~rle_vbyte_bwt_reader() noexcept {
        if (!closed) {
            close();
        }
    }
};

struct vbyte_rle_bwt_writer {

private:
    size_t last_sym=0;
    size_t last_len=0;
    size_t n_runs=0;
    size_t buff_len=0;
    size_t buff_pos=0;
    size_t seq_size=0;
    uint8_t *buffer=nullptr;
    bool closed=false;
    std::ofstream ofs{};

    void compress_previous() {
        if (n_runs==0) return;

        const size_t vlen_sym = vbyte_len(last_sym);
        const size_t vlen_len = vbyte_len(last_len);
        assert(vlen_sym + vlen_len <= buff_len);

        //store to file if we exceed buffer size
        if((buff_pos+vlen_sym+vlen_len)>=buff_len) {
            ofs.write(reinterpret_cast<char *>(buffer),
                      static_cast<std::streamsize>(buff_pos));
            if (!ofs) {
                LOG_INFO("Cannot write to file");
                exit(1);
            }
            buff_pos = 0;
        }

        //write the vbytes and move the offset
        buff_pos+=vbyte::write(buffer+buff_pos, last_sym);
        buff_pos+=vbyte::write(buffer+buff_pos, last_len);
    }

public:

    void push_back(const size_t sym_, const size_t len_) {

        if (n_runs > 0 && sym_ == last_sym) {//only update the length if the symbol is the same
            last_len += len_;
            seq_size += len_;
            return;
        }

        compress_previous();
        last_sym = sym_;
        last_len = len_;
        seq_size+=len_;
        n_runs++;
    }

    [[nodiscard]] size_t sym() const {
        return last_sym;
    }

    [[nodiscard]] size_t len() const {
        return last_len;
    }

    [[nodiscard]] size_t size() const {
        return seq_size;
    }

    [[nodiscard]] size_t tot_runs() const {
        return n_runs;
    }

    [[nodiscard]] double avg_len() const {
        if (n_runs == 0) return 0;
        return static_cast<double>(seq_size)/static_cast<double>(n_runs);
    }

    void update_sym(size_t sym_) {
        last_sym = sym_;
    }

    void increase_len(size_t len_) {
        last_len += len_;
        seq_size += len_;
    }

    void decrease_len(size_t len_) {
        assert(len_<last_len);
        last_len -= len_;
        seq_size -= len_;
    }

    explicit vbyte_rle_bwt_writer(const std::string& output_file, size_t buff_bytes= 1024 * 1024) {
        ofs.open(output_file,  std::ios::binary);
        if (!ofs) {
            LOG_INFO("Cannot open file "+output_file);
            exit(1);
        }
        buff_len = std::max<size_t>(16, buff_bytes);
        buffer = mem::allocate<uint8_t>(buff_len);
    }

    vbyte_rle_bwt_writer(const vbyte_rle_bwt_writer&) = delete;

    vbyte_rle_bwt_writer& operator=(const vbyte_rle_bwt_writer&) = delete;

    void close() {
        if (!closed) {
            compress_previous();
            if(buff_pos>0) {
                ofs.write(reinterpret_cast<char *>(buffer), static_cast<std::streamsize>(buff_pos));
                buff_pos=0;
                ofs.flush();
            }
            if (!ofs) {
                LOG_INFO("Cannot write to file");
                exit(1);
            }
            closed = true;
        }
    }

    ~vbyte_rle_bwt_writer() noexcept {
        close();
        buff_len = 0;
        ofs.close();
        mem::deallocate(buffer);
    }
};

inline void plain2vbyte_rle(const std::string& plain_bwt, const std::string& vbyte_rle_bwt){

    vbyte_rle_bwt_writer bwt_out(vbyte_rle_bwt);
    std::ifstream owned;
    if(plain_bwt!="-"){
        owned.open(plain_bwt, std::ios::binary);
        if(!owned) {
            LOG_INFO("Cannot open "+plain_bwt);
            exit(1);
        }
    }
    plain_bwt_reader_t bwt_in(plain_bwt=="-" ? static_cast<std::istream&>(std::cin) : owned);

    uint64_t len;
    uint8_t sym;
    while (bwt_in.next_run(sym, len)) {
        bwt_out.push_back(sym, len);
    }
    bwt_in.close();
    bwt_out.close();

    if(bwt_out.tot_runs()==0){
        LOG_INFO("The BWT file "+plain_bwt+" is empty");
        exit(1);
    }
}
#endif //VLBT_BWT_STREAMS_H