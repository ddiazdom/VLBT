//
// Created by Diaz, Diego on 23.1.2022.
//

#ifndef GBWT_BWT_IO_H
#define GBWT_BWT_IO_H
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <functional>
#include <cstring>
#include <cassert>
#include <cstdint>
#include <filesystem>

class bwt_buff_reader {

    char * buffer =nullptr;
    char * sec_buff=nullptr;

    size_t bpr=0; //bytes per run
    size_t sb=0; // bytes for a run symbol
    size_t fb=0; // bytes for a run freq
    size_t block_bg=0; //start of the current block
    size_t tot_runs =0; //total number of runs in the file
    size_t buffer_size=0; //size of the buffer
    size_t offset =sizeof(size_t)*2;
    std::ifstream ifs;
    std::string file;

    void read_block(){
        ifs.seekg((std::streamoff)(offset+block_bg));
        ifs.read((char*) buffer, (std::streamsize)buffer_size);
        if(ifs.gcount()<(std::streamsize)buffer_size){
            ifs.clear();
        }
        assert(ifs.good());
    }

public:

    explicit bwt_buff_reader(const std::string& input_file, size_t buff_size=1024 * 1024){
        file = input_file;
        ifs.open(file, std::ifstream::binary);
        assert(ifs.good());

        buffer = (char *)malloc(buff_size);
        memset(buffer, 0, buff_size);
        ifs.read((char *)&sb, sizeof(size_t));
        ifs.read((char *)&fb, sizeof(size_t));
        buffer_size = buff_size;
        bpr = sb + fb;
        sec_buff = (char *) malloc(bpr);
        memset(sec_buff, 0, bpr);


        ifs.seekg (0, std::ifstream::end);
        tot_runs = ifs.tellg();
        tot_runs = (tot_runs-offset)/bpr;
        read_block();
    }

    size_t read_sym(size_t i) {
        assert(i<tot_runs);
        size_t sym = 0;
        size_t start = (i*bpr);
        size_t end = start + sb-1;
        size_t buff_start = (start & (buffer_size-1));

        if(start<block_bg || (block_bg+buffer_size)<=end) {

            size_t b_start = start / buffer_size;
            size_t b_end = end / buffer_size;

            if (b_start < b_end) {
                size_t new_block_bg = (start/buffer_size)*buffer_size;
                if(new_block_bg!=block_bg){
                    block_bg = new_block_bg;
                    read_block();
                }

                size_t left = buffer_size-buff_start;
                auto *ptr = (char *)&sym;
                memcpy(ptr, buffer + buff_start, left);
                block_bg = (end / buffer_size) * buffer_size;
                read_block();
                memcpy(ptr+left, buffer, sb-left);
                return sym;
            }else{
                block_bg = (end / buffer_size) * buffer_size;
                read_block();
            }
        }
        memcpy(&sym, buffer+buff_start, sb);
        return sym;
    }

    void read_run(size_t i, size_t& sym, size_t& freq) {

        assert(i<tot_runs);
        sym = 0;
        freq = 0;
        size_t start = (i*bpr);
        size_t end = start + bpr-1;
        size_t buff_start = (start & (buffer_size-1));

        if(start<block_bg || (block_bg+buffer_size)<=end) {

            size_t b_start = start / buffer_size;
            size_t b_end = end / buffer_size;

            if (b_start < b_end) {

                size_t new_block_bg = (start/buffer_size)*buffer_size;
                if(new_block_bg!=block_bg){
                    block_bg = new_block_bg;
                    read_block();
                }

                size_t left = buffer_size-buff_start;
                memcpy(sec_buff, buffer + buff_start, left);

                block_bg = (end/buffer_size)*buffer_size;
                read_block();
                memcpy(sec_buff+left, buffer, bpr-left);

                memcpy(&sym, sec_buff, sb);
                memcpy(&freq, sec_buff+sb, fb);
                return;
            }else{
                block_bg = (end/buffer_size)*buffer_size;
                read_block();
            }
        }
        memcpy(&sym, buffer+buff_start, sb);
        memcpy(&freq, buffer+buff_start+sb, fb);
    }

    size_t bytes_per_rsym() const {
        return sb;
    }

    size_t bytes_per_rlen() const {
        return fb;
    }

    size_t size() const {
       return tot_runs;
    }

    void close(bool del_file=false){

        if(buffer!= nullptr){
            free(buffer);
            buffer = nullptr;
        }
        if(sec_buff!= nullptr){
            free(sec_buff);
            sec_buff = nullptr;
        }
        if(ifs.is_open()) ifs.close();

        if(del_file){
            if(remove(file.c_str())){
                std::cout<<"Error trying to remove the file"<<std::endl;
            }
        }
    }

    ~bwt_buff_reader(){
        close();
    }
};

class bwt_buff_writer {

    char * buffer =nullptr;
    char * sec_buff =nullptr;
    size_t bpr=0; //bytes per run
    size_t sb=0; // bytes for a run symbol
    size_t fb=0; // bytes for a run freq
    size_t block_bg=0; //start of the current block
    size_t tot_runs =0; //total number of runs in the file
    size_t buffer_size=0; //size of the buffer
    size_t offset = sizeof(size_t)*2;//the first offset bytes store the values for fb and sb
    size_t l_sym=0; //symbol of the last run
    size_t l_acc_sym=0; //symbol of the last accessed run
    size_t max_sym=0;
    size_t max_freq=0;
    std::string file;
    std::ofstream ofs;
    std::ifstream ifs;
    bool modified=false;
    std::function<size_t(size_t, size_t)> sub = [](size_t a, size_t b){ return a-b;};
    std::function<size_t(size_t, size_t)> add = [](size_t a, size_t b){ return a+b;};

    inline void mod_freq(size_t idx, size_t new_freq, std::function<size_t(size_t, size_t)>& op) {

        size_t start = (idx*bpr)+sb;
        size_t end = start + fb-1;
        size_t buff_start = (start & (buffer_size - 1));
        size_t freq=0;

        if(start<block_bg || (block_bg + buffer_size) <= end) {

            size_t b_start = start/buffer_size;
            size_t b_end = end/buffer_size;

            if(modified){
                write_block(buffer_size);
            }

            if(b_start<b_end){
                size_t new_block_bg = (start/buffer_size)*buffer_size;
                if(new_block_bg!=block_bg){
                    block_bg = new_block_bg;
                    read_block();
                }

                //retrieve the frequency
                size_t left = buffer_size-buff_start;
                memcpy(&freq, buffer + buff_start, left);

                block_bg = (end/buffer_size)*buffer_size;
                read_block();
                auto *ptr = (char *) &freq;
                memcpy(ptr + left, buffer, (fb-left));
                freq=op(freq, new_freq);
                assert(freq<=max_freq);

                //store the updated frequency
                block_bg = (start/buffer_size)*buffer_size;
                read_block();
                memcpy(buffer + buff_start, &freq, left);
                write_block(buffer_size);
                block_bg =  (end/buffer_size)*buffer_size;
                read_block();

                ptr = (char *) &freq;
                memcpy(buffer, ptr+left ,(fb-left));
                modified = true;
                return;
            }else{
                block_bg = (end/buffer_size)*buffer_size;
                read_block();
            }
        }
        memcpy(&freq, buffer+buff_start, fb);
        freq=op(freq, new_freq);
        assert(freq<=max_freq);
        memcpy(buffer+buff_start, &freq,  fb);
        modified = true;
    }

    inline void write(size_t idx, uint8_t off1, uint8_t off2, size_t new_val) {

        size_t start = (idx*bpr)+off1;
        size_t end = start + off2-1;
        size_t buff_start = (start & (buffer_size - 1));

        if(start<block_bg || (block_bg + buffer_size) <= end) {

            size_t b_start = start/buffer_size;
            size_t b_end = end/buffer_size;

            if(modified){
                write_block(buffer_size);
            }

            if(b_start<b_end){

                size_t new_block_bg = (start/buffer_size)*buffer_size;
                if(new_block_bg!=block_bg){
                    block_bg = new_block_bg;
                    read_block();
                }

                size_t left = buffer_size-buff_start;
                memcpy(buffer + buff_start, &new_val, left);
                write_block(buffer_size);
                block_bg =  (end/buffer_size)*buffer_size;
                read_block();

                auto *ptr = (char *) &new_val;
                memcpy(buffer, ptr+left ,(off2-left));
                modified = true;
                return;
            }else{
                block_bg = (end/buffer_size)*buffer_size;
                read_block();
            }
        }
        memcpy(buffer+buff_start, &new_val,  off2);
        modified = true;
    }


    inline size_t read(size_t idx, uint8_t off1, uint8_t off2) {

        size_t start = (idx*bpr)+off1;
        size_t end = start + off2-1;
        size_t buff_start = (start & (buffer_size - 1));
        size_t val=0;

        if(start<block_bg || (block_bg + buffer_size) <= end) {

            size_t b_start = start/buffer_size;
            size_t b_end = end/buffer_size;

            if(modified){
                write_block(buffer_size);
            }

            if(b_start<b_end){
                size_t new_block_bg = (start/buffer_size)*buffer_size;
                if(new_block_bg!=block_bg){
                    block_bg = new_block_bg;
                    read_block();
                }

                size_t left = buffer_size-buff_start;
                memcpy(&val, buffer + buff_start, left);
                block_bg = (end/buffer_size)*buffer_size;
                read_block();
                auto *ptr = (char *)&val;
                assert(left<off2);
                memcpy(ptr + left, buffer, (off2-left));
                return val;
            }else{
                block_bg = (end/buffer_size)*buffer_size;
                read_block();
            }
        }
        memcpy((char *)&val, buffer+buff_start, off2);
        return val;
    }

    void write_block(size_t n_bytes){
        ofs.seekp((std::streamoff)(offset+block_bg));
        ofs.write((char*) buffer, (std::streamsize) n_bytes);
        //TODO caveat: it might produce desynchronization with the ifs in some platforms (OSx)
        // if I don't flush the ofs buffer. However, flusing the buffer is expensive.
        // Maybe I should define an if macro here.
        ofs.flush();
        assert(ofs.good());
        modified=false;
    }

    void read_block(){
        assert(!modified);
        ifs.seekg((std::streamoff)(offset+block_bg));
        ifs.read((char*) buffer, (std::streamsize)buffer_size);
        if(ifs.gcount()<(std::streamsize)buffer_size){
            ifs.clear();
        }
        assert(ifs.good());
        modified=false;
    }
public:

    explicit bwt_buff_writer(const std::string& input_file, std::ios::openmode mode=std::ios::out,
                             uint8_t _sb=8, uint8_t _fb=8, size_t buff_size=1024 * 1024){

        file = input_file;

        ofs.open(file, mode | std::ios::out | std::ios::binary);
        assert(ofs.good());

        ifs.open(file,  std::ios::binary);
        assert(ifs.good());

        buffer_size = buff_size;
        buffer = (char *)malloc(buffer_size);
        memset(buffer, 0, buffer_size);

        if(mode==std::ios::out){
            sb = _sb;
            fb = _fb;
            ofs.write((char *)&sb, sizeof(size_t));
            ofs.write((char *)&fb, sizeof(size_t));
            bpr = sb+fb;
        }else{
            ifs.read((char *)&sb, sizeof(size_t));
            ifs.read((char *)&fb, sizeof(size_t));
            bpr = sb+fb;

            ifs.seekg (0, std::ifstream::end);
            tot_runs = ifs.tellg();
            tot_runs = (tot_runs-offset)/bpr;
            read_block();
        }

        sec_buff = (char *)malloc(bpr);
        memset(sec_buff, 0, bpr);

        max_freq = (1UL << (fb*8)) - 1;
        max_sym = (1UL << (sb*8)) - 1;
    }

    ~bwt_buff_writer(){
        close();
    };

    void read_run(size_t i, size_t& sym, size_t& freq) {

        assert(i<tot_runs);
        sym = 0;
        freq = 0;
        size_t start = (i*bpr);
        size_t end = start + bpr-1;
        size_t buff_start = (start & (buffer_size-1));

        if(start<block_bg || (block_bg+buffer_size)<=end) {

            if(modified) write_block(buffer_size);

            size_t b_start = start / buffer_size;
            size_t b_end = end / buffer_size;

            if (b_start < b_end) {

                size_t new_block_bg = (start/buffer_size)*buffer_size;
                if(new_block_bg!=block_bg){
                    block_bg = new_block_bg;
                    read_block();
                }

                size_t left = buffer_size-buff_start;
                memcpy(sec_buff, buffer + buff_start, left);

                block_bg = (end/buffer_size)*buffer_size;
                read_block();
                memcpy(sec_buff+left, buffer, bpr-left);

                memcpy(&sym, sec_buff, sb);
                memcpy(&freq, sec_buff+sb, fb);
                return;
            }else{
                block_bg = (end/buffer_size)*buffer_size;
                read_block();
            }
        }
        memcpy(&sym, buffer+buff_start, sb);
        memcpy(&freq, buffer+buff_start+sb, fb);
    }

    void push_back(size_t sym, size_t freq) {

        assert(sym<=max_sym && freq<=max_freq);

        size_t start = (tot_runs*bpr);
        size_t end = start + bpr-1;
        size_t buff_start = (start & (buffer_size - 1));
        tot_runs++;
        l_sym = sym;
        l_acc_sym = sym;

        if(start<block_bg || (block_bg + buffer_size) <= end){
            if(modified) write_block(buffer_size);

            size_t b_start = start/buffer_size;
            size_t b_end = end/buffer_size;

            if(b_start<b_end) {
                size_t new_block_bg = (start/buffer_size)*buffer_size;
                if(new_block_bg!=block_bg){
                    block_bg = new_block_bg;
                    read_block();
                }

                memcpy(sec_buff, &sym, sb);
                memcpy(sec_buff + sb, &freq, fb);
                size_t left = buffer_size - buff_start;
                memcpy(buffer+buff_start, sec_buff, left);
                write_block(buffer_size);

                block_bg = (end/buffer_size)*buffer_size;
                memcpy(buffer, sec_buff+left, bpr-left);
                modified = true;
                return;
            }else{
                block_bg = (end/buffer_size)*buffer_size;
                read_block();
            }
        }
        memcpy(buffer+buff_start, &sym, sb);
        memcpy(buffer+buff_start+sb, &freq, fb);
        modified=true;
    }

    void inc_freq(size_t idx, size_t val){
        mod_freq(idx, val, add);
    }

    void inc_freq_last(size_t val){
        mod_freq(tot_runs-1, val, add);
    }

    void dec_freq(size_t idx, size_t val){
        mod_freq(idx, val, sub);
    }

    size_t read_sym(size_t idx) {
        assert(idx<tot_runs);
        l_acc_sym = read(idx, 0, sb);
        return l_acc_sym;
    }

    size_t read_freq(size_t idx) {
        assert(idx<tot_runs);
        return read(idx, sb, fb);
    }

    void write_sym(size_t idx, size_t new_sym) {
        assert(idx<tot_runs && new_sym<=max_sym);
        write(idx, 0, sb, new_sym);
    }

    void write_freq(size_t idx, size_t new_freq) {
        assert(idx<tot_runs && new_freq<=max_freq);
        write(idx, sb, fb, new_freq);
    }

    size_t size() const {
        return tot_runs;
    }

    size_t last_sym() const {
        return l_sym;
    }

    size_t last_acc_sym() const {
        return l_acc_sym;
    }

    size_t last_freq() {
        return read_freq(tot_runs-1);
    }

    void close(bool del_file=false){
        if(modified) write_block(buffer_size);

        if(ofs.is_open() || ifs.is_open()){
            ofs.close();
            ifs.close();
            auto length = offset + (tot_runs*bpr);
            auto res = truncate(file.c_str(), (off_t)length);
            assert(res==0);
        }

        if(buffer!= nullptr){
            free(buffer);
            buffer = nullptr;
        }

        if(sec_buff!= nullptr){
            free(sec_buff);
            sec_buff = nullptr;
        }

        if(del_file){
            if(remove(file.c_str())){
                std::cout<<"Error trying to remove the file"<<std::endl;
            }
        }
    };
};

//bprs is bytes per run symbol (1, assuming the text has a byte alphabet (0-255))
//bprl is bytes per run len (5, assuming the length of the runs fits 5 bytes)
inline void plain2grlbwt(const std::string& orig_bwt, const std::string& new_bwt,
                         uint8_t bprs=1, uint8_t bprl=5){

    static constexpr size_t buff_size = 8 * 1024 * 1024;
    std::ifstream bwt_in(orig_bwt, std::ios::binary);
    bwt_buff_writer bwt_out(new_bwt, std::ios::out, bprs, bprl);
    const auto buffer = static_cast<uint8_t *>(malloc(buff_size));
    assert(buffer!=nullptr);
    size_t fsz = std::filesystem::file_size(orig_bwt);

    //read the first block separately to set the initial symbol
    bwt_in.read(reinterpret_cast<char *>(buffer), buff_size);
    size_t read_bytes = bwt_in.gcount();
    uint8_t sym = buffer[0];
    size_t len = 1;
    for (size_t i = 1; i < read_bytes; i++) {
        if(sym!=buffer[i]){
            bwt_out.push_back(sym, len);
            sym = buffer[i];
            len = 0;
        }
        len++;
    }
    fsz -= read_bytes;
    //

    while (fsz>0) {
        bwt_in.read(reinterpret_cast<char *>(buffer), buff_size);
        read_bytes = bwt_in.gcount();
        for (size_t i = 0; i < read_bytes; i++) {
            if(sym!=buffer[i]){
                bwt_out.push_back(sym, len);
                sym = buffer[i];
                len = 0;
            }
            len++;
        }
        fsz -= read_bytes;
    }
    assert(len>0);
    bwt_out.push_back(sym, len);
    free(buffer);
    bwt_in.close();
    bwt_out.close();
}

inline void grl2plain(const std::string& rl_file, const std::string& output_plain_file){
    std::ofstream ofs(output_plain_file, std::ios::out | std::ios::binary);
    uint8_t buffer[1024]={0};
    bwt_buff_reader bwt_reader(rl_file);
    size_t sym, freq, k=0;
    size_t sym_freqs[256]={0};
    for(size_t i=0;i<bwt_reader.size();i++){
        bwt_reader.read_run(i, sym, freq);
        for(size_t j=0;j<freq;j++){
            buffer[k++] = sym;
            if(k==1024){
                ofs.write(reinterpret_cast<char *>(buffer), 1024);
                k=0;
            }
        }
        sym_freqs[sym]+=freq;
    }
    if(k!=0){
        ofs.write(reinterpret_cast<char *>(buffer), static_cast<std::streamsize>(k));
    }
    bwt_reader.close();
    ofs.close();
}

#endif //GBWT_BWT_IO_H
