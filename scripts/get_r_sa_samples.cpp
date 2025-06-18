//
// Created by Diaz, Diego on 21.5.2025.
//
#include<iostream>
#include "../include/bwt_io.h"
#include<thread>
#include "../include/utils.h"
#include "custom_wt_rlmn.hpp"

enum run_annot_t{
    NONE=0,
    HEAD=1,
    TAIL=2,
};

struct sa_samp_type{
    uint32_t str=0;
    uint64_t pos=0;
    uint64_t run_id=0;
    run_annot_t annotation=NONE;

    sa_samp_type(uint32_t str_, uint64_t pos_, uint64_t run_id_, run_annot_t annot_): str(str_),
                                                                                        pos(pos_),
                                                                                        run_id(run_id_),
                                                                                        annotation(annot_){}
};

typedef sdsl::custom_wt_rlmn<> index_type;

void compute_samples(size_t n_threads, index_type& bwt, std::string& rsa_file, std::string& str_ranges_file) {

    std::vector<uint64_t> str_lens(bwt.n_strings()+1);
    std::vector<std::vector<sa_samp_type>> thread_sa_samples(n_threads);

    auto lambda_worker = [&](size_t start, size_t end, size_t t) -> void {

        // I am exploiting the fact that the strings in the BCR BWT can be decoded in text order.
        std::vector<sa_samp_type> tmp;
        bool is_tail;
        size_t run;

        for(size_t i=start;i<=end;i++){

            size_t bwt_pos = i;
            uint64_t len = 1;

            while(true){

                run = bwt.pos2run(bwt_pos);

                is_tail = bwt.is_run_tail(bwt_pos);
                if(is_tail){
                    tmp.emplace_back(i, len, run, TAIL);
                }

                bool is_head = bwt.is_run_head(bwt_pos);
                if(is_head){
                    tmp.emplace_back(i, len, run, HEAD);
                }

                auto res = bwt.lf(bwt_pos);
                if(res.first==bwt.sep_sym()){

                    if(!is_tail){
                        tmp.emplace_back(i, len, run, TAIL);
                    }
                    if(!is_head){
                        tmp.emplace_back(i, len, run, HEAD);
                    }

                    assert(res.second<bwt.n_strings());
                    str_lens[i] = len;
                    break;
                }
                len++;
                bwt_pos = res.second;
            }

            std::reverse(tmp.begin(), tmp.end());
            for(auto &samp : tmp){
                samp.pos = len-samp.pos;
                thread_sa_samples[t].push_back(samp);
            }
            tmp.clear();
        }
    };

    size_t n_strings = bwt.n_strings();
    size_t strings_per_thread = INT_CEIL((n_strings), n_threads);
    std::vector<std::thread> working_threads;
    for (size_t i = 0; i < n_threads; i++) {
        size_t start = strings_per_thread * i;
        size_t end = std::min<size_t>(strings_per_thread * (i + 1), n_strings)-1;

        working_threads.emplace_back(lambda_worker, start, end, i);
        if(end == (n_strings-1)) break;
    }

    std::cout<<"Computing the SA samples from the BWT using "<<working_threads.size()<<" working threads"<<std::endl;
    for (auto &working_thread: working_threads) {
        working_thread.join();
    }

    std::vector<sa_samp_type> sa_samples;
    for (auto& vec : thread_sa_samples) {
        sa_samples.insert(sa_samples.end(),
                          std::make_move_iterator(vec.begin()),
                          std::make_move_iterator(vec.end()));
    }

    size_t acc=0, tmp;
    for(auto &len : str_lens){
        tmp = len;
        len = acc;
        acc+=tmp;
    }
    str_lens[n_strings] = acc;

    for(auto &elm : sa_samples){
        elm.pos+=str_lens[elm.str];
        //std::cout<<"pos:"<<elm.pos<<" run:"<<elm.run_id<<" annot:"<<elm.annotation<<std::endl;
    }

    std::cout<<"Sorting the samples"<<std::endl;
    std::sort(sa_samples.begin(), sa_samples.end(), [](const sa_samp_type& a, const sa_samp_type& b) {
        if(a.run_id!=b.run_id){
            return a.run_id<b.run_id;
        }
        if(a.pos!=b.pos){
            return a.pos < b.pos;
        }

        return a.annotation < b.annotation;
    });

    /*for(size_t j=0; j<(sa_samples.size()-1);j++){
        //std::cout<<sa_samples[j].str<<" pos:"<<sa_samples[j].pos<<" run:"<<sa_samples[j].run_id<<" annot:"<<sa_samples[j].annotation<<" next_samp:"<<sa_samples[j].next_samp<<" str_len:"<<str_lens[sa_samples[j].str]<<std::endl;
        if(sa_samples[j].annotation==TAIL){
            assert(sa_samples[j+1].annotation==HEAD);
            sa_samples[j].next_samp = sa_samples[j+1].pos;
        }
    }*/

    //assert((sa_samples.size()/2)==bwt.n_runs());
    //we use 9 bytes per entry : 8 for the position, and 1 for the annotation (HEAD, TAIL, STR_START)

    off_t buff_size = sizeof(uint64_t)*4096;
    auto *buffer = (uint64_t *) malloc(buff_size);
    std::ofstream ofs(rsa_file, std::ios::binary);
    off_t b_pos=0;

    size_t k=0;
    for(auto & sa_sample : sa_samples){
        /*if(k>=117142022 && k<=117142024){
            std::cout<<k<<" "<<sa_sample.pos<<" "<<sa_sample.run_id<<" "<<sa_sample.annotation<<std::endl;
        }*/
        buffer[b_pos] = sa_sample.pos;
        b_pos++;
        if(b_pos==4096){
            ofs.write((char *)buffer, buff_size);
            b_pos=0;
        }
        k++;
    }
    if(b_pos!=0){
        buff_size = off_t(b_pos*sizeof(uint64_t));
        ofs.write((char *)buffer, buff_size);
    }
    ofs.close();
    free(buffer);

    std::ofstream ofs2(str_ranges_file, std::ios::binary);
    ofs2.write((char *)str_lens.data(), off_t(sizeof(uint64_t)*str_lens.size()));

    std::cout<<"Total number of SA samples collected "<<sa_samples.size()<<std::endl;
    std::cout<<"Extra samples collected due to BCR  "<<sa_samples.size()-(bwt.n_runs()*2)<<std::endl;
}

int main(int argc, char** argv) {
    if(argc!=4){
        std::cout<<"usage: ./get_r_sa_samples file.rlbwt n_threads output_file\n\n"
                    "file.rlbwt: the BCR BWT file generated by grlBWT\n"
                    "n_threads: number of threads to compute the samples\n"
                    "output_file: file where the samples will be stored"<< std::endl;
        exit(0);
    }
    errno = 0;
    std::string input_file = std::string(argv[1]);
    char *endptr;
    uint64_t n_threads = strtol(argv[2], &endptr, 10);
    if (*endptr != '\0') {
        fprintf(stderr, "Invalid number of threads: %s\n", argv[1]);
        return 1;
    }
    std::string rsa_file = std::string(argv[3])+".sa_samples";
    std::string str_ranges_file = std::string(argv[3])+".str_ranges";
    std::cout<<"Building the RLBWT ..."<<std::endl;
    index_type bwt(input_file);
    std::cout<<"Input file:"<<input_file<<" has "<<bwt.n_strings()<<" strings and "<<bwt.n_runs()<<" BWT runs"<<std::endl;
    compute_samples(n_threads, bwt, rsa_file, str_ranges_file);
    std::cout<<"SA samples were stored in "<<rsa_file<<std::endl;
    std::cout<<"The ranges of the strings in the text were stored in "<<str_ranges_file<<std::endl;
}
