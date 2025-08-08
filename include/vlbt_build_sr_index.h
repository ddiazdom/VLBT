//
// Created by Diaz, Diego on 29.5.2025.
//

#ifndef VLBT_BUILD_SR_INDEX_H
#define VLBT_BUILD_SR_INDEX_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <filesystem>

#include "vlbt_build_phi.h"
#include "vlbt_build_bwt.h"

struct sample_type {
    uint64_t head_val;
    uint64_t prev_tail_val;
    uint64_t next_head_val;
    uint64_t run_id:56;
    bool is_head_sampled:8;
};

template<class size_type>
struct buff_writer {

    std::ofstream ofs;
    size_t buff_pos=0;
    std::vector<size_type> buffer;
    size_t buffer_size = 1024*1024;
    size_t sz=0;
    bool is_open = true;

    explicit buff_writer(const std::string& output_file) {
        ofs = std::ofstream(output_file, std::ios::binary);
        buffer.resize(buffer_size);
    }

    void push_back(size_t diff, bool is_diff_neg, size_t valid_area, size_t len) {
        assert(is_open);
        buffer[buff_pos] = diff;
        buffer[buff_pos] = (buffer[buff_pos]<<1) | is_diff_neg;//we use 1 bit to encode if the difference is negative
        buffer[buff_pos] = (buffer[buff_pos]<<15) | valid_area;//we use 15 bits to encode the valid area
        buff_pos++;
        buffer[buff_pos++] = len;
        if(buff_pos==buffer_size){
            ofs.write((char *)buffer.data(), sizeof(size_type)*buffer_size);
            buff_pos=0;
        }
        sz+=len;
    }

    inline size_t size() const {
        return sz;
    }

    void close() {
        if (buff_pos!=0) {
            ofs.write((char *)buffer.data(), sizeof(size_type)*buff_pos);
        }
        ofs.close();
        is_open = false;
    }

    ~buff_writer() {
        close();
    }
};

static inline uint64_t get_diff(const uint64_t first, const uint64_t second) {
    const uint64_t abs_diff = (first > second) ? (first - second) : (second - first);
    assert(abs_diff<=INT64_MAX);
    return abs_diff;
}

template<class size_type>
void get_head_samples(std::vector<sample_type>& samples, std::string& str_ranges_file,
                      size_t ssamp_val, std::string& ssamp_heads_file){

    std::cout<<"Sampling the heads"<<std::endl;

    //sort the samples by the text position of the tails
    std::sort(samples.begin(), samples.end(), [](auto const& a, auto const&b){
        return a.prev_tail_val<b.prev_tail_val;
    });

    //sort the samples by text position
    std::sort(samples.begin(), samples.end(), [](auto const& a, auto const&b){
        return a.head_val<b.head_val;
    });
    //

    //compute and store the run head subsamples for phi
    size_t f_size = std::filesystem::file_size(str_ranges_file);
    size_t n_elements = f_size/sizeof(size_type);
    std::vector<size_type> str_ranges(n_elements, 0);
    std::ifstream ifs_str_ranges(str_ranges_file, std::ios::binary);
    ifs_str_ranges.read((char *)str_ranges.data(), off_t(f_size));

    size_t n_strings = n_elements-1, last_sampled, n_samp=0, len, acc_len=0;
    size_type str_boundary;
    size_t s_pos=0;
    for(size_t str=0;str<n_strings;str++){

        assert(samples[s_pos].head_val==str_ranges[str]);
        str_boundary = str_ranges[str+1]-1;
        last_sampled = s_pos;
        s_pos++;
        assert(s_pos<samples.size());

        while((s_pos+1)<samples.size() && samples[s_pos+1].head_val<=str_boundary){


            if((samples[s_pos+1].head_val-samples[last_sampled].head_val>ssamp_val)){
                len = samples[s_pos].head_val-samples[last_sampled].head_val;

                //std::cout<<"sampled "<<samples[last_sampled].head_val<<std::endl;
                //if (samples[last_sampled].head_val==430241310) {
                //    std::cout<<"holaa"<<std::endl;
                //}

                //samples[last_sampled].valid_area = 0;
                //if(last_sampled<(s_pos-1)) {
                //  samples[last_sampled].valid_area = samples[last_sampled+1].head_val-samples[last_sampled].head_val;
                    //TODO the value below is the real "invalid suffix" anything after this suffix is valid
                    //std::cout<<samples[s_pos].head_val-samples[last_sampled+1].head_val<<std::endl;
                //}

                //for(size_t k=last_sampled+1;k<s_pos;k++){
                //    std::cout<<"Not sampled "<<samples[k].head_val<<" "<<samples[k].head_val-samples[last_sampled].head_val<<" "<<std::endl;;
                //}

                samples[last_sampled].is_head_sampled = true;
                acc_len+=len;
                last_sampled = s_pos;
                n_samp++;
            }
            s_pos++;
        }

        if(samples[s_pos].head_val<=str_boundary){
            if((str_boundary+1)-samples[last_sampled].head_val){
                len = samples[s_pos].head_val-samples[last_sampled].head_val;
                samples[last_sampled].is_head_sampled = true;
                acc_len+=len;
                last_sampled = s_pos;
                n_samp++;
            }
            s_pos++;
        }

        samples[last_sampled].is_head_sampled = true;
        len = (str_boundary+1)-samples[last_sampled].head_val;
        //samples[last_sampled].valid_area = 0;
        acc_len += len;
        n_samp++;

        assert(acc_len==str_ranges[str+1]);
    }

    //sort by run
    std::sort(samples.begin(), samples.end(), [](sample_type const& a, sample_type const& b){
        return a.run_id < b.run_id;
    });
    std::ofstream ofs_ssamp_heads(ssamp_heads_file, std::ios::binary);

    size_t buffer_size = 1024*1024, rem = samples.size(), n_blocks = samples.size()/buffer_size;
    std::vector<size_type> buffer(buffer_size, 0);
    s_pos = 0;
    size_type discard_mark = std::numeric_limits<uint64_t>::max();

    for(size_t i=0;i<n_blocks;i++){
        for(size_t j=0;j<buffer_size;j++){
            buffer[j] = samples[s_pos].is_head_sampled ? samples[s_pos].head_val : discard_mark;
            s_pos++;
        }
        rem-=buffer_size;
        ofs_ssamp_heads.write((char *)buffer.data(), buffer_size*sizeof(size_type));
    }
    assert(rem<buffer_size);
    for(size_t j=0;j<rem;j++){
        buffer[j] = samples[s_pos].is_head_sampled ? samples[s_pos].head_val : discard_mark;
        s_pos++;
    }
    assert(s_pos==samples.size());
    ofs_ssamp_heads.write((char *)buffer.data(), rem*sizeof(size_type));
    assert(ofs_ssamp_heads.tellp()/sizeof(size_type)==samples.size());
    ofs_ssamp_heads.close();

    //a technical hack: samples[0].prev_tail_val contains the SA value of the last run in the BWT.
    //the difference of this value with the value of the following head should be 0 because there is no head.
    samples[0].head_val = samples[0].prev_tail_val;
    std::cout<<"We subsampled "<<n_samp<<" elements out of "<<samples.size()<<" ("<<double(n_samp)/double(samples.size())*100<<"%)"<<std::endl;
    //
}

template<class size_type>
void get_tail_samples(std::vector<sample_type>& samples,
                      std::string& str_ranges_file,
                      size_t subsamp_step,
                      std::string& ssamp_tails_file){

    std::cout<<"Sampling the tails"<<std::endl;

    //sort by the run tails
    std::sort(samples.begin(), samples.end(), [](sample_type const& a, sample_type const& b){
        return a.prev_tail_val < b.prev_tail_val;
    });

    //compute and store the run head subsamples for phi
    size_t f_size = std::filesystem::file_size(str_ranges_file);
    size_t n_elements = f_size/sizeof(size_type);
    std::vector<size_type> str_ranges(n_elements, 0);
    std::ifstream ifs_str_ranges(str_ranges_file, std::ios::binary);
    ifs_str_ranges.read((char *)str_ranges.data(), f_size);

    buff_writer<size_type> output_buffer(ssamp_tails_file);
    const size_t n_strings = n_elements-1;
    size_t len, s_pos=0, diff, valid_area;
    size_type str_boundary;
    bool is_diff_neg;

    for(size_t str=0;str<n_strings;str++){

        assert(samples[s_pos].prev_tail_val==str_ranges[str]);
        str_boundary = str_ranges[str+1]-1;
        size_t last_sampled = s_pos;
        s_pos++;
        assert(s_pos<samples.size());

        while(s_pos<samples.size() && samples[s_pos].prev_tail_val<=str_boundary){

            if(samples[s_pos].is_head_sampled){

                assert(s_pos>last_sampled);
                const size_t n_blocks = s_pos-last_sampled;

                valid_area = samples[last_sampled+1].prev_tail_val-samples[last_sampled].prev_tail_val;
                len = samples[s_pos].prev_tail_val-samples[last_sampled].prev_tail_val;
                diff = get_diff(samples[last_sampled].prev_tail_val, samples[last_sampled].head_val);
                is_diff_neg = samples[last_sampled].prev_tail_val>samples[last_sampled].head_val;

                //if (samples[last_sampled].prev_tail_val<=72440562 && 72440561<samples[last_sampled].prev_tail_val+len) {
                //    std::cout<<"holaa"<<std::endl;
                //}

                if (n_blocks==1) {
                    output_buffer.push_back(diff, is_diff_neg, 0, len);
                } else if (valid_area<subsamp_step) {
                    output_buffer.push_back(diff, is_diff_neg, valid_area, len);
                } else {
                    //we break the block into two pieces
                    len = samples[last_sampled+1].prev_tail_val-samples[last_sampled].prev_tail_val;
                    assert(len>=subsamp_step && valid_area==len);
                    diff = get_diff(samples[last_sampled].prev_tail_val, samples[last_sampled].head_val);
                    is_diff_neg = samples[last_sampled].prev_tail_val>samples[last_sampled].head_val;
                    output_buffer.push_back(diff, is_diff_neg, 0, len);

                    last_sampled++;
                    valid_area = samples[last_sampled+1].prev_tail_val-samples[last_sampled].prev_tail_val;
                    assert(valid_area<subsamp_step);
                    len = samples[s_pos].prev_tail_val-samples[last_sampled].prev_tail_val;
                    diff = get_diff(samples[last_sampled].prev_tail_val, samples[last_sampled].head_val);
                    is_diff_neg = samples[last_sampled].prev_tail_val>samples[last_sampled].head_val;
                    output_buffer.push_back(diff, is_diff_neg, valid_area, len);
                }
                last_sampled = s_pos;
            }
            s_pos++;
        }

        const size_t n_blocks = s_pos-last_sampled;
        valid_area = samples[last_sampled+1].prev_tail_val-samples[last_sampled].prev_tail_val;
        len = (str_boundary+1)-samples[last_sampled].prev_tail_val;
        diff = get_diff(samples[last_sampled].prev_tail_val, samples[last_sampled].head_val);
        is_diff_neg = samples[last_sampled].prev_tail_val>samples[last_sampled].head_val;

        if (n_blocks==1) {
            output_buffer.push_back(diff, is_diff_neg, 0, len);
        } else if(valid_area<subsamp_step) {
            output_buffer.push_back(diff, is_diff_neg, valid_area, len);
        } else {
            len = samples[last_sampled+1].prev_tail_val-samples[last_sampled].prev_tail_val;
            assert(len>=subsamp_step && valid_area==len);
            diff = get_diff(samples[last_sampled].prev_tail_val, samples[last_sampled].head_val);
            is_diff_neg = samples[last_sampled].prev_tail_val>samples[last_sampled].head_val;
            output_buffer.push_back(diff, is_diff_neg, 0, len);

            last_sampled++;
            valid_area = samples[last_sampled+1].prev_tail_val-samples[last_sampled].prev_tail_val;
            assert(valid_area<subsamp_step);
            len = (str_boundary+1)-samples[last_sampled].prev_tail_val;
            diff = get_diff(samples[last_sampled].prev_tail_val, samples[last_sampled].head_val);
            is_diff_neg = samples[last_sampled].prev_tail_val>samples[last_sampled].head_val;
            output_buffer.push_back(diff, is_diff_neg, valid_area, len);
        }
        assert(output_buffer.size()==str_ranges[str+1]);
    }
}

template<class size_type>
void subsample_sa_samples(std::string& sa_samples_file, std::string& str_ranges_file,
                          size_t ssamp_step, std::string& ssamp_heads_file,
                          std::string& ssamp_tail_file){

    size_t n_elements = std::filesystem::file_size(sa_samples_file)/sizeof(size_type);
    std::cout<<"Subsampling SA samples"<<std::endl;
    std::vector<sample_type> samples(n_elements/2);
    size_t s_pos=0;

    std::ifstream ifs_orig_samples(sa_samples_file, std::ios::binary);
    static constexpr size_t buffer_size = 1024*1024;
    std::vector<size_type> buffer(buffer_size+1, 0);
    size_t n_blocks = n_elements/buffer_size;
    size_t rem = n_elements;

    //read samples from the disk and reorganize them
    size_type discard_mark = std::numeric_limits<size_type>::max();
    size_type prev_tail_val = discard_mark;
    for(size_t i=0;i<n_blocks;i++){
        ifs_orig_samples.read((char *)buffer.data(), sizeof(size_type) * buffer_size);
        for(size_t j=0;j<buffer_size;j+=2){
            samples[s_pos].prev_tail_val = prev_tail_val;
            samples[s_pos].head_val = buffer[j];
            samples[s_pos].run_id = s_pos;
            ++s_pos;
            prev_tail_val = buffer[j+1];
        }
        rem -=buffer_size;
    }

    // We artificially associate the tail of the last run in the BWT with the head of the first run.
    // It is just a technicality to have the value of that tail in the vector samples
    if(rem>0){
        ifs_orig_samples.read((char *)buffer.data(), off_t(sizeof(size_type)*rem));
        for(size_t j=0;j<rem;j+=2){
            samples[s_pos].prev_tail_val = prev_tail_val;
            samples[s_pos].head_val = buffer[j];
            samples[s_pos].run_id = s_pos;
            ++s_pos;
            prev_tail_val = buffer[j+1];
        }
        samples[0].prev_tail_val = buffer[rem-1];
    }else{
        samples[0].prev_tail_val = buffer[buffer_size-1];
    }

    assert(s_pos==samples.size());
    assert(s_pos==(n_elements/2));
    ifs_orig_samples.close();
    //

    get_head_samples<size_type>(samples, str_ranges_file, ssamp_step,  ssamp_heads_file);
    get_tail_samples<size_type>(samples, str_ranges_file, ssamp_step, ssamp_tail_file);
}

template<class sr_index_type, class sa_samp_type>
void build_sr_index(sr_index_type& index, std::string& input_prefix, size_t subsamp_step, std::string& output_prefix){

    using bwt_type = typename sr_index_type::bwt_t;
    using phi_type = typename sr_index_type::phi_t;

    std::string samp_sa_file = input_prefix+".sa_samples";
    std::string str_ranges_file = input_prefix+".str_ranges";
    std::string bwt_file = input_prefix+".ebwt";

    std::string ssamp_heads_file = output_prefix+".ssamp_heads";
    std::string ssamp_tails_file = output_prefix+".ssamp_tails";

    index.phi.subsamp_step = subsamp_step;

    subsample_sa_samples<sa_samp_type>(samp_sa_file, str_ranges_file, subsamp_step, ssamp_heads_file, ssamp_tails_file);

    build_bwt_th<bwt_type, sa_samp_type>(index.bwt, bwt_file, subsamp_step, GRL_BWT, ssamp_heads_file);
    build_phi<phi_type, sa_samp_type>(index.phi, ssamp_tails_file);
}
#endif //VLBT_BUILD_SR_INDEX_H