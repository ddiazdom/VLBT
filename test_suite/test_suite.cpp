//
// Created by Diaz, Diego on 17.10.2022.
//

#include <iostream>
#include <ostream>

//==== VLBT framework
#include "../include/vlbt_build_bwt.h"
#include "../include/vlbt_build_phi.h"
#include "../include/vlbt_build_sr_index.h"
#include "../include/vlbt_bwt.h"
#include "../include/vlbt_phi.h"
#include "../include/vlbt_sr_index.h"
//=====

#include "simple_rlbwt.h"
#include "r-index/internal/r_index.hpp"
#include "../scripts/utils.h"
#include <unordered_set>
#include <random>

void print_histogram(std::vector<double>& times)  {

    std::sort(times.begin(), times.end());

    const auto q1 = times[times.size() / 4];
    const auto q3 = times[(3 * times.size()) / 4];
    const double iqr = q3 - q1;
    const double upper = q3 + 1.5 * iqr;

    std::vector<double> filtered;
    for (double t : times) {
        if (t <= upper) {
            filtered.push_back(t);
        }
    }

    // Find range
    const double min_val = filtered[0];
    const double max_val = filtered.back();

    // Create histogram bins
    constexpr int num_bins = 20;
    std::vector bins(num_bins, 0);
    const double bin_width = (max_val - min_val) / num_bins;

    for (double t : filtered) {
        int bin = std::min(int((t - min_val) / bin_width), num_bins - 1);
        bins[bin]++;
    }

    // Print ASCII histogram
    std::cout << "Histogram of runtimes (microseconds):\n";
    for (int i = 0; i < num_bins; ++i) {
        double bin_start = min_val + i * bin_width;
        double bin_end = bin_start + bin_width;
        std::cout << "[" << bin_start << ", " << bin_end << "): ";

        int count = bins[i];
        for (int j = 0; j < static_cast<int>(count * 50 / filtered.size()); ++j) { // scale to max width 50
            std::cout << '#';
        }
        std::cout << " (" << count << ")\n";
    }
}

std::vector<uint64_t> sample_unique(uint64_t n, uint64_t x) {
    if (x > n) throw std::invalid_argument("x cannot be larger than n");

    std::unordered_set<uint64_t> seen;
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist(0, n - 1);

    while (seen.size() < x) {
        uint64_t val = dist(rng);
        seen.insert(val);  // insert is a no-op if val already exists
    }

    return {seen.begin(), seen.end()};
}

namespace std {
    template<typename X, typename Y>
    struct hash<std::pair<X, Y>> {
        std::size_t operator()(const std::pair<X, Y> &pair) const {
            return std::hash<X>()(pair.first) ^ std::hash<Y>()(pair.second);
        }
    };
}

std::vector<std::pair<uint64_t, uint8_t>> compute_random_rank_queries(uint64_t text_size, uint8_t alphabet_size, uint64_t n_samps) {

    if(n_samps > text_size*alphabet_size) throw std::invalid_argument("invalid sample size");

    std::unordered_set<std::pair<uint64_t, uint8_t>> seen;
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist_idx(0, text_size - 1);
    std::uniform_int_distribution<uint64_t> dist_sym(0, alphabet_size- 1);

    while (seen.size() < n_samps) {
        uint64_t idx = dist_idx(rng);
        uint64_t symbol = dist_sym(rng);
        seen.insert({idx, symbol});
    }
    return {seen.begin(), seen.end()};
}

#define MEASURE(query, time_answer, query_answer, time_unit) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query_answer = query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<time_unit>( t2 - t1 ).count();\
}

template<class my_bwt_type, class other_bwt_type>
void test_count_with_sa_head(my_bwt_type& my_bwt, const std::string& my_dt_name,
                             other_bwt_type& other_index, const std::string& other_index_name,
                             const std::string& pat_file) {

    std::cout<<"Testing count with head (microsecs/pat and microsecs/occ)"<<std::endl;

    ulint n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);
    std::cout<<"\tSearching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    double other_acc_time=0;
    size_t other_acc_count=0;
    size_t j=0;
    std::vector<ulint> head(n_pats);
    for(auto &p : pat_list) {
        std::vector<ulint> tmp_ans;
        MEASURE(other_index.locate_all(p), other_acc_time, tmp_ans, std::chrono::microseconds)
        other_acc_count+=tmp_ans.size();
        head[j] = tmp_ans.back();
        j++;
    }

    double my_acc_time=0;
    size_t my_acc_count=0;
    j=0;
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> my_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(my_bwt.count_with_head(p), my_acc_time, my_ans[j], std::chrono::microseconds)
        my_acc_count+=std::get<1>(my_ans[j])-std::get<0>(my_ans[j])+1;
        j++;
    }
    assert(my_acc_count==other_acc_count);

    std::cout<<"\t"<<my_dt_name<<": ("<<my_acc_time/double(n_pats)<<", "<<my_acc_time/double(my_acc_count)<<"), ";
    //std::cout<<"\t"<<other_index_name<<": ("<<other_acc_time/double(n_pats)<<", "<<other_acc_time/double(other_acc_count)<<")"<<std::endl;;
    std::cout<<"\tTotal occurrences: "<<my_acc_count<<"="<<other_acc_count<<std::endl;

    for(size_t i=0;i<pat_list.size();i++){
        if(std::get<2>(my_ans[i])!=head[i]){
            std::cout<<"Pattern["<<i<<"]: \""<<pat_list[i]<<"\" SA val:"<<std::get<2>(my_ans[i])<<"!="<<head[i]<<std::endl;
            exit(1);
        }
    }
}

template<class my_dt_type, class other_dt_type>
void test_count(my_dt_type& my_dt, const std::string& my_dt_name,
                other_dt_type& other_dt, const std::string& other_dt_name,
                const std::string& pat_file){

    std::cout<<"Testing count (microsecs/pat and microsecs/occ)"<<std::endl;

    ulint n_pats, pat_len;
    const std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);
    std::cout<<"\tSearching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;

    double other_acc_time=0;
    size_t other_acc_count=0;
    size_t j=0;
    std::vector<std::pair<uint64_t, uint64_t>> other_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(other_dt.count(p), other_acc_time, other_ans[j], std::chrono::microseconds)
        other_acc_count+=other_ans[j].second-other_ans[j].first+1;
        j++;
    }

    double my_acc_time=0;
    size_t my_acc_count=0;
    j=0;
    std::vector<std::pair<uint64_t, uint64_t>> my_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(my_dt.count(p), my_acc_time, my_ans[j], std::chrono::microseconds)
        my_acc_count+=my_ans[j].second-my_ans[j].first+1;
        j++;
    }
    std::cout<<"\t"<<my_dt_name<<": ("<<my_acc_time/double(n_pats)<<", "<<my_acc_time/double(my_acc_count)<<"), ";
    std::cout<<"\t"<<other_dt_name<<": ("<<other_acc_time/double(n_pats)<<", "<<other_acc_time/double(other_acc_count)<<")"<<std::endl;;
    std::cout<<"\tTotal occurrences: "<<my_acc_count<<"="<<other_acc_count<<std::endl;

    for(size_t i=0;i<pat_list.size();i++){
        if(my_ans[i].first!=other_ans[i].first || my_ans[i].second!=other_ans[i].second){
            std::cout<<"Pattern["<<i<<"]: \""<<pat_list[i]<<"\" coords:"<<my_ans[i].first<<"!="<<other_ans[i].first <<" or "<<my_ans[i].second<<"!="<<other_ans[i].second<<std::endl;
            exit(1);
        }
    }
}

template<class my_index_type, class other_index_type>
void test_locate(my_index_type& my_dt, const std::string& my_dt_name,
                 other_index_type& other_index, const std::string& other_index_name,
                 const std::string& pat_file){

    std::cout<<"Testing locate (microsecs/pat and microsecs/occ)"<<std::endl;
    ulint n_pats, pat_len;
    const std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);
    std::cout<<"\tSearching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;

    double other_acc_time=0;
    size_t other_acc_count=0;
    size_t j=0;
    std::vector<std::vector<ulint>> other_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(other_index.locate_all(p), other_acc_time, other_ans[j], std::chrono::microseconds)
        other_acc_count+=other_ans[j].size();
        std::reverse(other_ans[j].begin(), other_ans[j].end());
        j++;
    }

    double my_acc_time=0;
    size_t my_acc_count=0;
    j=0;
    std::vector<std::vector<uint64_t>> my_ans(n_pats);
    for(const std::string& p : pat_list) {
        MEASURE(my_dt.locate(p), my_acc_time, my_ans[j], std::chrono::microseconds)
        my_acc_count+=my_ans[j].size();
        j++;
    }
    std::cout<<"\t"<<my_dt_name<<": ("<<my_acc_time/double(n_pats)<<", "<<my_acc_time/double(my_acc_count)<<"), ";
    std::cout<<"\t"<<other_index_name<<": ("<<other_acc_time/double(n_pats)<<", "<<other_acc_time/double(other_acc_count)<<")"<<std::endl;;
    std::cout<<"\tTotal occurrences: "<<my_acc_count<<"="<<other_acc_count<<std::endl;

    for(size_t i=0;i<pat_list.size();i++){
        assert(my_ans[i].size()==other_ans[i].size());
        size_t n_errors=0;
        for(size_t k=0;k<my_ans[i].size();k++){
            if (my_ans[i][k]!=other_ans[i][k]) {
                std::cout<<"Pattern["<<i<<"]: \""<<pat_list[i]<<"\" occ_idx="<<k<<", my_ans="<<my_ans[i][k]<<", other_ans="<<other_ans[i][k]<<", n_occ="<<my_ans[i].size()<<std::endl;
                n_errors++;
            }
        }
        if (n_errors>0) {
            std::cout<<""<<std::endl;
        }
    }
}

template<class my_bwt_type, class other_dt_type>
void test_access(my_bwt_type& my_dt, std::string my_dt_name, other_dt_type& other_dt, std::string other_dt_name){

    std::cout<<"Testing access (avg_time in nanoseconds)"<<std::endl;
    constexpr size_t samp_size = 1000000;
    std::vector<uint64_t> samples = sample_unique(other_dt.size(), samp_size);

    double acc_time=0;
    std::vector<uint8_t> my_dt_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(my_dt[samples[j]], acc_time, my_dt_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"\t"<<my_dt_name<<": "<<acc_time/static_cast<double>(samp_size)<<", ";

    acc_time=0;
    std::vector<uint8_t> other_dt_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(other_dt[samples[j]], acc_time, other_dt_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"\t"<<other_dt_name<<": "<<acc_time/static_cast<double>(samp_size)<<std::endl;

    for(size_t j=0;j<samples.size();j++){
        if(other_dt_ans[j]!=my_dt_ans[j]){
            std::cout<<my_dt_name<<" idx: "<<samples[j]<<" |\t sym: "<<int(other_dt_ans[j])<<std::endl;
            std::cout<<other_dt_name<<" idx: "<<samples[j]<<" |\t sym: "<<int(my_dt_ans[j])<<std::endl;
        }
        assert(other_dt_ans[j]==my_dt_ans[j]);
    }
}

template<class my_bwt_type, class other_dt_type>
void test_rank(my_bwt_type& my_dt, const std::string& my_dt_name, other_dt_type other_dt, const std::string other_dt_name){

    std::cout<<"Testing rank (avg_time in nanoseconds)"<<std::endl;

    constexpr size_t samp_size = 1000000;
    std::vector<std::pair<uint64_t, uint8_t>> tests = compute_random_rank_queries(other_dt.size(), other_dt.alphabet(), samp_size);

    double acc_time=0;
    std::vector<int64_t> my_dt_ans(samp_size);
    for(size_t j=0;j<tests.size();j++){
        MEASURE(my_dt.rank(tests[j].first, my_dt.eff2byte(tests[j].second)), acc_time, my_dt_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"\t"<<my_dt_name<<": "<<acc_time/static_cast<double>(samp_size)<<", ";

    acc_time=0;
    std::vector<int64_t> other_dt_ans(samp_size);
    for(size_t j=0;j<tests.size();j++){
        MEASURE(other_dt.rank(tests[j].first, my_dt.eff2byte(tests[j].second)), acc_time, other_dt_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"\t"<<other_dt_name<<": "<<acc_time/static_cast<double>(samp_size)<<std::endl;

    for(size_t j=0;j<tests.size();j++){
        if(my_dt_ans[j]<0) continue;
        if(other_dt_ans[j]!=my_dt_ans[j]){
            std::cout<<"query:  idx:"<<tests[j].first<<", sym:"<<int(my_dt.eff2byte(tests[j].second))<<", test_id:"<<j<<std::endl;
            std::cout<<my_dt_name+" rank answer: "<<other_dt_ans[j]<<std::endl;
            std::cout<<other_dt_name+" rank answer: "<<my_dt_ans[j]<<"\n"<<std::endl;
        }
        assert(other_dt_ans[j]==my_dt_ans[j]);
    }
}

template<class my_bwt_type, class other_bwt_type>
void test_inverse_select(my_bwt_type& my_dt, std::string my_dt_name, other_bwt_type& other_dt, std::string other_dt_name){

    std::cout<<"Testing inverse select (avg_time in nanoseconds)"<<std::endl;

    constexpr size_t samp_size = 1000000;
    std::vector<uint64_t> samples = sample_unique(other_dt.size(), samp_size);

    double acc_time=0;
    std::vector<std::pair<uint8_t, uint64_t>> my_dt_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(my_dt.inverse_select(samples[j]), acc_time, my_dt_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"\t"<<my_dt_name<<": "<<acc_time/static_cast<double>(samp_size)<<", ";

    acc_time=0;
    std::vector<std::pair<uint8_t, uint64_t>> other_dt_ans(samp_size);
    for(size_t j=0;j<samples.size();j++){
        MEASURE(other_dt.inverse_select(samples[j]), acc_time, other_dt_ans[j], std::chrono::nanoseconds);
    }
    std::cout<<"\t"<<other_dt_name<<": "<<acc_time/static_cast<double>(samp_size)<<std::endl;

    for(size_t j=0;j<samples.size();j++){
        if(other_dt_ans[j].first!=my_dt_ans[j].first ||
           other_dt_ans[j].second!=my_dt.eff2byte(my_dt_ans[j].second)){
            std::cout<<other_dt_name+" idx: "<<samples[j]<<" |\t sym: "<<int(other_dt_ans[j].second)<<" rank: "<<int(other_dt_ans[j].first)<<std::endl;
            std::cout<<my_dt_name+" idx: "<<samples[j]<<" |\t sym: "<<int(my_dt.eff2byte(my_dt_ans[j].second))<<" rank: "<<int(my_dt_ans[j].first)<<"\n"<<std::endl;
        }
        assert(other_dt_ans[j].first==my_dt_ans[j].first &&
               other_dt_ans[j].second==my_dt.eff2byte(my_dt_ans[j].second));
    }
}

void test_bwt(const std::string& input_prefix, const BWT_FORMAT bwt_file_fmt, const std::string& output_prefix){

    vlbt_rlbwt<4096> bwt_dt;
    const std::string bwt_file = input_prefix+".bwt";

    std::cout<<"Building VLBT from input BWT "<<bwt_file<<std::endl;
    build_bwt(bwt_dt, bwt_file, bwt_file_fmt);
    const std::string output_file = output_prefix+".vlbt_bwt";
    const size_t written_bytes = store_to_file(output_file, bwt_dt);
    std::cout<<"We store "<<written_bytes<<" bytes in "<<output_file<<"\n"<<std::endl;

    //=====
    //Create the RLBWT in case it does not exist
    std::string rlwt_file = input_prefix+".wt_rlmn";
    if(!std::filesystem::exists(rlwt_file)){
        std::cout<<"Creating the RLBWT for "<<bwt_file<<" to compare against VLBT "<<std::endl;
        sdsl::wt_rlmn<> wt;
        sdsl::construct(wt, bwt_file, 1);
        simple_rlbwt simple_rlbwt(wt);
        sdsl::store_to_file(simple_rlbwt, rlwt_file);
        std::cout<<"Wavelet trees created"<<std::endl;
    }
    std::cout<<"Loading RLBWT from "<<rlwt_file<<std::endl;
    simple_rlbwt<sdsl::wt_rlmn<>> other_bwt;
    sdsl::load_from_file(other_bwt, rlwt_file);

    std::string query_pat_file = input_prefix+".pats";
    assert(std::filesystem::exists(query_pat_file));
    //=====

    test_rank(bwt_dt,"rlbwt_vlbt", other_bwt, "wt_rlmn");
    test_count(bwt_dt, "rlbwt_vlbt", other_bwt, "wt_rlmn", query_pat_file);
    test_access(bwt_dt, "rlbwt_vlbt", other_bwt, "wt_rlmn");
    test_inverse_select(bwt_dt, "rlbwt_vlbt", other_bwt, "wt_rlmn");
}

template<class sa_samp_type>
void test_bwt_th(std::string& input_prefix, const BWT_FORMAT bwt_file_fmt, size_t sri_samp_val, const std::string& output_prefix){

    std::cout<<"Testing VLBT BWT with toeholds"<<std::endl;

    vlbt_rlbwt_th<65536> bwt_th_dt;
    build_bwt_th<sa_samp_type>(bwt_th_dt, input_prefix, bwt_file_fmt, sri_samp_val);

    std::string output_file = output_prefix+".rlbwt_th_vlbt";
    size_t written_bytes = store_to_file(output_file, bwt_th_dt);
    std::cout<<"We store "<<written_bytes<<" bytes in "<<output_file<<std::endl;

    //=====
    //Create the RLBWT in case it does not exist
    std::string rlwt_file = input_prefix+".wt_rlmn";
    if(!std::filesystem::exists(rlwt_file)){
        std::string bwt_file = input_prefix+".bwt";
        std::cout<<"Creating the RLBWT for "<<bwt_file<<" to compare against VLBT "<<std::endl;
        sdsl::wt_rlmn<> wt;
        sdsl::construct(wt, bwt_file, 1);
        simple_rlbwt simple_rlbwt(wt);
        sdsl::store_to_file(simple_rlbwt, rlwt_file);
        std::cout<<"Wavelet trees created"<<std::endl;
    }
    simple_rlbwt<sdsl::wt_rlmn<>> other_bwt;
    sdsl::load_from_file(other_bwt, rlwt_file);

    std::string query_pat_file = input_prefix+".pats";
    assert(std::filesystem::exists(query_pat_file));
    //=====

    test_rank(bwt_th_dt,"rlbwt_th_vlbt", other_bwt, "wt_rlmn");
    test_count(bwt_th_dt, "rlbwt_th_vlbt", other_bwt, "wt_rlmn", query_pat_file);
    test_access(bwt_th_dt, "rlbwt_th_vlbt", other_bwt, "wt_rlmn");
    test_inverse_select(bwt_th_dt, "rlbwt_th_vlbt", other_bwt, "wt_rlmn");

    //create the r-index in case it does not exist
    const std::string ri_file = input_prefix+".ri";
    bool fast = false;
    if(!std::filesystem::exists(ri_file)) {
        std::cout<<"Creating the r-index for "<<input_prefix<<" to compare against VLBT "<<std::endl;
        auto idx = ri::r_index(input_prefix, true);

        std::ofstream out(ri_file);
        out.write((char*)&fast,sizeof(fast));
        idx.serialize(out);
        out.close();
        std::cout<<"r-index created"<<std::endl;
    }

    std::ifstream in(ri_file);
    in.read((char*)&fast,sizeof(fast));
    ri::r_index<> ri;
    ri.load(in);
    in.close();

    test_count_with_sa_head(bwt_th_dt, "rlbwt_th_vlbt", ri, "r_index", query_pat_file);
}

template<class size_type>
void test_phi(std::string& input_prefix, size_t ssamp_step, std::string& output_prefix){

    std::string samp_sa_file = input_prefix+".sa_samples";
    std::string str_ranges_file = input_prefix+".str_ranges";

    std::string ssamp_heads_file = output_prefix+".ssamps_heads";
    std::string ssamp_tails_file = output_prefix+".ssamps_tails";

    subsample_sa_samples<size_type>(samp_sa_file, str_ranges_file, ssamp_step, ssamp_heads_file, ssamp_tails_file);

    tmp_workspace tmp_ws("./", true);

    using phi_type = vlbt_phi<NO_VALID_AREA, 65536>;
    phi_type phi;
    phi.subsamp_step = ssamp_step;
    build_phi<phi_type, uint64_t>(phi, ssamp_tails_file, tmp_ws);

    std::string output_file = output_prefix+".vlbt_phi";
    size_t written_bytes = store_to_file(output_file, phi);
    std::cout<<"We store "<<written_bytes<<" in "<<output_file<<std::endl;
}

template<class sa_samp_type>
void test_sr_index(const std::string& input_prefix, BWT_FORMAT bwt_file_fmt, size_t sri_samp_val, const std::string& output_prefix){

    vlbt_sri_va<65536, 65536> sr_index;
    build_sr_index<sa_samp_type>(sr_index, input_prefix, bwt_file_fmt, sri_samp_val);
    std::string output_file = output_prefix+".sri_vlbt";
    size_t written_bytes = store_to_file(output_file, sr_index);
    std::cout<<"We store "<<written_bytes<<" bytes ("<< double(written_bytes*8)/double(sr_index.size())<<" bps) in "<<output_file<<std::endl;

    //=====
    //Create the RLBWT in case it does not exist
    std::string rlwt_file = input_prefix+".wt_rlmn";
    if(!std::filesystem::exists(rlwt_file)){
        std::string bwt_file = input_prefix+".bwt";
        std::cout<<"Creating the RLBWT for "<<bwt_file<<" to compare against VLBT "<<std::endl;
        sdsl::wt_rlmn<> wt;
        sdsl::construct(wt, bwt_file, 1);
        simple_rlbwt simple_rlbwt(wt);
        sdsl::store_to_file(simple_rlbwt, rlwt_file);
        std::cout<<"Wavelet trees created"<<std::endl;
    }
    simple_rlbwt<sdsl::wt_rlmn<>> other_bwt;
    sdsl::load_from_file(other_bwt, rlwt_file);

    std::string query_pat_file = input_prefix+".pats";
    assert(std::filesystem::exists(query_pat_file));
    //=====

    test_rank(sr_index,"rlbwt_sri", other_bwt, "wt_rlmn");
    test_access(sr_index, "rlbwt_sri", other_bwt, "wt_rlmn");
    test_inverse_select(sr_index, "rlbwt_sri", other_bwt, "wt_rlmn");

    //create the r-index in case it does not exist
    const std::string ri_file = input_prefix+".ri";
    bool fast = false;
    if(!std::filesystem::exists(ri_file)) {
        std::cout<<"Creating the r-index for "<<input_prefix<<" to compare against VLBT "<<std::endl;
        std::string text = input_prefix;
        auto idx = ri::r_index(text, true);

        std::ofstream out(ri_file);
        out.write((char*)&fast,sizeof(fast));
        idx.serialize(out);
        out.close();
        std::cout<<"r-index created"<<std::endl;
    }

    std::ifstream in(ri_file);
    in.read((char*)&fast,sizeof(fast));
    ri::r_index<> ri;
    ri.load(in);
    in.close();

    test_count(sr_index, "rlbwt_sri", ri, "r_index", query_pat_file);
    test_count_with_sa_head(sr_index, "rlbwt_sri", ri, "r_index", query_pat_file);
    test_locate(sr_index, "rlbwt_sri", ri, "r_index", query_pat_file);
}

int main(int argc, char** argv) {

    if(argc!=3){
        std::cout<<"usage: ./test_vlbt input_prefix output_prefix"<<std::endl;
        exit(1);
    }

    auto input_text = std::string(argv[1]);
    const auto output_prefix = std::string(argv[2]);

    std::cout<<"Testing RLBWT"<<std::endl;
    test_bwt(input_text, PLAIN, output_prefix);

    //std::cout<<"Testing RLBWT with toeholds"<<std::endl;
    //test_bwt_th<uint64_t>(input_text, PLAIN, 4, output_prefix);

    //std::cout<<"Testing sr-index with valid area"<<std::endl;
    //test_sr_index<uint64_t>(input_text, PLAIN, 4, output_prefix);
    //test_phi<uint64_t>(input_prefix, 4, output_prefix);
}