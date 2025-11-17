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

#include "fm_index.h"
#include "simple_rlbwt.h"
#include "custom_wt_rlmn.hpp"
#include "r-index/internal/r_index.hpp"
#include "../scripts/utils.h"
#include <unordered_set>
#include <vector>
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

#define build_dt(dt, suffix) \
{                            \
dt instance;                 \
auto t1 = std::chrono::high_resolution_clock::now();\
sdsl::construct(instance, plain_input_file, 1);\
auto t2 = std::chrono::high_resolution_clock::now();\
sdsl::store_to_file(instance, input_prefix+"."+suffix); \
std::cout<<suffix<<" "<<report_time(t1, t2, 0)<<",  space_usage: "<<double(sdsl::size_in_bytes(instance))/1000000<<" MB ("<<float(sdsl::size_in_bytes(instance)*8)/float(instance.size())<<" bps)"<<std::endl;\
}\

#define TESTED_DTS \
build_dt(sdsl::wt_rlmn<>, "wt_rlmn");\

#define MEASURE_VECTOR(query, time_answer, time_vec, query_answer, time_unit) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query_answer = query;\
auto t2 = std::chrono::high_resolution_clock::now();\
size_t time = std::chrono::duration_cast<time_unit>( t2 - t1 ).count();\
time_answer+=time;\
time_vec.push_back(time);\
}\

#define MEASURE(query, time_answer, query_answer, time_unit) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query_answer = query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<time_unit>( t2 - t1 ).count();\
}

template<class my_bwt_type, class other_bwt_type>
void test_count(my_bwt_type& my_bwt, const std::string& my_dt_name,
                other_bwt_type& other_bwt, std::string other_dt_name, std::string pat_file){

    std::cout<<"Testing count (microsecs/pat and microsecs/occ)"<<std::endl;

    ulint n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    double other_acc_time=0;
    size_t other_acc_count=0;
    size_t j=0;
    std::vector<std::pair<uint64_t, uint64_t>> other_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(other_bwt.backward_search(p), other_acc_time, other_ans[j], std::chrono::microseconds)
        other_acc_count+=other_ans[j].second-other_ans[j].first+1;
        j++;
    }

    double my_acc_time=0;
    size_t my_acc_count=0;
    j=0;
    std::vector<std::pair<uint64_t, uint64_t>> my_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(my_bwt.count(p), my_acc_time, my_ans[j], std::chrono::microseconds)
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

template<class bwt_type>
void test_locate(bwt_type& my_dt, std::string& input_prefix, const std::string& my_dt_name){

    std::cout<<"Testing locate (microsecs/pat and microsecs/occ)"<<std::endl;

    std::string pat_file = input_prefix+".pats";
    std::ifstream ifs(pat_file);
    std::string header;
    std::getline(ifs, header);
    ulint n_pats = get_number_of_patterns(header);
    ulint pat_len = get_patterns_length(header);
    std::cout<<"\t"<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    std::vector<std::string> pat_list(n_pats);
    for(ulint i=0;i<n_pats;++i){
        pat_list[i].reserve(pat_len);
        for(ulint j=0;j<pat_len;++j){
            char c;
            ifs.get(c);
            pat_list[i].push_back(c);
        }
    }

    /*std::string sa_file = input_prefix+".sa";
    size_t f_size = util::file_size(sa_file);
    std::ifstream sa_ifs(sa_file);
    std::vector<uint32_t> sa;
    sa.resize(f_size/sizeof(uint32_t));
    sa_ifs.read((char *)sa.data(), f_size);
    assert(sa.size()==my_dt.size());*/

    double my_acc_time=0;
    //size_t ans;
    //size_t n_valid=0;
    //for (size_t i=0;i<(sa.size()-1);i++) {
    //    MEASURE(my_dt.phi(sa[i]), my_acc_time, ans, std::chrono::nanoseconds)
    //    if(ans==-1) continue;//these are computed differently
    //    if(sa[i+1]!=ans) {
    //        std::cout<<"i:"<<i<<" sa_val:"<<sa[i]<<" -> correct:"<<sa[i+1]<<" / my_answer:"<<ans<<std::endl;
    //    }
    //    assert(sa[i+1]==ans);
    //    n_valid++;
    //}
    //std::cout<<"\t"<<my_dt_name<<": ("<<my_acc_time/double(sa.size()-1)<<", "<<my_acc_time/double(sa.size()-1)<<"), "<<n_valid<<"/"<<(sa.size()-1)<<std::endl;

    //std::string rindex_file = input_prefix+".ri";
    //std::ifstream rindex_ifs(rindex_file);
    //ri::r_index<> rindex;
    //rindex.load(rindex_ifs);

    /*custom_wt_rlmn<> wt_rlmn;
    sdsl::load_from_file(wt_rlmn, input_prefix+".wt_rlmn");
    std::vector<uint64_t> C(wt_rlmn.sigma+1, 0);
    for(size_t sym=0;sym<wt_rlmn.sigma;sym++){
        C[sym] = wt_rlmn.rank(wt_rlmn.size(), my_dt.eff2byte(sym));
    }
    size_t acc=0, tmp;
    for(size_t i=0;i<wt_rlmn.sigma;i++){
        tmp = C[i];
        C[i] = acc;
        acc+=tmp;
    }
    C[wt_rlmn.sigma] = acc;
    std::string samp_sa_file = input_prefix+".sa_samples";
    fm_index<custom_wt_rlmn<>, true> csa_rlmn(wt_rlmn, C, samp_sa_file, my_dt.get_packed_alpha(), my_dt.get_unpacked_alpha());
    size_t acc_count=0;
    size_t j=0;
    double rlmn_acc_time=0;
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> rlmn_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(csa_rlmn.count_with_head(p), rlmn_acc_time, rlmn_ans[j], std::chrono::nanoseconds)
        acc_count+=std::get<1>(rlmn_ans[j])-std::get<0>(rlmn_ans[j])+1;
        assert(std::get<2>(rlmn_ans[j])==sa[std::get<0>(rlmn_ans[j])]);
        j++;
    }*/

    my_acc_time=0;
    size_t my_acc_count=0;
    size_t j=0;
    std::vector<std::vector<uint64_t>> my_ans(n_pats);
    std::vector<double> times;
    for(const std::string& p : pat_list) {
        MEASURE_VECTOR(my_dt.locate(p), my_acc_time, times, my_ans[j], std::chrono::microseconds)
        my_acc_count+=my_ans[j].size();

        /*auto res = my_dt.count(p);
        assert((res.second-res.first+1)==my_ans[j].size());
        for(size_t i=res.first, k=0;i<=res.second;++i,++k){
            if (my_ans[j][k]!=sa[i]) {
                std::cout<<j<<" pattern:"<<p<<" my_ans:"<<my_ans[j][k]<<" != real:"<<sa[i]<<std::endl;
            }
            assert(my_ans[j][k]==sa[i]);
        }*/

        j++;
    }
    print_histogram(times);
    std::cout<<"\t"<<my_dt_name<<": ("<<my_acc_time/double(n_pats)<<", "<<my_acc_time/double(my_acc_count)<<"), tot. occ: "<<my_acc_count<<std::endl;;

    //std::cout<<"find leaf:"<<double(my_dt.phi.acc_time_a)/double(my_dt.phi.acc_time_b+my_dt.phi.acc_time_a)<<" scan leaf:"<<double(my_dt.phi.acc_time_b)/double(my_dt.phi.acc_time_b+my_dt.phi.acc_time_a)<<std::endl;
    //std::cout<<"wt_rlmn: ("<<rlmn_acc_time/double(n_pats)<<", "<<rlmn_acc_time/double(acc_count)<<"), tot. occ: "<<my_acc_count<<std::endl;
    /*for(size_t i=0;i<pat_list.size();i++){
        if(std::get<0>(my_ans[i])!=std::get<0>(rlmn_ans[i]) ||
           std::get<1>(my_ans[i])!=std::get<1>(rlmn_ans[i]) ||
           std::get<2>(my_ans[i])!=std::get<2>(rlmn_ans[i])){
            std::cout<<"Pattern["<<i<<"]: \""<<pat_list[i]<<"\" coords:"<<std::get<0>(my_ans[i])<<"!="<<std::get<0>(rlmn_ans[i])<<" or "
                                                                        <<std::get<1>(my_ans[i])<<"!="<<std::get<1>(rlmn_ans[i])<<" or "
                                                                        <<std::get<2>(my_ans[i])<<"!="<<std::get<2>(rlmn_ans[i])<<std::endl;
            exit(1);
        }
    }*/
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
void test_rank(my_bwt_type& my_dt, std::string my_dt_name, other_dt_type other_dt, std::string other_dt_name){

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

void test_bwt(const std::string& input_prefix, const BWT_FORMAT fmt, const std::string& output_prefix){

    std::cout<<"Testing VLBT BWT"<<std::endl;
    using my_bwt_type = vlbt_bwt<RLBWT, 262144>;
    my_bwt_type bwt_dt;

    std::string input_bwt = input_prefix+".bwt";
    std::cout<<"Building VLBT from input BWT "<<input_bwt<<std::endl;
    build_bwt<my_bwt_type>(bwt_dt, input_bwt, fmt);
    const std::string output_file = output_prefix+".vlbt_bwt";
    const size_t written_bytes = store_to_file(output_file, bwt_dt);
    std::cout<<"We store "<<written_bytes<<" bytes in "<<output_file<<"\n"<<std::endl;
    //std::cout<<"Loading "<<output_file<<std::endl;
    //load_from_file(bwt_dt, output_file);

    //=====
    //Create the RLBWT in case it does not exist
    std::string rlwt_file = input_prefix+".wt_rlmn";
    if(!std::filesystem::exists(rlwt_file)){
        std::cout<<"Creating the RLBWT for "<<input_bwt<<" to compare against VLBT "<<std::endl;
        sdsl::wt_rlmn<> wt;
        sdsl::construct(wt, input_bwt, 1);
        simple_rlbwt simple_rlbwt(wt);
        sdsl::store_to_file(simple_rlbwt, rlwt_file);
        std::cout<<"Wavelet trees created"<<std::endl;
    }
    simple_rlbwt<sdsl::wt_rlmn<>> other_bwt;
    sdsl::load_from_file(other_bwt, rlwt_file);

    std::string query_pat_file = input_prefix+".pats";
    assert(std::filesystem::exists(query_pat_file));
    //=====

    test_rank(bwt_dt,"vlbt_bwt", other_bwt, "wt_rlmn");
    test_count(bwt_dt, "vlbt_bwt", other_bwt, "wt_rlmn", query_pat_file);
    test_access(bwt_dt, "vlbt_bwt", other_bwt, "wt_rlmn");
    test_inverse_select(bwt_dt, "vlbt_bwt", other_bwt, "wt_rlmn");
}

template<class sa_samp_type>
void test_bwt_th(std::string& input_prefix, size_t subsamp_val, std::string& output_prefix){

    std::cout<<"Testing VLBT BWT with toeholds"<<std::endl;
    std::string bwt_file = input_prefix+".ebwt";
    using bwt_th_type = vlbt_bwt<RLBWT_WITH_TOEHOLDS, 65536>;
    bwt_th_type bwt_dt;

    std::string samp_sa_file = input_prefix+".sa_samples";
    std::string str_ranges_file = input_prefix+".str_ranges";
    std::string out_ssamp_heads_file = output_prefix+".ssamp_heads";
    std::string out_ssamp_tails_file = output_prefix+".ssamp_tails";
    subsample_sa_samples<sa_samp_type>(samp_sa_file, str_ranges_file, subsamp_val, out_ssamp_heads_file, out_ssamp_tails_file);

    build_bwt_th<bwt_th_type, sa_samp_type>(bwt_dt, bwt_file, GRL_BWT, out_ssamp_heads_file);
    std::string output_file = output_prefix+".vlbt_bwt";
    size_t written_bytes = store_to_file(output_file, bwt_dt);
    std::cout<<"We store "<<written_bytes<<" in "<<output_file<<"\n"<<std::endl;

    //test_inverse_select(bwt_dt, input_prefix, "vlbt_bwt_th");
    //test_access(bwt_dt, input_prefix, "vlbt_bwt");
    //test_rank(bwt_dt, input_prefix, "vlbt_bwt");
    //test_count(bwt_dt, input_prefix, "vlbt_bwt");
    //test_locate(bwt_dt, input_prefix, "vlbt_bwt");
}

template<class size_type>
void test_phi(std::string& input_prefix, size_t ssamp_step, std::string& output_prefix){

    std::string samp_sa_file = input_prefix+".sa_samples";
    std::string str_ranges_file = input_prefix+".str_ranges";

    std::string ssamp_heads_file = output_prefix+".ssamps_heads";
    std::string ssamp_tails_file = output_prefix+".ssamps_tails";

    subsample_sa_samples<size_type>(samp_sa_file, str_ranges_file, ssamp_step, ssamp_heads_file, ssamp_tails_file);

    using phi_type = vlbt_phi<NO_VALID_AREA, 65536, 64, 4>;
    phi_type phi;
    phi.subsamp_step = ssamp_step;
    build_phi<phi_type, uint64_t>(phi, ssamp_tails_file);

    std::string output_file = output_prefix+".vlbt_phi";
    size_t written_bytes = store_to_file(output_file, phi);
    std::cout<<"We store "<<written_bytes<<" in "<<output_file<<std::endl;
}

template<class size_type>
void test_sr_index(std::string& input_prefix, size_t subsamp_step, std::string& output_prefix){

    using bwt_th_type = vlbt_bwt<RLBWT_WITH_TOEHOLDS, 65536>;
    using phi_type = vlbt_phi<WITH_VALID_AREA, 65536>;
    using sr_index_type = vlbt_sr_index<bwt_th_type, phi_type>;
    sr_index_type sr_index;
    build_sr_index<sr_index_type , size_type>(sr_index, input_prefix, subsamp_step, output_prefix);

    std::string output_sr_index_file = output_prefix+".sr_index";
    size_t written_bytes = store_to_file(output_sr_index_file, sr_index);
    std::cout<<"Final sr-index uses "<<written_bytes<<" bytes ("<< double(written_bytes*8)/double(sr_index.size())<<" bps)"<<std::endl;

    std::string ri_file = input_prefix+".ri";
    if (!std::filesystem::exists(ri_file)) {
        std::cout<<"Creating the r-index for "<<input_prefix<<std::endl;
        auto idx = ri::r_index<>(input_prefix, true);
        std::ofstream out(ri_file);
        idx.serialize(out);
    }

    //test_access(sr_index, input_prefix, "sr_index");
    //test_rank(sr_index, input_prefix, "sr_index");
    //test_inverse_select(sr_index, input_prefix, "sr_index");
    //test_count(sr_index, input_prefix, "sr_index");
    //test_locate(sr_index, input_prefix, "sr_index");
}

int main(int argc, char** argv) {

    if(argc!=3){
        std::cout<<"usage: ./test_vlbt input_prefix output_prefix"<<std::endl;
        exit(1);
    }

    std::string input_text = std::string(argv[1]);
    const auto output_prefix = std::string(argv[2]);
    test_bwt(input_text, PLAIN, output_prefix);

    //test_bwt_th<uint64_t>(input_prefix, 4, output_prefix);
    //test_phi<uint64_t>(input_prefix, 4, output_prefix);
    //test_sr_index<uint64_t>(input_prefix, 16, output_prefix);
}