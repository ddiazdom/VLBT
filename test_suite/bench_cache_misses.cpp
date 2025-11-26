//
// Created by Diaz, Diego on 25.11.2025.
//
#include <iostream>
#include <ostream>
#include <algorithm>

#include "../include/vlbt_bwt.h"
#include "../include/vlbt_sr_index.h"

#include <unordered_set>
#include <random>

//this function discards any data in the LD1 and LD2 caches
void flush_cache() {
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<uint8_t> dist_sym(0, std::numeric_limits<uint8_t>::max());
    static constexpr size_t SIZE = 64*1024*1024; // 64 MB
    volatile uint8_t* buf = new uint8_t[SIZE];
    for (size_t i = 0; i < SIZE; i++) buf[i] = dist_sym(rng);
    delete[] buf;
}

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

    for (const double t : filtered) {
        const int bin = std::min(int((t - min_val) / bin_width), num_bins - 1);
        bins[bin]++;
    }

    // Print ASCII histogram
    std::cout << "Histogram of runtimes (microseconds):\n";
    for (int i = 0; i < num_bins; ++i) {
        const double bin_start = min_val + i * bin_width;
        const double bin_end = bin_start + bin_width;
        std::cout << "[" << bin_start << ", " << bin_end << "): ";

        const int count = bins[i];
        for (int j = 0; j < static_cast<int>(count * 50 / filtered.size()); ++j) { // scale to max width 50
            std::cout << '#';
        }
        std::cout << " (" << count << ")\n";
    }
}

std::vector<uint64_t> sample_unique(const uint64_t n, uint64_t x) {
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

template<class dt_type>
std::vector<std::pair<uint64_t, uint8_t>> compute_random_rank_queries(const dt_type& dt, uint64_t n_samps) {

    const size_t text_size = dt.size();
    const size_t alphabet_size = dt.alphabet_size();

    if(n_samps > text_size*alphabet_size) throw std::invalid_argument("invalid sample size");

    std::unordered_set<std::pair<uint64_t, uint8_t>> seen;
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist_idx(0, text_size - 1);
    std::uniform_int_distribution<uint64_t> dist_sym(0, alphabet_size- 1);

    while (seen.size() < n_samps) {
        uint64_t idx = dist_idx(rng);
        uint8_t symbol = dt.eff2byte(dist_sym(rng));
        if (dt.rank(idx, symbol)<0)  continue;
        seen.insert({idx, symbol});
    }
    return {seen.begin(), seen.end()};
}


/*template<class my_bwt_type, class other_bwt_type>
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

    size_t my_acc_time=0;
    size_t my_acc_count=0;
    size_t j=0;
    std::vector<std::pair<uint64_t, uint64_t>> my_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(my_dt.count(p), my_acc_time, my_ans[j], std::chrono::nanoseconds)
        my_acc_count+=my_ans[j].second-my_ans[j].first+1;
        j++;
    }

    size_t other_acc_time=0;
    size_t other_acc_count=0;
    j=0;
    std::vector<std::pair<uint64_t, uint64_t>> other_ans(n_pats);
    for(auto const& p : pat_list) {
        MEASURE(other_dt.count(p), other_acc_time, other_ans[j], std::chrono::nanoseconds)
        other_acc_count+=other_ans[j].second-other_ans[j].first+1;
        j++;
    }

    const double my_ns_per_pat = double(my_acc_time)/double(n_pats);
    const double my_ns_per_occ = double(my_acc_time)/double(my_acc_count);

    const double other_ns_per_pat = double(other_acc_time)/double(n_pats);
    const double other_ns_per_occ = double(other_acc_time)/double(my_acc_count);
    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"\t"<<my_dt_name<<": ("<<my_ns_per_pat<<", "<<my_ns_per_occ<<"), ";
    std::cout<<"\t"<<other_dt_name<<": ("<<other_ns_per_pat<<", "<<other_ns_per_occ<<")"<<std::endl;;
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
}*/

//isolate the function to count the number of LD1 and LD2 cache misses
template<class dt_type>
__attribute__((noinline)) void bench_iso_rank(dt_type& dt, const std::vector<std::pair<uint64_t, uint8_t>>& queries, size_t& dummy) {
    for(const auto &[pos ,sym] : queries) {
        dummy += dt.rank(pos, sym);
    }
}
template<class dt_type>
void benchmark_rank(dt_type &dt, const size_t n) {

    std::vector<std::pair<uint64_t, uint8_t>> queries = compute_random_rank_queries(dt, n);
    size_t dummy = 0;
    //NOTE: no need to warmup because the computation of random queries calls ``dt.rank'' to check the answer is >0
    //pollute the cache so the dt starts cold
    flush_cache();

    //perform the benchmarks
    bench_iso_rank(dt, queries, dummy);
    std::cout<<"rank dummy: "<<dummy<<std::endl;//print it to avoid optimizations
}

//isolate the function to count the number of LD1 and LD2 cache misses
template<class dt_type>
__attribute__((noinline)) void bench_iso_access(dt_type& dt, const std::vector<uint64_t>& queries, size_t& dummy) {
    for(unsigned long long query : queries) {
        dummy += dt[query];
    }
}
template<class dt_type>
void benchmark_access(dt_type &dt, const size_t n) {

    //fill in with random positions
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist(0, n-1);
    std::vector<uint64_t> query_pos(n);
    size_t dummy = 0;
    for (int k = 0; k < n; k++) {
        query_pos[k] = dist(rng);
    }
    //

    //warmup
    size_t warmup = std::min<size_t>(5000, n);
    for (int k = 0; k < warmup; k++) {
        dummy+=dt[dist(rng)];
    }
    //

    //pollute the cache so the dt starts cold
    flush_cache();
    bench_iso_access(dt, query_pos, dummy);
    std::cout<<"access dummy: "<<dummy<<std::endl;//print it to avoid optimizations
}

template<class dt_type>
void bench_int2(const std::string& input_index, const size_t n_samp) {
    dt_type bwt_th_dt;
    load_from_file(input_index, bwt_th_dt);
    benchmark_access(bwt_th_dt, std::min<size_t>(n_samp, bwt_th_dt.size()));
    benchmark_rank(bwt_th_dt, std::min<size_t>(n_samp, bwt_th_dt.size()));
}

template<size_t b_size>
void bench_int(const std::string& input_index, const size_t n_samp, VLBT_TYPE& tag) {

    std::cout<<"Using "<<n_samp<<" samples"<<std::endl;

    switch (tag) {
        case RLBWT:
            std::cout<<"Testing VLBT RLBWT with block size "<<b_size<<std::endl;
            bench_int2<vlbt_rlbwt<b_size>>(input_index, n_samp);
            break;
        case RLBWT_WITH_TOEHOLDS:
            std::cout<<"Testing VLBT RLBWT with toeholds and block size "<<b_size<<std::endl;
            bench_int2<vlbt_rlbwt_th<b_size>>(input_index, n_samp);
            break;
        case SRI_VALID_AREA:
            std::cout<<"Testing VLBT sr-index with valid area and block size "<<b_size<<std::endl;
            bench_int2<vlbt_sri_va<b_size, b_size>>(input_index, n_samp);
            break;
        default:
            std::cerr<<"Unknown index_type"<<std::endl;
    }
}

int main(int argc, char** argv) {

    if(argc!=3){
        std::cout<<"usage: ./bench_cmiss_vlbt input_dt n_samples"<<std::endl;
        exit(1);
    }

    const auto input_index = std::string(argv[1]);
    char *pend;
    long int n_samp = strtol(argv[2], &pend, 10);

    temp_param_t tp = read_template_param(input_index);
    switch (tp.b_size) {
        case 1024:
            bench_int<1024>(input_index, n_samp, tp.tag);
            break;
        case 4096:
            bench_int<4096>(input_index, n_samp, tp.tag);
            break;
        case 16384:
            bench_int<16384>(input_index, n_samp, tp.tag);
            break;
        case 65536:
            bench_int<65536>(input_index, n_samp, tp.tag);
            break;
        case 262144:
            bench_int<262144>(input_index, n_samp, tp.tag);
            break;
        case 1048576:
            bench_int<1048576>(input_index, n_samp, tp.tag);
            break;
        default:
            exit(1);
    }
}