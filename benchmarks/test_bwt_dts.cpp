#include <iostream>
#include "../simple_rl_bwt.h"
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/wt_algorithm.hpp"
#include "fb_wt/wt-fbb-0.1.0/wt_fbb.hpp"
#include <vector>

std::vector<std::string> wt_dt = {"my_simple_bwt", "wt_huff_bv", "wt_rlmn", "wt_fbb_hyb"};

#define LOAD(dt, var, file_name) \
dt var;\
sdsl::load_from_file(var, file_name); \


#define MY_LOAD(var, file_name, alphabet) \
simple_rl_bwt<alphabet> var;\
load_from_file(file_name, var); \

#define MEASURE(query, time_answer, query_answer) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query_answer = query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();\
}                                                 \

#define MEASURE_VOID_OUT(query, time_answer) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();\
}\

template<uint8_t alphabet>
void test_inverse_select(std::string& input_file){

    std::vector<double> times(wt_dt.size(), 0);
    std::vector<std::pair<size_t, uint8_t>> answers(wt_dt.size());

    MY_LOAD(srlbwt, input_file+"."+wt_dt[0], alphabet);
    LOAD(sdsl::wt_huff<>, wt_huff, input_file+"."+wt_dt[1])
    LOAD(sdsl::wt_rlmn<>, wt_rlmn, input_file+"."+wt_dt[2])
    LOAD(wt_fbb<sdsl::hyb_vector<>>, wt_fbb_hyb, input_file+"."+wt_dt[3])

    size_t samp_size = (srlbwt.size()*10)/100;

    //for(size_t i=0;i<samp_size;i++){
    for(size_t j=0;j<srlbwt.size();j++){

        //size_t j = rand() % srlbwt.size();

        MEASURE(srlbwt.inverse_select(j), times[0], answers[0]);
        MEASURE(wt_huff.inverse_select(j), times[1], answers[1]);
        MEASURE(wt_rlmn.inverse_select(j), times[2], answers[2]);
        MEASURE(wt_fbb_hyb.inverse_select(j), times[3], answers[3]);

        size_t n_eqs=0;
        for(size_t k=1;k<4;k++){
            n_eqs+=answers[k].first==answers[0].first && answers[k].second==answers[0].second;
        }

        if(n_eqs<(wt_dt.size()-1)){
            std::cout<<"? pos:"<<j<<","<<" -> rank:"<<answers[0].first<<" sym:"<<int(answers[0].second)<<" -> rank:"
                      <<answers[1].first<<" sym:"<<answers[1].second<<std::endl;
        }
        assert(n_eqs==(wt_dt.size()-1));
    }

    std::cout<<"inverse_select ";
    for(size_t i=0;i<wt_dt.size();i++){
        std::cout<<times[i]/double(samp_size)<<" ";
    }
    std::cout<<""<<std::endl;
}

template<uint8_t alphabet>
void test_access(std::string& input_file){

    std::vector<double> times(wt_dt.size(), 0);
    std::vector<uint8_t> answers(wt_dt.size());

    MY_LOAD(srlbwt, input_file+"."+wt_dt[0], alphabet);
    LOAD(sdsl::wt_huff<>, wt_huff, input_file+"."+wt_dt[1])
    LOAD(sdsl::wt_rlmn<>, wt_rlmn, input_file+"."+wt_dt[2])
    LOAD(wt_fbb<sdsl::hyb_vector<>>, wt_fbb_hyb, input_file+"."+wt_dt[3])

    size_t samp_size = (srlbwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){
    //for(size_t j=0;j<srlbwt.size();j++){

        size_t j = rand() % srlbwt.size();

        MEASURE(srlbwt[j], times[0], answers[0]);
        MEASURE(wt_huff[j], times[1], answers[1]);
        MEASURE(wt_rlmn[j], times[2], answers[2]);
        MEASURE(wt_fbb_hyb[j], times[3], answers[3]);

        size_t n_eqs=0;
        for(size_t k=1;k<4;k++){
            n_eqs+=answers[k]==answers[0];
        }

        if(n_eqs<(wt_dt.size()-1)){
            std::cout<<"? pos:"<<j<<","<<" -> sym:"<<int(answers[0])<<" -> sym:"<<int(answers[1])<<std::endl;
        }
        assert(n_eqs==(wt_dt.size()-1));
    }

    std::cout<<"access ";
    for(size_t i=0;i<wt_dt.size();i++){
        std::cout<<times[i]/double(samp_size)<<" ";
    }
    std::cout<<""<<std::endl;
}

template<uint8_t alphabet>
void test_rank(std::string& input_file) {

    std::vector<double> times(wt_dt.size(), 0);
    std::vector<size_t> answers(wt_dt.size(), 0);

    MY_LOAD(srlbwt, input_file+"."+wt_dt[0], alphabet);
    LOAD(sdsl::wt_huff<>, wt_huff, input_file+"."+wt_dt[1])
    LOAD(sdsl::wt_rlmn<>, wt_rlmn, input_file+"."+wt_dt[2])
    LOAD(wt_fbb<sdsl::hyb_vector<>>, wt_fbb_hyb, input_file+"."+wt_dt[3])

    size_t samp_size = (srlbwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){
    //for(size_t u=0;u<samp_size;u++){

        size_t u = rand() % srlbwt.size();
        uint8_t s = rand() % srlbwt.alphabet;

        MEASURE(srlbwt.rank(u, srlbwt.sym_inv_map[s]), times[0], answers[0]);
        MEASURE(wt_huff.rank(u, srlbwt.sym_inv_map[s]), times[1], answers[1]);
        MEASURE(wt_rlmn.rank(u, srlbwt.sym_inv_map[s]), times[2], answers[2]);
        MEASURE(wt_fbb_hyb.rank(u, srlbwt.sym_inv_map[s]), times[3], answers[3]);

        size_t n_eqs=0;
        for(size_t k=1;k<4;k++){
            n_eqs+=answers[k]==answers[0];
        }

        if(n_eqs<(wt_dt.size()-1)){
            std::cout<<"? position:"<<u<<", symbol:"<<int(srlbwt.sym_inv_map[s])<<std::endl;
            std::cout<<answers[0]<<" -> "<<answers[1]<<" "<<answers[2]<<" "<<answers[3]<<std::endl;
        }
        assert(n_eqs==(wt_dt.size()-1));
    }

    std::cout<<"rank ";
    for(size_t i=0;i<wt_dt.size();i++){
        std::cout<<times[i]/double(samp_size)<<" ";
    }
    std::cout<<""<<std::endl;
}

template<uint8_t alphabet>
void test_select(std::string& input_file) {

    std::vector<double> times(wt_dt.size(), 0);
    std::vector<size_t> answers(wt_dt.size(), 0);

    MY_LOAD(srlbwt, input_file+"."+wt_dt[0], alphabet);
    LOAD(sdsl::wt_huff<>, wt_huff, input_file+"."+wt_dt[1])
    LOAD(sdsl::wt_rlmn<>, wt_rlmn, input_file+"."+wt_dt[2])
    LOAD(wt_fbb<sdsl::hyb_vector<>>, wt_fbb_il, input_file+"."+wt_dt[3])

    //get the max rank of each symbol to avoid failed asserts
    size_t rank_answers[16] ={0};
    for(size_t i=0;i<srlbwt.alphabet;i++){
        rank_answers[i] = wt_huff.rank(wt_huff.size(), srlbwt.sym_inv_map[i]);
    }


    size_t samp_size = (srlbwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        uint8_t s = rand() % srlbwt.alphabet;
        size_t u = 1+(rand() % rank_answers[s]);

        MEASURE(srlbwt.select(u, srlbwt.sym_inv_map[s]), times[0], answers[0]);
        MEASURE(wt_huff.select(u, srlbwt.sym_inv_map[s]), times[1], answers[1]);
        MEASURE(wt_rlmn.select(u, srlbwt.sym_inv_map[s]), times[2], answers[2]);
        //MEASURE(wt_fbb_il.select(u, srlbwt.sym_inv_map[s]), times[3], answers[3]);

        size_t n_eqs=0;
        for(size_t k=1;k<4;k++){
            n_eqs+=answers[k]==answers[0];
        }

        if(n_eqs<(wt_dt.size()-2)){
            std::cout<<"? position:"<<u<<", symbol:"<<int(srlbwt.sym_inv_map[s])<<std::endl;
            std::cout<<answers[0]<<" -> "<<answers[1]<<" "<<answers[2]<<" "<<std::endl;
        }
        assert(n_eqs==(wt_dt.size()-2));
    }

    std::cout<<"select ";
    for(size_t i=0;i<wt_dt.size();i++){
        if(wt_dt[i]=="wt_fbb_hyb"){
            std::cout<<"- ";
        }else{
            std::cout<<times[i]/double(samp_size)<<" ";
        }
    }
    std::cout<<""<<std::endl;
}

template<class size_type>
struct int_symbol_data{
    std::vector<uint8_t> cs;
    std::vector<size_type> rank_c_i;
    std::vector<size_type> rank_c_j;
    size_type k=0;

    int_symbol_data(): cs(16, 0), rank_c_i(16, 0), rank_c_j(16, 0){}

    bool operator==(int_symbol_data& other) const{
        return k==other.k && rank_c_i==other.rank_c_i && rank_c_j==other.rank_c_j && cs==other.cs;
    }

    void lex_sort(){
        std::vector<std::tuple<size_t, size_t, size_t>> res_sorted(k);
        for(size_t u=0;u<k;u++){
            res_sorted[u] = {cs[u], rank_c_i[u], rank_c_j[u]};
        }

        std::sort(res_sorted.begin(), res_sorted.end(), [](auto &a, auto &b){
            return std::get<0>(a)<std::get<0>(b);
        });

        for(size_t u=0;u<k;u++){
            cs[u] = std::get<0>(res_sorted[u]);
            rank_c_i[u] = std::get<1>(res_sorted[u]);
            rank_c_j[u] = std::get<2>(res_sorted[u]);
        }
    }

    void print(){
        std::cout<<k<<": "<<std::endl;
        for(size_t u=0;u<k;u++){
            std::cout<<"sym: "<<int(cs[u])<<" before:"<<rank_c_i[u]<<" after:"<<rank_c_j[u]<<std::endl;
        }
    }
};

template<uint8_t alphabet>
void test_interval_symbols(std::string& input_file){

    std::vector<double> times(wt_dt.size(), 0);

    typedef sdsl::wt_huff<>::size_type size_type;
    int_symbol_data<size_type> their_answer;
    int_symbol_data<uint64_t> my_answer;

    MY_LOAD(srlbwt, input_file+"."+wt_dt[0], alphabet);
    LOAD(sdsl::wt_huff<>, wt_huff, input_file+"."+wt_dt[1])

    size_t samp_size = (srlbwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t a = rand() % (srlbwt.size()-1);
        size_t b = rand() % (srlbwt.size()-1);

        if(a>b) std::swap(a, b);
        if(a==b) b++;
        assert(a<srlbwt.size() && b<=srlbwt.size());

        MEASURE_VOID_OUT(srlbwt.interval_symbols(a, b, my_answer.k, my_answer.cs, my_answer.rank_c_i, my_answer.rank_c_j), times[0])
        MEASURE_VOID_OUT(wt_huff.interval_symbols(a, b, their_answer.k, their_answer.cs, their_answer.rank_c_i, their_answer.rank_c_j), times[1])

        their_answer.lex_sort();

        assert(my_answer.k==their_answer.k);
        for(size_t k=0;k<my_answer.k;k++){
            if(my_answer.cs[k]!=their_answer.cs[k] ||
               my_answer.rank_c_i[k]!=their_answer.rank_c_i[k] ||
               my_answer.rank_c_j[k]!=their_answer.rank_c_j[k]){
                std::cout<<"range : "<<a<<", "<<b<<std::endl;

                std::cout<<"My results: "<<std::endl;
                my_answer.print();
                std::cout<<" "<<std::endl;
                std::cout<<"their results: "<<std::endl;
                their_answer.print();
                exit(1);
            }
        }

    }
    std::cout<<"interval_symbols ";
    for(size_t i=0;i<wt_dt.size();i++){
        if(i<2){
            std::cout<<times[i]/double(samp_size)<<" ";
        }else{
            std::cout<<"- ";
        }
    }
    std::cout<<""<<std::endl;
}

template<uint8_t alphabet>
void run_measurements(std::string& input_file){
    test_inverse_select<alphabet>(input_file);
    test_select<alphabet>(input_file);
    test_interval_symbols<alphabet>(input_file);
    test_access<alphabet>(input_file);
    test_rank<alphabet>(input_file);
}

int main(int argc,  char** argv) {

    if(argc!=3){
        std::cout<<"usage: ./test_bwt_dts prefix alphabet"<<std::endl;
        exit(1);
    }

    std::string input_file = std::string(argv[1]);
    char *pend;
    long int alphabet = strtol(argv[2], &pend, 10);
    assert(alphabet>2 && alphabet<=16);

    std::cout<<"query ";
    for(const auto & dt : wt_dt){
        std::cout<<dt<<" ";
    }
    std::cout<<""<<std::endl;

    if(alphabet==3){
        run_measurements<3>(input_file);
    }else if(alphabet==4){
        run_measurements<4>(input_file);
    }else if(alphabet==5){
        run_measurements<5>(input_file);
    }else if(alphabet==6){
        run_measurements<6>(input_file);
    }else if(alphabet==7){
        run_measurements<7>(input_file);
    }else if(alphabet==8){
        run_measurements<8>(input_file);
    }else if(alphabet==9){
        run_measurements<9>(input_file);
    }else if(alphabet==10){
        run_measurements<10>(input_file);
    }else if(alphabet==11){
        run_measurements<11>(input_file);
    }else if(alphabet==12){
        run_measurements<12>(input_file);
    }else if(alphabet==13){
        run_measurements<13>(input_file);
    }else if(alphabet==14){
        run_measurements<14>(input_file);
    }else if(alphabet==15){
        run_measurements<15>(input_file);
    }else if(alphabet==16){
        run_measurements<16>(input_file);
    }

    return 0;
}
