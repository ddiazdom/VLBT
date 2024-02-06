#include <iostream>
#include "simple_rl_bwt.h"
#include "fm_index.h"
#include "sdsl/wavelet_trees.hpp"

template<class fb_bwt, class competitor_bwt>
void test_inverse_select(competitor_bwt& their_bwt, fb_bwt& my_bwt){
    unsigned long my_time=0, their_time=0;
    std::cout<<"Testing inverse select"<<std::endl;
    size_t samp_size = (their_bwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t j = rand() % their_bwt.size();

        auto t1 = std::chrono::high_resolution_clock::now();
        auto my_res = my_bwt.inverse_select(j);
        auto t2 = std::chrono::high_resolution_clock::now();
        my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        t1 = std::chrono::high_resolution_clock::now();
        auto their_res = their_bwt.inverse_select(j);
        t2 = std::chrono::high_resolution_clock::now();
        their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        bool equal = my_res.first==their_res.first && my_res.second==their_res.second;

        if(!equal){
            std::cout<<"? pos:"<<i<<","<<" -> rank:"<<my_res.first<<" sym:"<<int(my_res.second)<<" -> rank:"<<their_res.first<<" sym:"<<int(their_res.second)<<std::endl;
        }
        assert(equal);
    }
    std::cout<<"my average time:    "<<double(my_time)/double(samp_size)<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(samp_size)<<" nano seconds "<<std::endl;
}

template<class fb_bwt, class competitor_bwt>
void test_access(competitor_bwt& their_bwt, fb_bwt& my_bwt){
    unsigned long my_time=0, their_time=0;
    std::cout<<"Testing access"<<std::endl;
    size_t samp_size = (their_bwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t j = rand() % their_bwt.size();

        auto t1 = std::chrono::high_resolution_clock::now();
        auto my_res = my_bwt[j];
        auto t2 = std::chrono::high_resolution_clock::now();
        my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        t1 = std::chrono::high_resolution_clock::now();
        auto their_res = their_bwt[j];
        t2 = std::chrono::high_resolution_clock::now();
        their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        bool equal = my_res==their_res;

        if(!equal){
            std::cout<<"? pos:"<<i<<","<<" -> sym:"<<int(my_res)<<" -> sym:"<<int(their_res)<<std::endl;
        }
        assert(equal);
    }
    std::cout<<"my average time:    "<<double(my_time)/double(samp_size)<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(samp_size)<<" nano seconds "<<std::endl;
}

template<class fb_bwt, class competitor_bwt>
void test_rank(competitor_bwt& their_bwt, fb_bwt& my_bwt){

    std::cout<<"Testing rank"<<std::endl;

    unsigned long my_time=0, their_time=0;
    size_t samp_size = (their_bwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t u = rand() % their_bwt.size();
        uint8_t s = rand() % their_bwt.sigma;

        auto t1 = std::chrono::high_resolution_clock::now();
        auto my_res = my_bwt.rank(u, my_bwt.sym_inv_map[s]);
        auto t2 = std::chrono::high_resolution_clock::now();
        my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        t1 = std::chrono::high_resolution_clock::now();
        auto their_res = their_bwt.rank(u, my_bwt.sym_inv_map[s]);
        t2 = std::chrono::high_resolution_clock::now();
        their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        bool equal = my_res==their_res;
        if(!equal){
            std::cout<<"? position:"<<i<<", symbol:"<<my_bwt.sym_inv_map[s]<<" -> "<<my_res<<" -> "<<their_res<<std::endl;
        }
        assert(equal);
    }
    std::cout<<"my average time:    "<<double(my_time)/double(samp_size)<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(samp_size)<<" nano seconds "<<std::endl;
}

template<class fb_bwt, class competitor_bwt>
void test_select(competitor_bwt& their_bwt, fb_bwt& my_bwt){

    std::cout<<"Testing select"<<std::endl;

    //get the max rank of each symbol to avoid failed asserts
    size_t rank_answers[16] ={0};
    for(size_t i=0;i<their_bwt.sigma;i++){
        rank_answers[i] = their_bwt.rank(their_bwt.size(), my_bwt.sym_inv_map[i]);
    }

    unsigned long my_time=0, their_time=0;
    size_t samp_size = (their_bwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t s = rand() % their_bwt.sigma;
        size_t rank = 1+(rand() % rank_answers[s]);

        auto t1 = std::chrono::high_resolution_clock::now();
        auto my_res = my_bwt.select(rank, my_bwt.sym_inv_map[s]);
        auto t2 = std::chrono::high_resolution_clock::now();
        my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        t1 = std::chrono::high_resolution_clock::now();
        auto their_res = their_bwt.select(rank, my_bwt.sym_inv_map[s]);
        t2 = std::chrono::high_resolution_clock::now();
        their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        bool equal = my_res==their_res;
        if(!equal){
            std::cout<<"? symbol: "<<my_bwt.sym_inv_map[s]<<" and rank: "<<rank<<" -> "<<my_res<<" -> "<<their_res<<std::endl;
        }
        assert(equal);
    }

    std::cout<<"my average time:    "<<double(my_time)/double(samp_size)<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(samp_size)<<" nano seconds "<<std::endl;
}

template<class fb_bwt, class competitor_bwt>
void test_interval_symbols(competitor_bwt& their_bwt, fb_bwt& my_bwt){

    std::cout<<"Testing interval symbols "<<std::endl;

    std::vector<uint8_t> my_cs(16, 0);
    std::vector<size_t> my_rank_c_i(16, 0);
    std::vector<size_t> my_rank_c_j(16, 0);
    size_t my_k;

    typedef sdsl::wt_huff<>::size_type size_type;
    std::vector<uint8_t> their_cs(16, 0);
    std::vector<size_type> their_rank_c_i(16, 0);
    std::vector<size_type> their_rank_c_j(16, 0);
    size_type their_k;

    unsigned long my_time=0, their_time=0;
    size_t samp_size = (their_bwt.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t a = rand() % their_bwt.size()-1;
        size_t b = rand() % their_bwt.size()-1;

        if(a>b) std::swap(a, b);
        if(a==b) b++;

        auto t1 = std::chrono::high_resolution_clock::now();
        their_bwt.interval_symbols(a, b, their_k, their_cs, their_rank_c_i, their_rank_c_j);
        auto t2 = std::chrono::high_resolution_clock::now();
        their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
        std::vector<std::pair<size_t, size_t>> res_sorted(their_k);
        for(size_t k=0;k<their_k;k++){
            res_sorted[k] = {k, their_cs[k]};
        }
        std::sort(res_sorted.begin(), res_sorted.end(), [](auto &a, auto &b){
            return a.second<b.second;
        });

        t1 = std::chrono::high_resolution_clock::now();
        my_bwt.interval_symbols(a, b, my_k, my_cs, my_rank_c_i, my_rank_c_j);
        t2 = std::chrono::high_resolution_clock::now();
        my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        if(their_k!=my_k){
            std::cout<<"My results: "<<std::endl;
            for(size_t k=0;k<my_k;k++){
                std::cout<<int(my_cs[k])<<" "<<my_rank_c_i[k]<<" "<<my_rank_c_j[k]<<std::endl;
            }
            std::cout<<" "<<std::endl;
            std::cout<<their_k<<" "<<my_k<<" "<<a<<" "<<b<<std::endl;
        }

        assert(their_k==my_k);
        for(size_t k=0;k<my_k;k++){
            if(their_cs[res_sorted[k].first]!=my_cs[k] ||
               their_rank_c_i[res_sorted[k].first]!=my_rank_c_i[k] ||
               their_rank_c_j[res_sorted[k].first]!=my_rank_c_j[k]){

                std::cout<<"My results: "<<std::endl;
                for(size_t p=0;p<my_k;p++){
                    std::cout<<int(my_cs[p])<<" "<<my_rank_c_i[p]<<" "<<my_rank_c_j[p]<<std::endl;
                }
                std::cout<<" "<<std::endl;
                std::cout<<"their results: "<<std::endl;
                for(size_t p=0;p<my_k;p++){
                    std::cout<<int(their_cs[res_sorted[p].first])<<" "<<their_rank_c_i[res_sorted[p].first]<<" "<<their_rank_c_j[res_sorted[p].first]<<std::endl;
                }
                std::cout<<" "<<std::endl;
                std::cout<<their_k<<" "<<my_k<<" "<<a<<" "<<b<<std::endl;
            }
            assert(their_cs[res_sorted[k].first]==my_cs[k]);
            assert(their_rank_c_i[res_sorted[k].first] == my_rank_c_i[k]);
            assert(their_rank_c_j[res_sorted[k].first] == my_rank_c_j[k]);
        }
    }
    std::cout<<"my average time:    "<<double(my_time)/double(samp_size)<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(samp_size)<<" nano seconds "<<std::endl;
}

int main() {
    //std::string file="/Users/ddiaz/CLionProjects/ryu/cmake-build-debug/ryu.idx.cbAMY2/rl_bwt_WQA";
    std::string file="/Users/ddiaz/CLionProjects/grlBWT/cmake-build-debug/SRR10971019-part.rl_bwt";

    std::cout<<"Building my rl BWT"<<std::endl;
    simple_rl_bwt<5> bwt(file);
    std::string output_file="resulting_bwt";
    size_t written_bytes = store_to_file(output_file, bwt);
    std::cout<<"It uses "<<written_bytes<<" bytes "<<std::endl;
    bwt.stats();

    std::cout<<"Now building the standard FM index"<<std::endl;
    fm_index fmi(file);
    std::cout<<"It uses "<<sdsl::size_in_mega_bytes(fmi.bwt)<<" MB "<<std::endl;

    test_access(fmi.bwt, bwt);
    test_rank(fmi.bwt, bwt);
    test_inverse_select(fmi.bwt, bwt);
    test_interval_symbols(fmi.bwt, bwt);
    test_select(fmi.bwt, bwt);

    return 0;
}
