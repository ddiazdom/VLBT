#include <iostream>
#include "simple_rl_bwt.h"
#include "fm_index.h"
#include "sdsl/wavelet_trees.hpp"

void test_inverse_select(fm_index& fmi, simple_rl_bwt& srlbwt){
    unsigned long my_time=0, their_time=0;
    std::cout<<"Testing inverse select"<<std::endl;
    size_t samp_size = (fmi.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t j = rand() % srlbwt.size();

        auto t1 = std::chrono::high_resolution_clock::now();
        auto my_res = srlbwt.inverse_select(j);
        auto t2 = std::chrono::high_resolution_clock::now();
        my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        t1 = std::chrono::high_resolution_clock::now();
        auto their_res = fmi.bwt.inverse_select(j);
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

void test_access(fm_index& fmi, simple_rl_bwt& srlbwt){
    unsigned long my_time=0, their_time=0;
    std::cout<<"Testing access"<<std::endl;
    size_t samp_size = (fmi.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t j = rand() % srlbwt.size();

        auto t1 = std::chrono::high_resolution_clock::now();
        auto my_res = srlbwt[j];
        auto t2 = std::chrono::high_resolution_clock::now();
        my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

        t1 = std::chrono::high_resolution_clock::now();
        auto their_res = fmi.bwt[j];
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

void test_rank(fm_index& fmi, simple_rl_bwt& srlbwt){

    std::cout<<"Testing rank"<<std::endl;

    unsigned long my_time=0, their_time=0, n_tries=0;
    size_t samp_size = (fmi.size()*10)/100;

    for(size_t i=0;i<samp_size;i++){

        size_t u = rand() % srlbwt.size();

        for(size_t j=0;j<srlbwt.sym_inv_map.size(); j++){
            auto t1 = std::chrono::high_resolution_clock::now();
            auto my_res = srlbwt.rank(u,srlbwt.sym_inv_map[j]);
            auto t2 = std::chrono::high_resolution_clock::now();
            my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

            t1 = std::chrono::high_resolution_clock::now();
            auto their_res = fmi.bwt.rank(u, srlbwt.sym_inv_map[j]);
            t2 = std::chrono::high_resolution_clock::now();
            their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

            bool equal = my_res==their_res;
            if(!equal){
                std::cout<<"? position:"<<i<<", symbol:"<<srlbwt.sym_inv_map[j]<<" -> "<<my_res<<" -> "<<their_res<<std::endl;
            }
            assert(equal);
            n_tries++;
        }
    }
    std::cout<<"my average time:    "<<double(my_time)/double(n_tries)<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(n_tries)<<" nano seconds "<<std::endl;
}

void test_interval_symbols(fm_index& fmi, simple_rl_bwt& srlbwt){

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

    unsigned long my_time=0, their_time=0, n_tries=0;

    for(size_t i=0;i<60000;i++) {
        size_t a = rand() % fmi.size()-1;
        for(size_t j=i+1;j<i+120;j++){
            //std::cout<<"We will try : "<<i<<" "<<j<<" "<<fmi.size()<<std::endl;

            /*for(size_t u=i;u<j;u++){
                std::cout<<srlbwt[u]<<","<<fmi.bwt[u]<<"   ";
            }
            std::cout<<""<<std::endl;*/
            size_t b = rand() % fmi.size()-1;

            if(a>b){
                std::swap(a, b);
            }

            if(a==b){
                b++;
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            fmi.bwt.interval_symbols(a, b, their_k, their_cs, their_rank_c_i, their_rank_c_j);
            auto t2 = std::chrono::high_resolution_clock::now();
            their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

            std::vector<std::pair<size_t, size_t>> res_sorted(their_k);
            for(size_t k=0;k<their_k;k++){
                res_sorted[k] = {k, their_cs[k]};
            }
            std::sort(res_sorted.begin(), res_sorted.end(), [](auto &a, auto &b){
                return a.second<b.second;
            });

            /*std::cout<<"Their results: "<<std::endl;
            for(size_t k=0;k<their_k;k++){
                std::cout<<int(their_cs[res_sorted[k].first])<<" "<<their_rank_c_i[res_sorted[k].first]<<" "<<their_rank_c_j[res_sorted[k].first]<<std::endl;
            }*/

            t1 = std::chrono::high_resolution_clock::now();
            srlbwt.interval_symbols(a, b, my_k, my_cs, my_rank_c_i, my_rank_c_j);
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
            n_tries++;
        }
    }
    std::cout<<"my average time:    "<<double(my_time)/double(n_tries)<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(n_tries)<<" nano seconds "<<std::endl;
}

int main() {
    std::string file="/Users/ddiaz/CLionProjects/ryu/cmake-build-debug/ryu.idx.cbAMY2/rl_bwt_WQA";

    std::cout<<"Building my rl BWT"<<std::endl;
    simple_rl_bwt bwt(file);
    std::string output_file="resulting_bwt";
    size_t written_bytes = store_to_file(output_file, bwt);
    std::cout<<"It uses "<<written_bytes<<" bytes "<<std::endl;
    bwt.stats();

    std::cout<<"Now building the standard FM index"<<std::endl;
    fm_index fmi(file);
    std::cout<<"It uses "<<sdsl::size_in_mega_bytes(fmi.bwt)<<" MB "<<std::endl;


    //test_access(fmi, bwt);
    //test_rank(fmi, bwt);
    //test_inverse_select(fmi, bwt);
    test_interval_symbols(fmi, bwt);
    //test_select(fmi, bwt);

    return 0;
}
