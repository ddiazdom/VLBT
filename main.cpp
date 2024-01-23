#include <iostream>
#include "simple_rl_bwt.h"
#include "fm_index.h"

void test_interval_symbols(fm_index& fmi, simple_rl_bwt& srlbwt){

    std::vector<uint8_t> my_cs(16, 0);
    std::vector<size_t> my_rank_c_i(16, 0);
    std::vector<size_t> my_rank_c_j(16, 0);
    size_t my_k;

    typedef unsigned long long size_type;
    std::vector<uint8_t> their_cs(16, 0);
    std::vector<size_type> their_rank_c_i(16, 0);
    std::vector<size_type> their_rank_c_j(16, 0);
    size_type their_k;

    //srlbwt.interval_symbols(0, 4097, my_k, my_cs, my_rank_c_i, my_rank_c_j);
    unsigned long my_time=0, their_time=0, n_tries=0;


    for(size_t i=0;i<1;i++) {
        for(size_t j=i+1;j<70894549;j++){
            std::cout<<"We will try : "<<i<<" "<<j<<" "<<fmi.size()<<std::endl;

            /*for(size_t u=i;u<j;u++){
                std::cout<<srlbwt[u]<<","<<fmi.bwt[u]<<"   ";
            }
            std::cout<<""<<std::endl;*/
            auto t1 = std::chrono::high_resolution_clock::now();
            fmi.bwt.interval_symbols(i, j, their_k, their_cs, their_rank_c_i, their_rank_c_j);
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
            srlbwt.interval_symbols(i, j, my_k, my_cs, my_rank_c_i, my_rank_c_j);
            t2 = std::chrono::high_resolution_clock::now();
            my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

            /*std::cout<<"My results: "<<std::endl;
            for(size_t k=0;k<my_k;k++){
                std::cout<<int(my_cs[k])<<" "<<my_rank_c_i[k]<<" "<<my_rank_c_j[k]<<std::endl;
            }
            std::cout<<" "<<std::endl;*/

            for(size_t k=0;k<my_k;k++){
                assert(their_k==my_k);
                assert(their_cs[res_sorted[k].first]==my_cs[k]);
                assert(their_rank_c_i[res_sorted[k].first] == my_rank_c_i[k]);
                assert(their_rank_c_j[res_sorted[k].first] == my_rank_c_j[k]);
            }
            n_tries++;
        }
    }
    std::cout<<"my average time:    "<<double(my_time)/double(fmi.size())<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(fmi.size())<<" nano seconds "<<std::endl;
}

int main() {
    std::string file="/Users/ddiaz/CLionProjects/grlBWT/cmake-build-debug/SRR10971019-part.rl_bwt";
    simple_rl_bwt bwt(file);
    std::string output_file="resulting_bwt";
    size_t written_bytes = store_to_file(output_file, bwt);

    std::cout<<"Now building the FM index"<<std::endl;
    fm_index fmi(file);

    test_interval_symbols(fmi, bwt);

    /*bwt_buff_reader bwt_buff(file);
    std::vector<size_t> tmp_ranks(256, 0);
    size_t pos=0, sym, len;
    unsigned long my_time=0, their_time=0;
    for(size_t i=0;i<bwt_buff.size();i++){
        bwt_buff.read_run(i, sym, len);

        if((i % 100000)==0){
            //std::cout<<" -> "<<sym<<" "<<len<<" "<<pos<<std::endl;
        }

        for(size_t k=0;k<len;k++){

            auto t1 = std::chrono::high_resolution_clock::now();
            //size_t my_res = bwt.rank(pos, sym);
            //size_t my_res = bwt[pos];
            auto my_res = bwt.inverse_select(pos);
            auto t2 = std::chrono::high_resolution_clock::now();
            my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();


            t1 = std::chrono::high_resolution_clock::now();
            //size_t their_rank = fmi.bwt.rank(pos, sym);
            //size_t their_rank = fmi.bwt[pos];
            auto their_res = bwt.inverse_select(pos);
            t2 = std::chrono::high_resolution_clock::now();
            their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

            bool equal = my_res==their_res;

            //std::cout<<"? "<<pos<<","<<sym<<" -> "<<my_res.first<<" "<<int(my_res.second)<<" -> "<<their_res.first<<" "<<int(their_res.second)<<std::endl;
            if(!equal){
                std::cout<<"? "<<pos<<","<<sym<<" -> "<<my_res.first<<" "<<int(my_res.second)<<" -> "<<their_res.first<<" "<<int(their_res.second)<<std::endl;
            }
            assert(equal);
            tmp_ranks[sym]++;
            pos++;
        }
    }

    std::cout<<sdsl::size_in_bytes(fmi.bwt)<<std::endl;
    bwt.stats();
    std::cout<<"The rl encoding uses "<<written_bytes<<" bytes "<<std::endl;
    std::cout<<"my average time:    "<<double(my_time)/double(fmi.size())<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(fmi.size())<<" nano seconds "<<std::endl;*/

    return 0;
}
