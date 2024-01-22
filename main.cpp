#include <iostream>
#include "simple_rl_bwt.h"
#include "fm_index.h"

int main() {
    std::string file="/Users/ddiaz/CLionProjects/grlBWT/cmake-build-debug/SRR10971019-part.rl_bwt";
    simple_rl_bwt bwt(file);
    std::string output_file="resulting_bwt";
    size_t written_bytes = store_to_file(output_file, bwt);
    std::cout<<"The rl encoding uses "<<written_bytes<<" bytes "<<std::endl;

    std::cout<<"Now building the FM index"<<std::endl;
    fm_index fmi(file);
    std::cout<<sdsl::size_in_bytes(fmi.bwt)<<std::endl;

    bwt_buff_reader bwt_buff(file);


    std::vector<size_t> tmp_ranks(256, 0);
    size_t pos=0, sym, len;
    unsigned long my_time=0, their_time=0;
    for(size_t i=0;i<bwt_buff.size();i++){
        bwt_buff.read_run(i, sym, len);

        if((i % 100000)==0){
            std::cout<<" -> "<<sym<<" "<<len<<" "<<pos<<std::endl;
        }

        for(size_t k=0;k<len;k++){

            auto t1 = std::chrono::high_resolution_clock::now();
            size_t my_rank =bwt.rank(pos, sym);
            auto t2 = std::chrono::high_resolution_clock::now();
            my_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();


            t1 = std::chrono::high_resolution_clock::now();
            size_t their_rank = fmi.bwt.rank(pos, sym);
            t2 = std::chrono::high_resolution_clock::now();
            their_time += std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

            bool equal = my_rank==their_rank;

            if(!equal){
                std::cout<<"? "<<pos<<","<<sym<<" -> "<<my_rank<<" "<<their_rank<<" "<<tmp_ranks[sym]<<std::endl;
                bwt.rank(pos, sym);
            }
            assert(equal);
            tmp_ranks[sym]++;
            pos++;
        }
    }
    std::cout<<"my average time:    "<<double(my_time)/double(fmi.size())<<" nano seconds "<<std::endl;
    std::cout<<"their average time: "<<double(their_time)/double(fmi.size())<<" nano seconds "<<std::endl;

    return 0;
}
