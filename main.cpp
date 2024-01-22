#include <iostream>
#include "simple_rl_bwt.h"

int main() {
    std::string file="/Users/ddiaz/CLionProjects/grlBWT/cmake-build-debug/SRR10971019-part.rl_bwt";
    simple_rl_bwt bwt(file);
    std::string output_file="resulting_bwt";
    size_t written_bytes = store_to_file(output_file, bwt);
    std::cout<<"The rl encoding uses "<<written_bytes<<" bytes "<<std::endl;
    return 0;
}
