#include "../include/bwt_io.h"
#include <filesystem>

int main(int argc, char** argv){

    if(argc!=3){
        std::cout<<"usage: ./plain2grlbwt <input_plain_bwt> <output_grlbwt_fmt_bwt>"<<std::endl;
        exit(1);
    }

    auto plain_bwt = std::string(argv[1]);
    auto grlbwt_fmt_bwt = std::filesystem::path(std::string(argv[2])).replace_extension(".grlbwt");
    plain2grlbwt(plain_bwt, grlbwt_fmt_bwt);
}