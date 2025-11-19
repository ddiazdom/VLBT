//
// Created by Diaz, Diego on 19.11.2025.
//

#include <iostream>
#include <fstream>
#include <string>
#include "utils.h"

int main(int argc, char** argv){

    if(argc!=3){
        std::cout<<"usage: ./pizza&chilli2fasta <pat_file> <output_fasta>"<<std::endl;
        exit(1);
    }

    auto pat_file = std::string(argv[1]);
    auto output_file = std::string(argv[2]);

    std::ifstream ifs(pat_file);
    std::ofstream ofs(output_file);

    std::string header;
    std::getline(ifs, header);

    size_t n_pats = get_n_patterns(header);
    size_t pat_len = get_patterns_len(header);

    std::string pat;
    size_t count=0;
    pat.resize(pat_len);

    std::cout<<"Parsing "<<n_pats<<" patterns of length "<<pat_len<<std::endl;

    for(ulint i=0;i<n_pats;++i){
        for(ulint j=0;j<pat_len;++j){
            char c;
            ifs.get(c);
            pat[j] = c;
        }
        count++;
        ofs << ">pat_"<<count<<'\n'<<pat<<"\n";
    }

    ofs.close();
    ifs.close();
}