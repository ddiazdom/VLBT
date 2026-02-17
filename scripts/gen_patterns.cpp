/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */

#include <cassert>
#include <string>
#include <iostream>
#include <random>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <algorithm>

std::vector<uint64_t> sample_random_positions(std::unordered_set<uint64_t>& seen_positions, uint64_t n, uint64_t x) {

    if (x > n) throw std::invalid_argument("x cannot be larger than n");

    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist(0, n - 1);
    std::unordered_set<uint64_t> new_positions;

    while(new_positions.size() < x) {
        uint64_t val = dist(rng);
        if(seen_positions.find(val)==seen_positions.end()){
            new_positions.insert(val); // insert is a no-op if val already exists
        }
    }
    return {new_positions.begin(), new_positions.end()};
}

int main(int argc, char** argv){

    if(argc!=6){
        std::cout<<"usage: ./gen_patterns <input_file> <length> <number> <patterns file> <forbidden symbol>"<<std::endl;
        exit(1);
    }

    auto input_file = std::string(argv[1]);

    char *pend;
    long int len = strtol(argv[2], &pend, 10);
    long int number = strtol(argv[3], &pend, 10);
    auto output_file = std::string(argv[4]);

    long int fb_sym_int = strtol(argv[5], &pend, 10);
    assert(fb_sym_int<256);
    auto fb_sym = static_cast<uint8_t>(fb_sym_int);
    auto f_size = static_cast<long int>(std::filesystem::file_size(input_file));

    number = std::min(f_size-len+1, number);
    len = std::min(f_size, len);
    size_t rem = number;

    std::unordered_set<uint64_t> seen_positions;

    std::basic_string<uint8_t> buffer(len, 0);
    std::ifstream ifs(input_file, std::ios::binary);
    std::ofstream ofs(output_file, std::ios::out);

    std::string file_name = std::filesystem::path(input_file).filename();
    std::string header = "# number="+std::to_string(number)+" length="+std::to_string(len)+" file="+file_name+" forbidden="+std::to_string((int)fb_sym)+"\n";
    ofs.write(header.data(), static_cast<std::streamsize>(header.size()));

    while(rem>0){

        std::vector<uint64_t> rd_pos = sample_random_positions(seen_positions, f_size-len+1, rem);
        std::sort(rd_pos.begin(), rd_pos.end());

        for(auto const& pos : rd_pos){
            ifs.seekg(static_cast<long long>(pos));
            ifs.read(reinterpret_cast<char *>(buffer.data()), len);
            if(buffer.find(fb_sym)==std::string::npos){
                ofs.write(reinterpret_cast<char *>(buffer.data()), len);
                rem--;
            }
            seen_positions.insert(pos);
        }
    }

    ifs.close();
    ofs.close();
    std::cout<<"We extracted "<<number<<" random patterns of length "<<len<<" from file "<<file_name<<std::endl;
}
