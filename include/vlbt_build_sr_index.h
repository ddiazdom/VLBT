//
// Created by Diaz, Diego on 29.5.2025.
//

#ifndef VLBT_BUILD_SR_INDEX_H
#define VLBT_BUILD_SR_INDEX_H

#include <cassert>
#include <string>
#include <filesystem>

#include "vlbt_build_phi.h"
#include "vlbt_build_bwt.h"
#include "subsamping.h"

template<class sa_samp_type, class sr_index_type>
void build_sr_index(sr_index_type& index, const std::string& input_prefix,
                    BWT_FORMAT bwt_file_fmt, size_t sri_samp_val,
                    std::string tmp_dir="./"){

    using bwt_th_type = typename sr_index_type::bwt_t;
    using phi_type = typename sr_index_type::phi_t;

    const std::string bwt_file = input_prefix+".bwt";
    const std::string sa_heads_file = input_prefix+".ssa";
    const std::string sa_tails_file = input_prefix+".esa";
    assert(std::filesystem::exists(bwt_file));
    assert(std::filesystem::exists(sa_heads_file));
    assert(std::filesystem::exists(sa_tails_file));
    assert(std::filesystem::file_size(sa_heads_file)==std::filesystem::file_size(sa_tails_file));

    tmp_workspace twd(tmp_dir, true, "vlbt_sri");

    const std::string sa_heads_subsamp_file = twd.get_file("ssa_subsamp");
    const std::string sa_tails_subsamp_file = twd.get_file("esa_subsamp");

    std::cout<<"Subsapling "<<std::endl;
    uint64_t n = std::filesystem::file_size(bwt_file);//number of symbols in the BWT

    //TODO testing
    //uint64_t n = 267410983471;
    //std::cout<<"Fixed the number of symbols, retore it"<<std::endl;
    //

    subsample_sa_samples<sa_samp_type>(sa_heads_file, sa_tails_file, sri_samp_val, sa_heads_subsamp_file, sa_tails_subsamp_file, n);

    std::cout<<"Building the data structure "<<std::endl;
    build_bwt_th_int<bwt_th_type, sa_samp_type>(index.bwt, bwt_file, bwt_file_fmt,
                                                   sri_samp_val, sa_heads_subsamp_file, twd);

    /*std::string samp_sa_file = input_prefix+".sa_samples";
    std::string str_ranges_file = input_prefix+".str_ranges";
    std::string bwt_file = input_prefix+".ebwt";
    std::string ssamp_heads_file = output_prefix+".ssamp_heads";
    std::string ssamp_tails_file = output_prefix+".ssamp_tails";
    subsample_sa_samples<sa_samp_type>(samp_sa_file, str_ranges_file, sri_samp_val, ssamp_heads_file, ssamp_tails_file);*/
    //build_bwt_th<bwt_type, sa_samp_type>(index.bwt, bwt_file, sri_samp_val, GRL_BWT, ssamp_heads_file, tmp_dir);

    index.phi.subsamp_step = sri_samp_val;
    build_phi<phi_type, sa_samp_type>(index.phi, sa_tails_subsamp_file, twd);
}
#endif //VLBT_BUILD_SR_INDEX_H