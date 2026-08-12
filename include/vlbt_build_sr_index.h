/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */

#ifndef VLBT_BUILD_SR_INDEX_H
#define VLBT_BUILD_SR_INDEX_H

#include <cassert>
#include <string>
#include <filesystem>

#include "vlbt_build_phi.h"
#include "vlbt_build_bwt.h"
#include "subsamping.h"
#include "logger.h"

template<class sa_samp_type, class sr_index_type>
void build_sr_index(sr_index_type& index, const std::string& bwt_file,
                    const std::string& sa_heads_file, const std::string& sa_tails_file,
                    BWT_FORMAT bwt_file_fmt, size_t sri_samp_val,
                    std::string tmp_dir="./"){

    using bwt_th_type = typename sr_index_type::bwt_t;
    using phi_type = typename sr_index_type::phi_t;

    if (!std::filesystem::exists(bwt_file)) {
        LOG_ERROR("File "+bwt_file+" does not exist");
        exit(1);
    }
    if (!std::filesystem::exists(sa_heads_file)) {
        LOG_ERROR("File "+sa_heads_file+" does not exist");
        exit(1);
    }
    if (!std::filesystem::exists(sa_tails_file)) {
        LOG_ERROR("File "+sa_tails_file+" does not exist");
        exit(1);
    }
    if (std::filesystem::file_size(sa_heads_file)!=std::filesystem::file_size(sa_tails_file)) {
        LOG_ERROR("The number of symbols in the SA files does not match");
        exit(1);
    }

    tmp_workspace twd(tmp_dir, true, "vlbt_sri");

    const std::string sa_heads_subsamp_file = twd.get_file("ssa_subsamp");
    const std::string sa_tails_subsamp_file = twd.get_file("esa_subsamp");

    uint64_t n = std::filesystem::file_size(bwt_file);//number of symbols in the BWT

    subsample_sa_samples<sa_samp_type>(sa_heads_file, sa_tails_file, sri_samp_val, sa_heads_subsamp_file, sa_tails_subsamp_file, n);

    LOG_INFO("Building the data structure");
    build_bwt_th_int<bwt_th_type, sa_samp_type>(index.bwt, bwt_file, bwt_file_fmt,
                                                   sri_samp_val, sa_heads_subsamp_file, twd);

    index.phi.subsamp_step = sri_samp_val;
    build_phi<phi_type, sa_samp_type>(index.phi, sa_tails_subsamp_file, twd);
}
#endif //VLBT_BUILD_SR_INDEX_H