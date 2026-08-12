/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */

#ifndef VLBT_VLBT_COMMON_H
#define VLBT_VLBT_COMMON_H

#include <fstream>
#include "logger.h"
#include "utils.h"

enum VLBT_TYPE{
    RLBWT=0,
    RLBWT_WITH_TOEHOLDS=1,
    SRI_VALID_AREA=2,
};

enum phi_variant {
    NO_VALID_AREA = 0,
    WITH_VALID_AREA = 1,
    NO_SUBSAMPLING=2
};

enum node_class {
    INTERNAL,
    LEAF
};

//template parameters that were serialized,
//we used them to ensure the correct load of class templates from the disk
struct temp_param_t {
    VLBT_TYPE tag;
    size_t b_size;
    size_t b_size2;
};

inline temp_param_t read_template_param(const std::string& str) {
    temp_param_t tp{};
    std::ifstream ifs(str, std::ios::binary);
    if (!ifs) {
        LOG_ERROR("Cannot open index file "+str);
        exit(1);
    }
    load_elm(ifs, tp.tag);
    load_elm(ifs, tp.b_size);
    tp.b_size2 = tp.b_size;
    if (tp.tag==SRI_VALID_AREA) {
        //b_size is the block size for the BWT
        //b_size2 is the block size for phi
        load_elm(ifs, tp.b_size2);
    }
    if (!ifs) {
        LOG_ERROR("File "+str+" is too short to be a valid index");
        exit(1);
    }
    ifs.close();
    return tp;
}
#endif //VLBT_VLBT_COMMON_H