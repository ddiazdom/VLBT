//
// Created by Diaz, Diego on 17.11.2025.
//

#ifndef VLBT_VLBT_COMMON_H
#define VLBT_VLBT_COMMON_H
enum VLBT_TYPE{
    RLBWT=0,
    RLBWT_WITH_TOEHOLDS=1,
    SRI_VALID_AREA=2,
};

//template parameters that were serialized
//we used them to ensure the correct load of class templates from disk
struct temp_param_t {
    VLBT_TYPE tag;
    size_t b_size;
    size_t b_size2;
};

inline temp_param_t read_template_param(const std::string& str) {
    temp_param_t tp;
    std::ifstream ifs(str);
    load_elm(ifs, tp.tag);
    load_elm(ifs, tp.b_size);
    tp.b_size2 = tp.b_size;
    if (tp.tag==SRI_VALID_AREA) {
        //b_size is the block size for the BWT
        //b_size2 is the block size for phi
        load_elm(ifs, tp.b_size2);
    }
    ifs.close();
    return tp;
}
#endif //VLBT_VLBT_COMMON_H