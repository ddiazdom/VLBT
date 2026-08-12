/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */

#ifndef VLBT_DEF_SCAN_H
#define VLBT_DEF_SCAN_H

#include <cstdint>

template<bool overflow16, bool overflow32=false>
static inline uint8_t access_scl_8(const uint16_t* stream, uint8_t sigma, uint64_t idx){
    return 0;
}
template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline uint8_t access_scl_16(const uint16_t* stream, uint8_t sigma, uint64_t idx){
    return 0;
}
template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline uint8_t access_scl_32(const uint16_t* stream, uint8_t sigma, uint64_t idx){
    return 0;
}
template<uint8_t bytes_per_run>
static inline uint8_t access_scl_64(const uint16_t* stream, uint8_t sigma, uint64_t idx){
    return 0;
    /*uint8_t alpha_bits = 4;
    uint8_t alpha_mask = 15;

    uint64_t r_len[2]={0};
    r_len[1] = stream[0]>>alpha_bits;
    size_t acc = r_len[1], i=0, rank=0;
    while(acc<idx){
        rank+= r_len[(stream[i] & alpha_mask)==sym];
        r_len[1] = stream[++i]>>alpha_bits;
        acc+= r_len[1];
    }
    return rank + (idx-(acc-r_len[1]))*((stream[i]&15)==sym);*/
}

template<bool overflow16, bool overflow32=false>
static inline uint64_t inv_select_scl_8(const uint16_t* stream, uint8_t sigma, uint64_t idx){
    return 0;
}
template<bool vbyte_compressed, bool overflow8, bool overflow16=false>
static inline uint64_t inv_select_scl_16(const uint16_t* stream, uint8_t sigma, uint64_t idx){
    return 0;
}
template<bool vbyte_compressed, uint8_t bytes_per_run>
static inline uint64_t inv_select_scl_32(const uint16_t* stream, uint8_t sigma, uint64_t idx){
    return 0;
}
template<uint8_t bytes_per_run>
static inline uint64_t inv_select_scl_64(const uint16_t* stream, uint8_t sigma, uint64_t idx){
    return 0;
}

template<bool overflow16, bool overflow32, bool check_head>
static inline uint64_t rank_scl_8(const uint16_t* stream, uint8_t sigma, uint64_t idx, uint8_t sym){
    return 0;
}
template<bool vbyte_compressed, bool overflow8, bool overflow16, bool check_head>
static inline uint64_t rank_scl_16(const uint16_t* stream, uint8_t sigma, uint64_t idx, uint8_t sym){
    return 0;
}
template<bool vbyte_compressed, uint8_t bytes_per_run, bool check_head>
static inline uint64_t rank_scl_32(const uint16_t* stream, uint8_t sigma, uint64_t idx, uint8_t sym){
    return 0;
}
template<uint8_t bytes_per_run, bool check_head>
static inline uint64_t rank_scl_64(const uint16_t* stream, uint8_t sigma, uint64_t idx, uint8_t sym){
    return 0;
}

//AArch64 compilers define __ARM_NEON (ACLE); __ARM_NEON__ is the legacy AArch32 spelling
#if defined(__ARM_NEON__) || defined(__ARM_NEON)
#include "scan_neon.h"

#define INV_SELECT_8 inv_select_neon_8x16
#define INV_SELECT_16 inv_select_neon_16x8
#define INV_SELECT_32 inv_select_neon_32x4
#define INV_SELECT_64 inv_select_neon_64x2

#define ACCESS_8 access_neon_8x16
#define ACCESS_16 access_neon_16x8
#define ACCESS_32 access_neon_32x4
#define ACCESS_64 access_neon_64x2

#define RANK_8 rank_neon_8x16
#define RANK_16 rank_neon_16x8
#define RANK_32 rank_neon_32x4
#define RANK_64 rank_neon_64x2

#define RANGE_RANK_8 range_rank_neon_8x16
#define RANGE_RANK_16 range_rank_neon_16x8
#define RANGE_RANK_32 range_rank_neon_32x4
#define RANGE_RANK_64 range_rank_neon_64x2

#define FIRST_RUN_8 first_run_neon_8x16
#define FIRST_RUN_16 first_run_neon_16x8
#define FIRST_RUN_32 first_run_neon_32x4
#define FIRST_RUN_64 first_run_neon_64x2

#define SUCC_8 succ_neon_8x16
#define SUCC_16 succ_neon_16x8
#define SUCC_32 succ_neon_32x4
#define SUCC_64 succ_neon_64x2

#define GET_PHI_RUN_8 get_phi_run_neon_8x16
#define GET_PHI_RUN_16 get_phi_run_neon_16x8
#define GET_PHI_RUN_32 get_phi_run_neon_32x4
#define GET_PHI_RUN_64 get_phi_run_neon_64x2

#elif defined(__AVX2__)
#include "scan_sse42.h"
#include "scan_avx2.h"

#define INV_SELECT_8 inv_select_sse42_8x16
#define INV_SELECT_16 inv_select_sse42_16x8
#define INV_SELECT_32 inv_select_sse42_32x4
#define INV_SELECT_64 inv_select_sse42_64x2

#define ACCESS_8 access_sse42_8x16
#define ACCESS_16 access_sse42_16x8
#define ACCESS_32 access_sse42_32x4
#define ACCESS_64 access_sse42_64x2

#define RANK_8 rank_sse42_8x16
#define RANK_16 rank_sse42_16x8
#define RANK_32 rank_sse42_32x4
#define RANK_64 rank_sse42_64x2

#define RANGE_RANK_8 range_rank_sse42_8x16
#define RANGE_RANK_16 range_rank_sse42_16x8
#define RANGE_RANK_32 range_rank_sse42_32x4
#define RANGE_RANK_64 range_rank_sse42_64x2

#define FIRST_RUN_8 first_run_sse42_8x16
#define FIRST_RUN_16 first_run_sse42_16x8
#define FIRST_RUN_32 first_run_sse42_32x4
#define FIRST_RUN_64 first_run_sse42_64x2

#define SUCC_8 succ_sse42_8x16
#define SUCC_16 succ_sse42_16x8
#define SUCC_32 succ_sse42_32x4
#define SUCC_64 succ_sse42_64x2

#define GET_PHI_RUN_8 get_phi_run_avx2_8x32
#define GET_PHI_RUN_16 get_phi_run_avx2_16x16
#define GET_PHI_RUN_32 get_phi_run_sse42_32x4
#define GET_PHI_RUN_64 get_phi_run_sse42_64x2

#elif defined(__SSE4_2__)
#include "scan_sse42.h"

#define INV_SELECT_8 inv_select_sse42_8x16
#define INV_SELECT_16 inv_select_sse42_16x8
#define INV_SELECT_32 inv_select_sse42_32x4
#define INV_SELECT_64 inv_select_sse42_64x2

#define ACCESS_8 access_sse42_8x16
#define ACCESS_16 access_sse42_16x8
#define ACCESS_32 access_sse42_32x4
#define ACCESS_64 access_sse42_64x2

#define RANK_8 rank_sse42_8x16
#define RANK_16 rank_sse42_16x8
#define RANK_32 rank_sse42_32x4
#define RANK_64 rank_sse42_64x2

#define RANGE_RANK_8 range_rank_sse42_8x16
#define RANGE_RANK_16 range_rank_sse42_16x8
#define RANGE_RANK_32 range_rank_sse42_32x4
#define RANGE_RANK_64 range_rank_sse42_64x2

#define FIRST_RUN_8 first_run_sse42_8x16
#define FIRST_RUN_16 first_run_sse42_16x8
#define FIRST_RUN_32 first_run_sse42_32x4
#define FIRST_RUN_64 first_run_sse42_64x2

#define SUCC_8 succ_sse42_8x16
#define SUCC_16 succ_sse42_16x8
#define SUCC_32 succ_sse42_32x4
#define SUCC_64 succ_sse42_64x2

#define GET_PHI_RUN_8 get_phi_run_sse42_8x16
#define GET_PHI_RUN_16 get_phi_run_sse42_16x8
#define GET_PHI_RUN_32 get_phi_run_sse42_32x4
#define GET_PHI_RUN_64 get_phi_run_sse42_64x2

#else
//TODO scalar fallback, pending implementation. The *_scl_* functions above are stubs that
//return 0, and RANGE_RANK_*, FIRST_RUN_*, SUCC_* and GET_PHI_RUN_* have no scalar
//counterpart yet, so enabling this path would answer every query with 0
//#define INV_SELECT_8 inv_select_scl_8
//#define INV_SELECT_16 inv_select_scl_16
//#define INV_SELECT_32 inv_select_scl_32
//#define INV_SELECT_64 inv_select_scl_64

//#define ACCESS_8 access_scl_8
//#define ACCESS_16 access_scl_16
//#define ACCESS_32 access_scl_32
//#define ACCESS_64 access_scl_64

//#define RANK_8 rank_scl_8
//#define RANK_16 rank_scl_16
//#define RANK_32 rank_scl_32
//#define RANK_64 rank_scl_64
#error "VLBT needs NEON, SSE4.2 or AVX2. Compile with -march=native (or an equivalent flag)."
#endif

#endif //VLBT_DEF_SCAN_H
