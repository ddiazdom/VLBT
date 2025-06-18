//
// Created by Diaz, Diego on 22.5.2025.
//

#ifndef VLBT_DEF_SCAN_H
#define VLBT_DEF_SCAN_H

#if defined(__ARM_NEON__)
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

/*#elif defined(__AVX2__)
#include "scan_avx2.h"

#define INV_SELECT_8 inv_select_avx2_8x32
#define INV_SELECT_16 inv_select_avx2_16x16
#define INV_SELECT_32 inv_select_avx2_32x8
#define INV_SELECT_64 inv_select_avx2_64x4

#define ACCESS_8 access_avx2_8x32
#define ACCESS_16 access_avx2_16x16
#define ACCESS_32 access_avx2_32x8
#define ACCESS_64 access_avx2_64x4

#define RANK_8 rank_avx2_8x16
#define RANK_16 rank_avx2_16x8
#define RANK_32 rank_avx2_32x4
#define RANK_64 rank_avx2_64x2
*/

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

#else
#define INV_SELECT_8 inv_select_scl_8
#define INV_SELECT_16 inv_select_scl_16
#define INV_SELECT_32 inv_select_scl_32
#define INV_SELECT_64 inv_select_scl_64

#define ACCESS_8 access_scl_8
#define ACCESS_16 access_scl_16
#define ACCESS_32 access_scl_32
#define ACCESS_64 access_scl_64

#define RANK_8 rank_scl_8
#define RANK_16 rank_scl_16
#define RANK_32 rank_scl_32
#define RANK_64 rank_scl_64
#endif
#endif //VLBT_DEF_SCAN_H
