#!/usr/bin/bash

export PAT_RT_PERFCTR=PAPI_L1_DCM,PAPI_L2_DCM
export PAT_RT_EXPDIR_NAME=cache_miss_exp_tmp_dir

module load perftools

#compile
echo "Compiling VLBT"
cc -o vlbt_bench_cmiss /path/to/folder/vlbt_experiments/VLBT/test_suite/vlbt_bench_cmisses.cpp -g -O3 -msse4.2 -fno-omit-frame-pointer
pat_build -f -w -T /bench_iso vlbt_bench_cmiss

echo "Compiling RLBWT"
cc -o rlbwt_bench_cmiss /path/to/folder/vlbt_experiments/VLBT/test_suite/rlbwt_bench_cmisses.cpp -g -O3 -msse4.2 -fno-omit-frame-pointer -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64
pat_build -f -w -T /bench_iso rlbwt_bench_cmiss

echo "Compiling the r-index"
cc -o ri_bench_cmiss /path/to/folder/vlbt_experiments/VLBT/test_suite/ri_bench_cmisses.cpp -g -O3 -msse4.2 -fno-omit-frame-pointer -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64
pat_build -f -w -T /bench_iso ri_bench_cmiss

echo "Compiling the sri-index"
cc -o sri_bench_cmiss /path/to/folder/vlbt_experiments/VLBT/test_suite/sri_bench_cmisses.cpp -g -O3 -msse4.2 -fno-omit-frame-pointer -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 -I /path/to/folder/vlbt_experiments/VLBT/test_suite/sr-index/include/json/include
pat_build -f -w -T /bench_iso sri_bench_cmiss
