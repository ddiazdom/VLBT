#!/usr/bin/bash

export PAT_RT_EXPDIR_NAME=cache_miss_exp_tmp_dir
module load perftools

#the name of the file containing query patterns
pat_prefix=( "30bac" "covid" "hum" "kernel" )

bench_vlbt_rlbwt() {
	#running the experiments for the VLBT RLBWT
	rlbwt_prefix=( "30bac" "covid" "40hum" "kernel" )
	block_sizes=( 4096 16384 65536 262144 )
	rlbwt_dir=/path/to/folder/vlbt_experiments/VLBT/build/rlbwt_dts
	for (( j=0; j<4; j++ ));
	do
	        for (( i=0; i<4; i++ ));
	        do
	                file_id=vlbt_rlbwt_${rlbwt_prefix[$j]}_${block_sizes[$i]}
	                input_file=${rlbwt_dir}/${rlbwt_prefix[$j]}_rlbwt_${block_sizes[$i]}.rlbwt_vlt
			file_bytes=$(wc -c ${input_file} | cut -f1 -d' ')
			echo "vlbt_rlbwt ${block_sizes[$i]} ${rlbwt_prefix[$j]} ${file_bytes}"
	        done
	done
}

bench_mn_kp_rlbwt() {
	#running the experiments for the RLBWTs
	rlbwt_prefix=( "30bac" "covid" "40humans" "linux_kernel" )
	rlbwt_types=( "mn" "kp" )
	rlbwt_dir=/path/to/folder/vlbt_experiments/rlbwts/build

	for (( j=0; j<4; j++ ));
	do
		for(( k=0; k<2; k++ ));
		do
			file_id=rlbwt_${rlbwt_types[$k]}
			input_file=${rlbwt_dir}/rlbwt_${rlbwt_types[$k]}_${rlbwt_prefix[$j]}.rl${rlbwt_types[$k]}
			file_bytes=$(wc -c ${input_file} | cut -f1 -d' ')
			echo ${file_id} ${rlbwt_types[$k]} ${rlbwt_prefix[$j]} ${file_bytes}
		done
	done
}

bench_vlbt_sri() {
	#running the experiments for the VLBT RLBWT
	sri_prefix=( "30bac" "covid" "40hum" "kernel" )
	block_sizes=( 4096 16384 65536 262144 )
	sri_subsamp=( 4 8 16 32 )
	sri_dir=/path/to/folder/vlbt_experiments/VLBT/build/sri_dts
	for (( j=0; j<4; j++ ));
	do
	        for (( i=0; i<4; i++ ));
	        do
			for (( l=0; l<4; l++ ));
			do
	                	file_id=vlbt_sri_${sri_prefix[$j]}_${block_sizes[$i]}_${sri_subsamp[${l}]}
	                	input_file=${sri_dir}/${sri_prefix[$j]}_sri_${block_sizes[$i]}_${sri_subsamp[${l}]}.sri_vlt
				file_bytes=$(wc -c ${input_file} | cut -f1 -d' ')
				echo "vlbt_sri ${block_sizes[$i]}_${sri_subsamp[$l]} ${sri_prefix[$j]} ${file_bytes}"
			done
	        done
	done
}

bench_ri() {
	#running the experiments for the RLBWTs
	ri_prefix=( "30bac" "covid" "40humans" "linux_kernel" )
	ri_dir=/path/to/folder/vlbt_experiments/r-index-exp/build
	for (( j=0; j<4; j++ ));
	do
		file_id=r_index_${ri_prefix[$j]}
		input_file=${ri_dir}/${file_id}.ri
		file_bytes=$(wc -c ${input_file} | cut -f1 -d' ')
		echo "r_index r_index ${ri_prefix[$j]} ${file_bytes}"
	done
}

bench_sri() {
	#running the experiments for the RLBWTs
	index_names=( "30bac" "covid" "40humans" "kernel" )
	s_value=( 8 12 16 20 )
	index_type=( 0 1 2 )
	ext=( "sri" "sri_vm" "sri_va" )
	tp=( "sri" "vm" "va" )
	sri_dir=/path/to/folder/vlbt_experiments/sr-index/build
	for (( i=0; i<4; i++ ));
	do
		for (( j=0; j<4; j++ ));
		do
			for (( k=0; k<3; k++ ));
			do
                		input_file=${sri_dir}/${index_names[$i]}_s_${s_value[$j]}_i_${index_type[$k]}.${ext[$k]}
				file_bytes=$(wc -c ${input_file} | cut -f1 -d' ')
				echo "sr-index ${s_value[$j]}_${tp[$k]} ${index_names[$i]} ${file_bytes}"
			done
		done
	done
}

bench_vlbt_rlbwt
bench_mn_kp_rlbwt
bench_vlbt_sri
bench_ri
bench_sri
