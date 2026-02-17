#!/usr/bin/bash

export PAT_RT_EXPDIR_NAME=cache_miss_exp_tmp_dir
module load perftools

#the name of the file containing query patterns
pat_prefix=( "30bac" "covid" "hum" "kernel" )

bench_vlbt_rlbwt() {
	#running the experiments for the VLBT RLBWT
	rlbwt_prefix=( "30bac" "covid" "40hum" "kernel" )
	block_sizes=( 4096 16384 65536 262144 )
	rlbwt_dir=/home/ddiaz/gsa/ddiaz/vlbt_experiments/VLBT/build/rlbwt_dts
	truncate -s 0 vlbt_rlbwt_experiments.txt
	for (( j=0; j<4; j++ ));
	do
	        #pat_file=/home/ddiaz/gsa/ddiaz/vlbt_experiments/datasets/patterns/${pat_prefix[$j]}_patlen_105_npats_fil_50k.pat
	        pat_file=/home/ddiaz/gsa/ddiaz/vlbt_experiments/VLBT/build/${pat_prefix[$j]}_patlen_105_npats_fil_50k.pat
	        for (( i=0; i<4; i++ ));
	        do
	                file_id=vlbt_rlbwt_${rlbwt_prefix[$j]}_${block_sizes[$i]}
	                input_file=${rlbwt_dir}/${rlbwt_prefix[$j]}_rlbwt_${block_sizes[$i]}.rlbwt_vlt
	                echo ./vlbt_bench_cmiss+pat ${input_file} 5000000 ${pat_file}
	                ./vlbt_bench_cmiss+pat ${input_file} 5000000 ${pat_file}
			pat_report -O export -d counters -b fu ${PAT_RT_EXPDIR_NAME} \
			       	| grep "bench_iso" \
				| awk -v name="${file_id}" '
					BEGIN{FS=","}{
						print name, $1, $2, $3;
					}' >> vlbt_rlbwt_experiments.txt
	                rm -rf ${PAT_RT_EXPDIR_NAME}
	        done
	done
}

bench_mn_kp_rlbwt() {
	#running the experiments for the RLBWTs
	rlbwt_prefix=( "30bac" "covid" "40humans" "linux_kernel" )
	rlbwt_types=( "mn" "kp" )
	rlbwt_dir=/home/ddiaz/gsa/ddiaz/vlbt_experiments/rlbwts/build
	truncate -s 0 mn_kp_rlbwt_experiments.txt
	for (( j=0; j<4; j++ ));
	do
	        #pat_file=/home/ddiaz/gsa/ddiaz/vlbt_experiments/datasets/patterns/${pat_prefix[$j]}_patlen_105_npats_fil_50k.pat
	        pat_file=/home/ddiaz/gsa/ddiaz/vlbt_experiments/VLBT/build/${pat_prefix[$j]}_patlen_105_npats_fil_50k.pat
		for(( k=0; k<2; k++ ));
		do
			file_id=rlbwt_${rlbwt_types[$k]}_${rlbwt_prefix[$j]}
			input_file=rlbwt_${rlbwt_types[$k]}_${rlbwt_prefix[$j]}.rl${rlbwt_types[$k]}
	                echo ./rlbwt_bench_cmiss+pat ${rlbwt_dir}/${input_file} ${k} 5000000 ${pat_file}
	                ./rlbwt_bench_cmiss+pat ${rlbwt_dir}/${input_file} ${k} 5000000 ${pat_file}
			pat_report -O export -d counters -b fu ${PAT_RT_EXPDIR_NAME} \
				| grep "bench_iso" \
				| awk -v name="${file_id}" '
					BEGIN{FS=","}{
						print name, $1, $2, $3;
					}' >> mn_kp_rlbwt_experiments.txt
	                rm -rf ${PAT_RT_EXPDIR_NAME}
		done
	done
}

bench_vlbt_sri() {
	#running the experiments for the VLBT RLBWT
	sri_prefix=( "30bac" "covid" "40hum" "kernel" )
	block_sizes=( 4096 16384 65536 262144 )
	sri_subsamp=( 4 8 16 32 )
	sri_dir=/home/ddiaz/gsa/ddiaz/vlbt_experiments/VLBT/build/sri_dts
	truncate -s 0 vlbt_sri_experiments.txt 
	for (( j=0; j<4; j++ ));
	do
	        pat_file=/home/ddiaz/gsa/ddiaz/vlbt_experiments/VLBT/build/${pat_prefix[$j]}_patlen_105_npats_fil_50k.pat
	        for (( i=0; i<4; i++ ));
	        do
			for (( l=0; l<4; l++ ));
			do
	                	file_id=vlbt_sri_${sri_prefix[$j]}_${block_sizes[$i]}_${sri_subsamp[${l}]}
	                	input_file=${sri_dir}/${sri_prefix[$j]}_sri_${block_sizes[$i]}_${sri_subsamp[${l}]}.sri_vlt
	                	echo ./vlbt_bench_cmiss+pat ${input_file} 5000000 ${pat_file}
	                	./vlbt_bench_cmiss+pat ${input_file} 5000000 ${pat_file}
				pat_report -O export -d counters -b fu ${PAT_RT_EXPDIR_NAME} \
				       	| grep "bench_iso" \
					| awk -v name="${file_id}" '
						BEGIN{FS=","}{
							print name, $1, $2, $3;
						}' >> vlbt_sri_experiments.txt
	                	rm -rf ${PAT_RT_EXPDIR_NAME}
			done
	        done
	done
	sed 's/_/ /g' vlbt_sri_experiments.txt | awk '{print $1"_"$2,$4"_"$5,$3,$6"_"$7"_"$8,$9,$10,$11}' | sort -k4,4 -t' ' > vlbt_sri_experiments_parsed.txt
}

bench_ri() {
	#running the experiments for the RLBWTs
	ri_prefix=( "30bac" "covid" "40humans" "linux_kernel" )
	ri_dir=/home/ddiaz/gsa/ddiaz/vlbt_experiments/r-index-exp/build
	truncate -s 0 ri_experiments.txt 
	for (( j=0; j<4; j++ ));
	do
	        #pat_file=/home/ddiaz/gsa/ddiaz/vlbt_experiments/datasets/patterns/${pat_prefix[$j]}_patlen_105_npats_fil_50k.pat
      		pat_file=/home/ddiaz/gsa/ddiaz/vlbt_experiments/VLBT/build/${pat_prefix[$j]}_patlen_105_npats_fil_50k.pat
		file_id=r_index_${ri_prefix[$j]}
		input_file=${file_id}.ri
	        echo ./ri_bench_cmiss+pat ${ri_dir}/${input_file} ${pat_file}
	        ./ri_bench_cmiss+pat ${ri_dir}/${input_file} ${pat_file}
		pat_report -O export -d counters -b fu ${PAT_RT_EXPDIR_NAME} \
			| grep "bench_iso" \
			| awk -v name="${file_id}" '
				BEGIN{FS=","}{
					print name, $1, $2, $3;
				}' >> ri_experiments.txt
	        rm -rf ${PAT_RT_EXPDIR_NAME}
	done
}

bench_sri() {
	#running the experiments for the RLBWTs
	index_names=( "30bac" "covid" "40humans" "kernel" )
	s_value=( 8 12 16 20 )
	index_type=( 0 1 2 )
	ext=( "sri" "sri_vm" "sri_va" )
	truncate -s 0 sri_experiments.txt 
	sri_dir=/home/ddiaz/gsa/ddiaz/vlbt_experiments/sr-index/build
	for (( i=0; i<4; i++ ));
	do
      		pat_file=/home/ddiaz/gsa/ddiaz/vlbt_experiments/VLBT/build/${pat_prefix[$i]}_patlen_105_npats_fil_50k.pat
		for (( j=0; j<4; j++ ));
		do
			for (( k=0; k<3; k++ ));
			do
                		input_file=${index_names[$i]}_s_${s_value[$j]}_i_${index_type[$k]}.${ext[$k]}
	        		echo ./sri_bench_cmiss+pat ${sri_dir}/${input_file} ${pat_file}
	        		./sri_bench_cmiss+pat ${sri_dir}/${input_file} ${pat_file}
				pat_report -O export -d counters -b fu ${PAT_RT_EXPDIR_NAME} \
					| grep "bench_iso" \
					| awk -v name="${input_file}" '
						BEGIN{FS=","}{
							print name, $1, $2, $3;
						}' >> sri_experiments.txt
	        		rm -rf ${PAT_RT_EXPDIR_NAME}
			done
		done
	done
	sed 's/2.sri_va/va/' sri_experiments.txt | sed 's/1.sri_vm/vm/' | sed 's/0.sri /sri /' | sed 's/40humans/40hum/' | sed 's/_/ /g' | awk '{print "sr-index", $3"_"$5,$1,$8,$9,$10}' | sort -k4,4 -t' ' > sri_experiments_parsed.txt
}


bench_vlbt_rlbwt
bench_mn_kp_rlbwt
bench_vlbt_sri
bench_ri
bench_sri
