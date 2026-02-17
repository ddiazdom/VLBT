## This script summarizes the experiments for VLBT

library(plyr)
library(dplyr)
library(stringr)

compare_ratio <- function(C, x) {
  ifelse(C >= x, C / x, -(x / C))
}

count_time_rlbwt <- read.csv("~/Desktop/vlbt-exp-results/count_exp_data_rlbwt.csv")
count_time_csa <- read.csv("~/Desktop/vlbt-exp-results/count_exp_data_sri.csv")

locate_time_csa <- read.csv("~/Desktop/vlbt-exp-results/locate_exp_data.csv")

count_cmiss_rlbwt_and_csa <- read.csv("~/Desktop/vlbt-exp-results/vlbt-count-cmiss-exp.csv")
locate_cmiss_csa <- read.csv("~/Desktop/vlbt-exp-results/vlbt-locate-cmiss-exp.csv")

vlbt_rlbwt_stacked_space <- read.delim("~/Desktop/vlbt-exp-results/vlbt_rlbwt_stacked_space_data.tsv")
vlbt_sri_stacked_space <- read.delim("~/Desktop/vlbt-exp-results/vlbt_sri_stacked_space_data.tsv")

#count data for rlBWTs
data <- count_time_rlbwt[count_time_rlbwt$Pattern_length==105,]
data$nanosecs_pat <- as.numeric(data$nanosecs_pat)
tapply(data,
       as.factor(data$Input_file),
       function(x){
         stats <- tapply(x,
                as.factor(x$Index_class),
                function(y){
                  c(round(min(y$nanosecs_pat/1000), digits=4),
                    round(max(y$nanosecs_pat/1000), digits=4),
                    round(min(y$bits_per_symbol), digits=4),
                    round(max(y$bits_per_symbol), digits=4))
                }
                )
         dt <- ldply(stats, rbind)
         colnames(dt) <- c("tool", "min_speed", "max_speed", "min_space", "max_space")
         dt$ratio_speed <- paste(round(dt$min_speed/dt$min_speed[nrow(dt)], digits=2),
                                 round(dt$min_speed/dt$max_speed[nrow(dt)], digits=2),
                                 sep = " - ")
         
         dt$ratio_space <- paste(round(dt$min_space/dt$min_space[nrow(dt)], digits=2),
                                 round(dt$min_space/dt$max_space[nrow(dt)], digits=2),
                                 sep = " - ")
         dt[,c(1,6,7)]
       }
)

#performance as we vary block size in RLBWTs
data2 <- count_time_rlbwt[count_time_rlbwt$Pattern_length==105 & count_time_rlbwt$Index_class=="vlbt-bwt",]
tapply(data2,
       as.factor(data2$Input_file),
       function(x){
         stats <- tapply(x,
                as.factor(x$Index_class),
                function(y){
                  y$Index_type <- as.numeric(str_split_i(y$Index_type, "_", 3))
                  y[order(y$Index_type), ]
                }
         )
         dt <- ldply(stats, rbind)[,c(1,4,5,9)]
         dt$bps_ratio <- dt$bits_per_symbol/dt$bits_per_symbol[1]
         dt$speed_ratio <- dt$nanosecs_pat/dt$nanosecs_pat[1]
         dt
       }
)

data2 <- count_time_csa[count_time_csa$Pattern_length==105 & count_time_csa$Index_class %in% c("vlbt-sri-va-6", "vlbt-sri-va-7", "vlbt-sri-va-8", "vlbt-sri-va-9"),]
tapply(data2,
       as.factor(data2$Input_file),
       function(x){
         stats <- tapply(x,
                         as.factor(x$Index_class),
                         function(y){
                           y$Index_typeA <- as.numeric(str_split_i(y$Index_type, "_", 5))
                           y$Index_typeB <- as.numeric(str_split_i(y$Index_type, "_", 7))
                           y <- y[order(y$Index_typeA, y$Index_typeB), ]
                           y$bps_ratio <- y$bits_per_symbol/y$bits_per_symbol[1]
                           y$speed_ratio <- y$nanosecs_pat/y$nanosecs_pat[1]
                           y
                         }
         )
         #dt <- ldply(stats, rbind)[,c(1,4,5,9)]
         dt <- ldply(stats, rbind)[, c(1,11,12,13,14,15)]
         dt
       }
)

#count data in compressed suffix arrays
data3 <- count_time_csa[count_time_csa$Pattern_length==105,]
data3$nanosecs_pat <- as.numeric(data3$nanosecs_pat)
tapply(data3,
       as.factor(data3$Input_file),
       function(x){
         stats <- tapply(x,
                         as.factor(x$Index_class),
                         function(y){
                           #this is to compare sr-indexes
                           #y <- y[!(y$Index_type %in% c("sri_s_12", "sri_s_20",
                           #                            "sri_valid_marks_s_12", "sri_valid_marks_s_20",
                           #                            "sri_valid_area_s_12", "sri_valid_area_s_20",
                           #                            "SRI_VALID_AREA_b_4096_s_4", "SRI_VALID_AREA_b_4096_s_32",
                           #                            "SRI_VALID_AREA_b_16384_s_4", "SRI_VALID_AREA_b_16384_s_32",
                           #                            "SRI_VALID_AREA_b_65536_s_4", "SRI_VALID_AREA_b_65536_s_32",
                           #                            "SRI_VALID_AREA_b_262144_s_4", "SRI_VALID_AREA_b_262144_s_32"
                           #                            )),]
                           c(round(min(y$nanosecs_pat/1000), digits=4),
                             round(max(y$nanosecs_pat/1000), digits=4),
                             round(mean(y$nanosecs_pat/1000), digits=4),
                             round(min(y$bits_per_symbol), digits=4),
                             round(max(y$bits_per_symbol), digits=4),
                             round(mean(y$bits_per_symbol), digits=4))
                         }
         )
         dt <- ldply(stats, rbind)
         colnames(dt) <- c("tool", "min_speed", "max_speed", "mean_speed", "min_space", "max_space", "mean_space")
         dt <- dt[dt$tool!="sri" & dt$tool!="sri-vm",]
         mean_speed_vlbt <- sort(dt[str_detect(dt$tool, "^vlbt-sri-va"),4])
         mean_space_vlbt <- sort(dt[str_detect(dt$tool, "^vlbt-sri-va"),7])
         
         dt$ratio_mean_speed <- paste(round(dt$mean_speed/mean_speed_vlbt[1], digits=2),
                                      round(dt$mean_speed/mean_speed_vlbt[length(mean_speed_vlbt)], digits=2),
                                      sep = " - ")
         dt$ratio_mean_space <- paste(round(dt$mean_space/mean_space_vlbt[1], digits=2),
                                      round(dt$mean_space/mean_space_vlbt[length(mean_space_vlbt)], digits=2),
                                      sep = " - ")
         #dt[,c(1, 8,9)]
         dt
       }
)

#locate time in compressed suffix arrays
tapply(locate_time_csa,
       as.factor(locate_time_csa$Input_file),
       function(x){
         stats <- tapply(x,
                         as.factor(x$Index_class),
                         function(y){
                           #this is only to compare sri-indexes
                           #y <- y[!(y$Index_type %in% c("sri_s_12", "sri_s_20",
                           #                             "sri_valid_marks_s_12", "sri_valid_marks_s_20",
                           #                             "sri_valid_area_s_12", "sri_valid_area_s_20", 
                           #                             "SRI_VALID_AREA_b_4096_s_4", "SRI_VALID_AREA_b_4096_s_32",
                           #                             "SRI_VALID_AREA_b_16384_s_4", "SRI_VALID_AREA_b_16384_s_32",
                           #                             "SRI_VALID_AREA_b_65536_s_4", "SRI_VALID_AREA_b_65536_s_32",
                           #                             "SRI_VALID_AREA_b_262144_s_4", "SRI_VALID_AREA_b_262144_s_32"
                           #)),]
                           c(round(min(y$nanosecs_occ/1000), digits=4),
                             round(max(y$nanosecs_occ/1000), digits=4),
                             round(mean(y$nanosecs_occ/1000), digits=4),
                             round(min(y$bits_per_symbol), digits=4),
                             round(max(y$bits_per_symbol), digits=4),
                             round(mean(y$bits_per_symbol), digits=4))
                         }
         )
         dt <- ldply(stats, rbind)
         colnames(dt) <- c("tool", "min_speed", "max_speed", "mean_speed", "min_space", "max_space", "mean_space")
         dt <- dt[dt$tool!="sri" & dt$tool!="sri-vm",]
         mean_speed_vlbt <- sort(dt[str_detect(dt$tool, "^vlbt-sri-va"),4])
         mean_space_vlbt <- sort(dt[str_detect(dt$tool, "^vlbt-sri-va"),7])
         
         dt$ratio_mean_speed <- paste(round(dt$mean_speed/mean_speed_vlbt[1], digits=2),
                                      round(dt$mean_speed/mean_speed_vlbt[length(mean_speed_vlbt)], digits=2),
                                      sep = " - ")
         dt$ratio_mean_space <- paste(round(dt$mean_space/mean_space_vlbt[1], digits=2),
                                      round(dt$mean_space/mean_space_vlbt[length(mean_space_vlbt)], digits=2),
                                      sep = " - ")
         dt
       }
)

#====== cache misses stats

# count in run-length BWTs and compressed suffix arrays
count_miss_stats <- tapply(count_cmiss_rlbwt_and_csa,
       as.factor(count_cmiss_rlbwt_and_csa$input_file),
       function(x){
         stats <- tapply(x,
                         as.factor(x$data_structure),
                         function(y){
                           
                           #this is only to compare sri-indexes
                           y <- y[!(y$parameter %in% c("12_sri","20_sri",
                                                        "12_vm", "20_vm",
                                                        "12_va", "20_va",
                                                        "4096_4", "4096_32",
                                                        "16384_4", "16384_32",
                                                        "65536_4", "65536_32",
                                                        "262144_4","262144_32"
                           )),]
                           #print(y)
                           c(min(y$LD1_data_misses_symbol),
                             max(y$LD1_data_misses_symbol),
                             mean(y$LD1_data_misses_symbol),
                             round(min(y$bits_per_symbol), digits=4),
                             round(max(y$bits_per_symbol), digits=4),
                             round(mean(y$bits_per_symbol), digits=4))
                         }
         )
         dt <- ldply(stats, rbind)
         colnames(dt) <- c("tool", "min_misses", "max_misses", "mean_misses", "min_space", "max_space", "mean_space")
         dt <- dt[dt$tool!="sr-index-sri" & dt$tool!="sr-index-vm",]
         
         min_misses_vlbt_bwt <- dt[str_detect(dt$tool, "^vlbt_rlbwt"),2]
         max_misses_vlbt_bwt <- dt[str_detect(dt$tool, "^vlbt_rlbwt"),3]
         
         min_space_vlbt_bwt <- dt[str_detect(dt$tool, "^vlbt_rlbwt"),5]
         max_space_vlbt_bwt <- dt[str_detect(dt$tool, "^vlbt_rlbwt"),6]
         
         rlbwt_data <- dt[dt$tool %in% c("rlbwt_kp", "rlbwt_mn", "vlbt_rlbwt"),]
         
         rlbwt_data$ratio_mean_misses <- paste(round(rlbwt_data$min_misses/min_misses_vlbt_bwt[1], digits=2),
                                               round(rlbwt_data$min_misses/max_misses_vlbt_bwt[1], digits=2),
                                               sep = " - ")
         rlbwt_data$ratio_mean_space <- paste(round(rlbwt_data$min_space/min_space_vlbt_bwt[1], digits=2),
                                              round(rlbwt_data$min_space/max_space_vlbt_bwt[1], digits=2),
                                              sep = " - ")
         
         mean_misses_vlbt_sri_va <- sort(dt[str_detect(dt$tool, "^vlbt_sri"),4])
         mean_space_vlbt_sri_va <- sort(dt[str_detect(dt$tool, "^vlbt_sri"),7])
         csa_data <- dt[dt$tool %in% c("r_index", "sr-index-va", "vlbt_sri_1", "vlbt_sri_2", "vlbt_sri_3", "vlbt_sri_4"),]
         
         csa_data$ratio_mean_misses <- paste(round(csa_data$mean_misses/mean_misses_vlbt_sri_va[1], digits=4),
                                             round(csa_data$mean_misses/mean_misses_vlbt_sri_va[length(mean_misses_vlbt_sri_va)], digits=4),
                                             sep = " - ")
         csa_data$ratio_mean_space <- paste(round(csa_data$mean_space/mean_space_vlbt_sri_va[1], digits=2),
                                            round(csa_data$mean_space/mean_space_vlbt_sri_va[length(mean_space_vlbt_sri_va)], digits=2),
                                            sep = " - ")
         rbind(rlbwt_data, csa_data)[,c(1,8,9)]
       }
)
count_miss_stats <- ldply(count_miss_stats, rbind)
count_miss_stats <- count_miss_stats[order(count_miss_stats$tool),]

#misses in locate queries in compressed suffix arrays
locate_miss_stats <- tapply(locate_cmiss_csa,
       as.factor(locate_cmiss_csa$input_file),
       function(x){
         stats <- tapply(x,
                         as.factor(x$data_structure),
                         function(y){
                           #this is only to compare sri-indexes
                           y <- y[!(y$parameter %in% c("12_sri","20_sri",
                                                       "12_vm", "20_vm",
                                                       "12_va", "20_va",
                                                       "4096_4", "4096_32",
                                                       "16384_4", "16384_32",
                                                       "65536_4", "65536_32",
                                                       "262144_4","262144_32"
                           )),]
                           #print(y)
                           c(min(y$LD1_data_misses_symbol),
                             max(y$LD1_data_misses_symbol),
                             mean(y$LD1_data_misses_symbol),
                             round(min(y$bits_per_symbol), digits=4),
                             round(max(y$bits_per_symbol), digits=4),
                             round(mean(y$bits_per_symbol), digits=4))
                         }
         )
         dt <- ldply(stats, rbind)
         colnames(dt) <- c("tool", "min_misses", "max_misses", "mean_misses", "min_space", "max_space", "mean_space")
         dt <- dt[dt$tool!="sr-index-sri" & dt$tool!="sr-index-vm",]
         mean_miss_vlbt <- sort(dt[str_detect(dt$tool, "^vlbt_sri_"),4])
         mean_space_vlbt <- sort(dt[str_detect(dt$tool, "^vlbt_sri_"),7])
         
         dt$ratio_mean_miss <- paste(round(dt$mean_misses/mean_miss_vlbt[1], digits=2),
                                     round(dt$mean_misses/mean_miss_vlbt[length(mean_miss_vlbt)], digits=2),
                                     sep = " - ")
         dt$ratio_mean_space <- paste(round(dt$mean_space/mean_space_vlbt[1], digits=2),
                                      round(dt$mean_space/mean_space_vlbt[length(mean_space_vlbt)], digits=2),
                                      sep = " - ")
         dt
       }
)
locate_miss_stats <- ldply(locate_miss_stats, rbind)
locate_miss_stats <- locate_miss_stats[order(locate_miss_stats$tool),]


#space breakdown stats
tapply(vlbt_rlbwt_stacked_space,
       as.factor(vlbt_rlbwt_stacked_space$type),
       function(x){
         x
       }
)

tapply(vlbt_sri_stacked_space,
       as.factor(vlbt_sri_stacked_space$type),
       function(x){
         stats <- tapply(x, as.factor(x$bwt_block_size),
                         function(y){
                           y$sa_space <- y$bwt_SA_heads + y$phi_D_array
                           y$reduction <- y$sa_space-y$sa_space[1]
                           y$others <- y$bwt_tree_data + y$phi_tree_data + y$phi_valid_area
                           y
                         }) 
         stats
       }
)