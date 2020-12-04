## prepare the following files for Step3 test
# demo/HIPPO_demo_20201204.RData

###################### demo script to test Step3
source('src/pipeline_function.R')
######################
load('demo/HIPPO_demo_20201204.RData')
######################
output_pdf  <- 'demo/demo.pdf'
output_file <- 'demo/demo.txt'
######################
res <- HIPPO.prepare_dataset(mat=demo_mat,mod_region=mod_region,
                             mod_info=mod_info,int_pos=int_pos)
all_mat <- res$all_mat
group_info <- res$group_info
group_name <- res$group_name
pre_define <- res$pre_define
###################### haplotype imputation
res_trace <- HIPPO.impute_Haplotype(all_mat=all_mat,
                                    group_info=group_info,
                                    top_each=3)
########################### plot output
pdf(output_pdf,width=5,height = 4)

HIPPO.plot_imputeHaplotype(res_trace,prob_thre=0.01,
                           group_name=group_name)
HIPPO.plot_adjacentCombination(all_mat,remove_char='.',
                               group_info=group_info,
                               group_name=group_name,only_use_group=FALSE,
                               count_link_thre=0,part_thre=0,
                               top_n=2)
HIPPO.plot_adjacentCombination(all_mat,remove_char='.',
                               group_info=group_info,
                               group_name=group_name,only_use_group=FALSE,
                               count_link_thre=0,part_thre=0,
                               top_n=2,strategy = 'count')
HIPPO.plot_readsComponent(all_mat,group_info=group_info,
                          remove_char='.',
                          group_name=group_name)
dev.off()
########################### text output
res_readsComponent <- HIPPO.summ_readsComponent(all_mat,group_info,group_name)
res_finalHaplotype <- res_trace[[max(length(res_trace))]]
HIPPO.output_result(res_readsComponent=res_readsComponent,
                    res_finalHaplotype=res_finalHaplotype,
                    group_name=group_name,
                    output_file=output_file)
