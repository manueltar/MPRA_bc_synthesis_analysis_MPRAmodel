## 1. Get the regular matrix, add by barcode and normalize

$ bash ~/Scripts/Wraper_scripts/95_MPRA_bc_synthesis_MPRAmodel_v3_PartONE.sh /group/soranzo/manuel.tardaguila/MPRA_bc_synthesis/ Analysis_filtered_Plus_Pseudocount 0.255 0.01

## 2. Perform MPRA model and meta-analysis

$ bash ~/Scripts/Wraper_scripts/95_MPRA_bc_synthesis_MPRAmodel_v3_PartTWO.sh /group/soranzo/manuel.tardaguila/MPRA_bc_synthesis/ Analysis_filtered_Plus_Pseudocount 0.255 0.01

## 3. Generate the final files and perform lasso regression

$ bash ~/Scripts/Wraper_scripts/114_paper_relaunch_MPRA_v3_lasso_only_cell_matched.sh /group/soranzo/manuel.tardaguila/MPRA_bc_synthesis/test/Analysis_filtered_Plus_Pseudocount/final_files/ MPRA_S_no_global_metanalysis /group/soranzo/manuel.tardaguila/MPRA_bc_synthesis/test/Analysis_filtered_Plus_Pseudocount/collapsed_results/MPRA_results_meta_analysis_collapsed_by_SNP_and_Cell_Type_Threshold_log2FC_0.255_Threshold_FC_meta_padj_0.01.rds /group/soranzo/manuel.tardaguila/MPRA_explore_results/Alingment_check/MPRA_bc_synthesis_check_reconstruct_from_alignment.rds /group/soranzo/manuel.tardaguila/MPRA_bc_synthesis/test/Analysis_filtered_Plus_Pseudocount/collapsed_results/MPRA_results_prior_to_meta_analysis.tsv


-------------> error: Log_files/outfile_11_plot_variables.log


 NULL

Error in rbind(deparse.level, ...) :
  numbers of columns of arguments do not match
Calls: system.time -> main -> violin_Activity -> rbind -> rbind


######################## SAVE everything when fixed and repush to github
