
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("VennDiagram", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("nVennR", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)





data_wrangling_SNP_by_CT = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
 
  
  #### Read MPRA_alignment_results----
  
  MPRA_alignment_results<-readRDS(file=opt$MPRA_alignment_results)
  
  cat("MPRA_alignment_results_0\n")
  cat(str(MPRA_alignment_results))
  cat("\n")
  cat(str(unique(MPRA_alignment_results$VAR)))
  cat("\n")
  
  MPRA_alignment_results_subset<-unique(MPRA_alignment_results[,which(colnames(MPRA_alignment_results)%in%c("Oligo_ID_ALT","SNP","comparison"))])
  
  cat("MPRA_alignment_results_subset_0\n")
  cat(str(MPRA_alignment_results_subset))
  cat("\n")
  cat(str(unique(MPRA_alignment_results_subset$VAR)))
  cat("\n")
  
  colnames(MPRA_alignment_results_subset)[which(colnames(MPRA_alignment_results_subset) == 'Oligo_ID_ALT')]<-'VAR'
  
  MPRA_alignment_results_subset$VAR<-gsub("^.+TILE_","",MPRA_alignment_results_subset$VAR)
  MPRA_alignment_results_subset$VAR<-gsub("^[^_]+__","",MPRA_alignment_results_subset$VAR)
  
  
  cat("MPRA_alignment_results_subset_1\n")
  cat(str(MPRA_alignment_results_subset))
  cat("\n")
  
  #### READ MPRA_S_result ----
  
  MPRA_S_result_by_CT<-as.data.frame(readRDS(file=opt$MPRA_result_SNP_and_CT), stringsAsFactors=F)
  
  
  cat("MPRA_S_result_by_CT_0\n")
  cat(str(MPRA_S_result_by_CT))
  cat("\n")
  
  
  MPRA_S_result_by_CT$Project<-factor(MPRA_S_result_by_CT$Project,
                                levels=c("MPRA_bc_synthesis","MPRA_bc_synthesis_expCtrl","MPRA_bc_synthesis_expCtrl_and_emVAR","MPRA_bc_synthesis_negCtrl"),
                                ordered=T)
  
  
  
 
  
 
  
  ######################### SAVE #3 --------------------------------------------------------
  
  MPRA_S_result_by_CT<-merge(MPRA_S_result_by_CT,
                       MPRA_alignment_results_subset,
                       by=c('SNP','comparison'),
                       all.x=T)
  
  
  cat("MPRA_S_result_by_CT_0\n")
  cat(str(MPRA_S_result_by_CT))
  cat("\n")
  cat(str(unique(MPRA_S_result_by_CT$VAR)))
  cat("\n")
  
  
  setwd(out)
  
  
  write.table(MPRA_S_result_by_CT, file=paste("MPRA_results_meta_analysis_collapsed_by_SNP_and_Cell_Type_Threshold_log2FC_0.255_Threshold_FC_meta_padj_0.01","_VAR_added",".tsv",sep=''),sep="\t",quote=F,row.names = F)
  saveRDS(MPRA_S_result_by_CT,file=paste("MPRA_results_meta_analysis_collapsed_by_SNP_and_Cell_Type_Threshold_log2FC_0.255_Threshold_FC_meta_padj_0.01","_VAR_added",".rds",sep=''))
  write.table(MPRA_S_result_by_CT, file=paste("Table_S4_MPRA_results_meta_analysis_collapsed_by_SNP_and_Cell_Type_Threshold_log2FC_0.255_Threshold_FC_meta_padj_0.01","_VAR_added",".tsv",sep=''),sep="\t",quote=F,row.names = F)
  
  
  
}

data_wrangling_NEW_TABLE_S6 = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ Table_S6 ----
  
  Table_S6<-readRDS(file=opt$Table_S6)
  
  
  # Table_S6<-droplevels(Table_S6[-which(Table_S6$Mechanistic_Class == 'No_RNA_Seq_HET_carriers'),])
  
  Table_S6$Manual_curation<-revalue(Table_S6$Manual_curation,
                                    c("R in candidate" = "GeneBass/lit. supported"))
  
  cat("Table_S6_0\n")
  cat(str(Table_S6))
  cat("\n")
  cat(str(unique(Table_S6$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$genIE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$Manual_curation))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$Multi_Lineage))))
  cat("\n")
  
  
  check.VARs<-Table_S6[which(Table_S6$VAR%in%c("chr20_55990370_A_T","chr15_65174438_C_A","chr15_65174494_A_G")),]
  
  cat("check.VARs_0\n")
  cat(str(check.VARs))
  cat("\n")
  
  
  #### READ MPRA_result_SNP_and_CT ----
  
  setwd(out)
  
  filename<-paste("MPRA_results_meta_analysis_collapsed_by_SNP_and_Cell_Type_Threshold_log2FC_0.255_Threshold_FC_meta_padj_0.01","_VAR_added",".rds",sep='')
  
  MPRA_result_SNP_and_CT<-as.data.frame(readRDS(file=filename), stringsAsFactors=F)
  
  
  cat("MPRA_result_SNP_and_CT_0\n")
  cat(str(MPRA_result_SNP_and_CT))
  cat("\n")
  
  
  MPRA_result_SNP_and_CT_positive<-unique(MPRA_result_SNP_and_CT[which(MPRA_result_SNP_and_CT$MPRA_CLASS == 'MPRA_positive'),])
  
  
  cat("MPRA_result_SNP_and_CT_positive_0\n")
  cat(str(MPRA_result_SNP_and_CT_positive))
  cat("\n")
  cat(str(unique(MPRA_result_SNP_and_CT_positive$VAR)))
  cat("\n")
  
  
  MPRA_result_SNP_and_CT_positive.dt<-data.table(MPRA_result_SNP_and_CT_positive, key=c("VAR","MPRA_CLASS"))
  
  
  MPRA_result_SNP_and_CT_positive_collapsed<-as.data.frame(MPRA_result_SNP_and_CT_positive.dt[,.(MPRA_positive_cells=unique(paste(unique(Cell_Type), collapse = '|'))), 
                                                                                              by=key(MPRA_result_SNP_and_CT_positive.dt)], stringsAsFactors=F)
  
  
  cat("MPRA_result_SNP_and_CT_positive_collapsed_0\n")
  cat(str(MPRA_result_SNP_and_CT_positive_collapsed))
  cat("\n")
  cat(str(unique(MPRA_result_SNP_and_CT_positive_collapsed$VAR)))
  cat("\n")
  
  
  restricted_MPRA_result_SNP_and_CT_positive_collapsed<-MPRA_result_SNP_and_CT_positive_collapsed[which(MPRA_result_SNP_and_CT_positive_collapsed$VAR%in%Table_S6$VAR),]
  
  
  cat("restricted_MPRA_result_SNP_and_CT_positive_collapsed_0\n")
  cat(str(restricted_MPRA_result_SNP_and_CT_positive_collapsed))
  cat("\n")
  cat(str(unique(restricted_MPRA_result_SNP_and_CT_positive_collapsed$VAR)))
  cat("\n")
  
  
  # setwd(out)
  # 
  # 
  # write.table(restricted_MPRA_result_SNP_and_CT_positive_collapsed, file='test.tsv', sep="\t",quote=F, row.names = F)
  
  
 
  
  #### NEW TABLE S6----
  
  NEW_Table_S6<-merge(restricted_MPRA_result_SNP_and_CT_positive_collapsed,
                      Table_S6[,-which(colnames(Table_S6) == 'MPRA_CLASS')],
                      by=c('VAR'),
                      all.y=T)
  
  
  NEW_Table_S6$MPRA_CLASS[is.na(NEW_Table_S6$MPRA_CLASS)]<-'MPRA_negative'
  
  
  NEW_Table_S6$MPRA_CLASS<-factor(NEW_Table_S6$MPRA_CLASS,
                                      levels=c("MPRA_negative","MPRA_positive"),
                                      ordered=T)
  
  
  
  
  cat("NEW_Table_S6_0\n")
  cat(str(NEW_Table_S6))
  cat("\n")
  cat(str(unique(NEW_Table_S6$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(NEW_Table_S6$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(NEW_Table_S6$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(NEW_Table_S6$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(NEW_Table_S6$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(NEW_Table_S6$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(summary(NEW_Table_S6$Manual_curation))))
  cat("\n")
  cat(sprintf(as.character(names(summary(NEW_Table_S6$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(NEW_Table_S6$Multi_Lineage))))
  cat("\n")
  
  
  
  check.VARs<-NEW_Table_S6[which(NEW_Table_S6$VAR%in%c("chr20_55990370_A_T","chr15_65174438_C_A","chr15_65174494_A_G")),]
  
  cat("check.VARs_3\n")
  cat(str(check.VARs))
  cat("\n")
  
  
  levels_MPRA_CLASS<-rev(levels(NEW_Table_S6$MPRA_CLASS))
  
  cat("levels_MPRA_CLASS\n")
  cat(str(levels_MPRA_CLASS))
  cat("\n")
  
  NEW_Table_S6$MPRA_CLASS<-factor(NEW_Table_S6$MPRA_CLASS,
                                  levels=levels_MPRA_CLASS,
                                  ordered=T)
  
  
  
  levels_Mechanistic_Class<-levels(NEW_Table_S6$Mechanistic_Class)
  
  cat("levels_Mechanistic_Class\n")
  cat(str(levels_Mechanistic_Class))
  cat("\n")
  
  levels_Manual_curation<-levels(NEW_Table_S6$Manual_curation)
  
  cat("levels_Manual_curation\n")
  cat(str(levels_Manual_curation))
  cat("\n")
  
  levels_Multi_Lineage<-levels(NEW_Table_S6$Multi_Lineage)
  
  cat("levels_Multi_Lineage\n")
  cat(str(levels_Multi_Lineage))
  cat("\n")
  
  pool_levels<-levels_Mechanistic_Class
  
  cat("pool_levels_0\n")
  cat(str(pool_levels))
  cat("\n")
  
  
  NEW_Table_S6$Mechanistic_Class_compressed<-NA
  
  NEW_Table_S6$Mechanistic_Class_compressed[which(NEW_Table_S6$Mechanistic_Class%in%pool_levels[c(1:3)])]<-'DE and or ATU'
  NEW_Table_S6$Mechanistic_Class_compressed[which(NEW_Table_S6$Mechanistic_Class%in%pool_levels[4])]<-pool_levels[4]
  
  NEW_Table_S6$Mechanistic_Class_compressed<-factor(NEW_Table_S6$Mechanistic_Class_compressed,
                                                    levels=c('DE and or ATU',pool_levels[4]),
                                                    ordered=T)
  
  cat(sprintf(as.character(names(summary(NEW_Table_S6$Mechanistic_Class_compressed)))))
  cat("\n")
  cat(sprintf(as.character(summary(NEW_Table_S6$Mechanistic_Class_compressed))))
  cat("\n")
  
  
  
  
  NEW_Table_S6$interaction<-droplevels(interaction(NEW_Table_S6$MPRA_CLASS,NEW_Table_S6$Mechanistic_Class_compressed, sep='|'))
  
  
  pool_levels_interaction<-levels(NEW_Table_S6$interaction)
  
  cat("NEW_Table_S6$interaction_0\n")
  cat(str(pool_levels_interaction))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(NEW_Table_S6$interaction)))))
  cat("\n")
  cat(sprintf(as.character(summary(NEW_Table_S6$interaction))))
  cat("\n")
  
  
  NEW_Table_S6_subset<-unique(NEW_Table_S6[,which(colnames(NEW_Table_S6)%in%c("VAR","rs","Mechanistic_Class","Manual_curation","Candidate_effector","OpenTargets_QTL","Whole_blood_DE_HGNC_string","Whole_blood_DTU_HGNC_string","Monocyte_DE_HGNC_string","Tcell_DE_HGNC_string","Neutrophil_DE_HGNC_string",
                                                                              "Neutrophil_DTU_HGNC_string","Tcell_DTU_HGNC_stringMonocyte_DTU_HGNC_string","Proxy_rsid_string","Replication_OT_QTL","VAR_38","Multi_Lineage","Lineage_string","phenotype_DEF_string","maf_origin",
                                                                              "VEP_DEF_LABELS_wCSQ","integration_category","Activity","MPRA_CLASS","MPRA_positive_cells","Mechanistic_Class_compressed","interaction"))])
  
  cat("NEW_Table_S6_subset_0\n")
  cat(str(NEW_Table_S6_subset))
  cat("\n")
  cat(str(unique(NEW_Table_S6_subset$VAR)))
  cat("\n")
  
  
  
  ########################## SAVE # 1 ------------------------------------------------------------------------
  
  setwd(out)
  
  
  write.table(NEW_Table_S6_subset, file='NEW_Table_S6.tsv',sep="\t",quote=F,row.names = F)
  saveRDS(NEW_Table_S6_subset,file="NEW_Table_S6.rds")
  
  
}

data_wrangling_bed_table = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  
  #### Read MPRA_alignment_results----
  
  MPRA_alignment_results<-readRDS(file=opt$MPRA_alignment_results)
  
  cat("MPRA_alignment_results_0\n")
  cat(str(MPRA_alignment_results))
  cat("\n")
  cat(str(unique(MPRA_alignment_results$VAR)))
  cat("\n")
  
  MPRA_alignment_results_subset<-unique(MPRA_alignment_results[,which(colnames(MPRA_alignment_results)%in%c("Oligo_ID_ALT","SNP","comparison","window","chr","start","end"))])
  
  cat("MPRA_alignment_results_subset_0\n")
  cat(str(MPRA_alignment_results_subset))
  cat("\n")
  cat(str(unique(MPRA_alignment_results_subset$VAR)))
  cat("\n")
  
  colnames(MPRA_alignment_results_subset)[which(colnames(MPRA_alignment_results_subset) == 'Oligo_ID_ALT')]<-'VAR'
  
  
  MPRA_alignment_results_subset$VAR<-gsub("^.+TILE_","",MPRA_alignment_results_subset$VAR)
  MPRA_alignment_results_subset$VAR<-gsub("^[^_]+__","",MPRA_alignment_results_subset$VAR)
  
  cat("MPRA_alignment_results_subset_1\n")
  cat(str(MPRA_alignment_results_subset))
  cat("\n")
  
  
  #### READ MPRA_S_result ----
  
  MPRA_prior_to_meta_analysis<-as.data.frame(fread(file=opt$MPRA_prior_to_meta_analysis, sep="\t",header=T), stringsAsFactors=F)
  
  
  cat("MPRA_prior_to_meta_analysis_0\n")
  cat(str(MPRA_prior_to_meta_analysis))
  cat("\n")
  

  
  MPRA_prior_to_meta_analysis$Cell_Type<-factor(MPRA_prior_to_meta_analysis$Cell_Type,
                                              levels=c("K562","CHRF","HL60","THP1"),
                                              ordered=T)
  
  
  MPRA_prior_to_meta_analysis$Project<-factor(MPRA_prior_to_meta_analysis$Project,
                                      levels=c("MPRA_bc_synthesis","MPRA_bc_synthesis_expCtrl","MPRA_bc_synthesis_expCtrl_and_emVAR","MPRA_bc_synthesis_negCtrl"),
                                      ordered=T)
  
  
  
  cat("MPRA_prior_to_meta_analysis_1\n")
  cat(str(MPRA_prior_to_meta_analysis))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_prior_to_meta_analysis$Cell_Type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_prior_to_meta_analysis$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_prior_to_meta_analysis$Project))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_prior_to_meta_analysis$Project)))))
  cat("\n")
  
  ######################### SAVE #3 --------------------------------------------------------
  
  MPRA_prior_to_meta_analysis<-merge(MPRA_alignment_results_subset,
                                     MPRA_prior_to_meta_analysis,
                             by=c('SNP','comparison','window'),
                             all=T)
  
  
  cat("MPRA_prior_to_meta_analysis_2\n")
  cat(str(MPRA_prior_to_meta_analysis))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_prior_to_meta_analysis$Cell_Type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_prior_to_meta_analysis$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_prior_to_meta_analysis$Project))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_prior_to_meta_analysis$Project)))))
  cat("\n")
  
  
  setwd(out)
  
  
  write.table(MPRA_prior_to_meta_analysis, file=paste("MPRA","_prior_to_metanalysis",".bed",sep=''),sep="\t",quote=F,row.names = F)
  write.table(MPRA_prior_to_meta_analysis, file=paste("Table_S3_MPRA_result_per_tile_and_cell_type",".tsv",sep=''),sep="\t",quote=F,row.names = F)
  
  
  
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--MPRA_alignment_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_result_SNP_and_CT"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_prior_to_meta_analysis"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Table_S6"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  
  data_wrangling_SNP_by_CT(opt)
  data_wrangling_NEW_TABLE_S6(opt)
  data_wrangling_bed_table(opt)

  
  
}


###########################################################################

system.time( main() )