
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("reshape2", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
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
suppressMessages(library("vroom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("DESeq2", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)


#### Meta analysis code from https://github.com/julirsch/finemapped_mpra/blob/main/code/preprocess/mpra_meta.R ----


meta_analysis_function = function(option_list)
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
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  path_results<-paste(out,'results','/',sep='')
  
  if (file.exists(path_results)){
  }else{
    
    dir.create(file.path(path_results))
  }
  
  path_collapsed_results<-paste(out,'collapsed_results','/',sep='')
  
  if (file.exists(path_collapsed_results)){
  }else{
    
    dir.create(file.path(path_collapsed_results))
  }
  
  #### READ and transform NCGR ----
  
  NCGR = unlist(strsplit(opt$NCGR, split=","))
  
  cat("NCGR_\n")
  cat(sprintf(as.character(NCGR)))
  cat("\n")
  
  #### READ and transform ASE_CTRL ----
  
  ASE_CTRL = unlist(strsplit(opt$ASE_CTRL, split=","))
  
  cat("ASE_CTRL_\n")
  cat(sprintf(as.character(ASE_CTRL)))
  cat("\n")
  
  #### READ and transform Enhancer_CTRL ----
  
  Enhancer_CTRL = unlist(strsplit(opt$Enhancer_CTRL, split=","))
  
  cat("Enhancer_CTRL_\n")
  cat(sprintf(as.character(Enhancer_CTRL)))
  cat("\n")
  
  
  #### READ and transform Threshold_log2FC_meta ----
  
  Threshold_log2FC_meta = opt$Threshold_log2FC_meta
  
  cat("Threshold_log2FC_meta\n")
  cat(sprintf(as.character(Threshold_log2FC_meta)))
  cat("\n")
  
  #### READ and transform Threshold_FC_meta_padj ----
  
  Threshold_FC_meta_padj = opt$Threshold_FC_meta_padj
  
  cat("Threshold_FC_meta_padj\n")
  cat(sprintf(as.character(Threshold_FC_meta_padj)))
  cat("\n")
  
  #### READ and transform Threshold_Skew_meta_padj ----
  
  Threshold_Skew_meta_padj = opt$Threshold_Skew_meta_padj
  
  cat("Threshold_Skew_meta_padj\n")
  cat(sprintf(as.character(Threshold_Skew_meta_padj)))
  cat("\n")
  
  #### Read in all the results of the meta-analisis ----
  
  file_list <- list.files(path=path_results, include.dirs = FALSE)
  
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  
  indexes_sel <- grep("_RNA_emVAR_glm\\.out",file_list)
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  colnames(file_list_sel)<-"file"
  
 
  cat("file_list_sel_0\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$Cell_Type<-gsub("__.+$","",file_list_sel$file)
  
  cat("file_list_sel_0.25\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$comparison<-gsub("_RNA_emVAR_glm\\.out","",file_list_sel$file)
  
  cat("file_list_sel_0.5\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$comparison<-gsub("^[^__]+__","",file_list_sel$comparison)
  
  cat("file_list_sel_0.75\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$comparison<-gsub("diplotypes_[0-9]+_","",file_list_sel$comparison)
  
  cat("file_list_sel_0.95\n")
  cat(str(file_list_sel))
  cat("\n")
  
  ############# LOOP --------------------------
  
  List_RESULTS<-list()
  
  
  # Results_DEF<-data.frame()
  
  DEBUG<-0
  
  for(i in 1:dim(file_list_sel)[1]){
    
    read_file_sel<-file_list_sel$file[i]
    Cell_Type_sel<-file_list_sel$Cell_Type[i]
    comparison_sel<-file_list_sel$comparison[i]
    
    
    cat("-------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(read_file_sel)))
    cat("\t")
    cat(sprintf(as.character(Cell_Type_sel)))
    cat("\t")
    cat(sprintf(as.character(comparison_sel)))
    cat("\n")
    
    setwd(path_results)
    
    
    
    SIZE_gate<-file.info(read_file_sel)$size
    
    if(DEBUG ==1)
    {
      cat("SIZE_gate\n")
      cat(str(SIZE_gate))
      cat("\n")
    }
    
    
    
    if(SIZE_gate> 0)
    {
      
      
      LINE_gate<-length(readLines(read_file_sel))
      
      if(DEBUG ==1)
      {
        cat("LINE_gate\n")
        cat(str(LINE_gate))
        cat("\n")
      }
      
      if(LINE_gate> 0)
      {
        
        
        
        df<-as.data.frame(fread(file=read_file_sel, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
        
        df$KEY<-gsub("__.+","",df$ID)
        
        df$Project<-NA
        
        df$Project<-'MPRA_bc_synthesis'
        
        df$Project[which(df$KEY%in%NCGR)]<-paste('MPRA_bc_synthesis','negCtrl', sep="_")
        df$Project[which(df$KEY%in%ASE_CTRL)]<-paste('MPRA_bc_synthesis','expCtrl_and_emVAR', sep="_")
        df$Project[which(df$KEY%in%Enhancer_CTRL)]<-paste('MPRA_bc_synthesis','expCtrl', sep="_")
        
        df$sample<-'assayed_variant'
        
        df$sample[which(df$KEY%in%NCGR)]<-'negCtrl'
        df$sample[which(df$KEY%in%ASE_CTRL)]<-'expCtrl_and_emVAR'
        df$sample[which(df$KEY%in%Enhancer_CTRL)]<-'expCtrl'
        
        if(DEBUG ==1)
        {
          cat("df_0\n")
          cat(str(df))
          cat("\n")
        }
        
        df$Cell_Type<-Cell_Type_sel
        
        if(DEBUG ==1)
        {
          cat("df_1\n")
          cat(str(df))
          cat("\n")
        }
        
        df$comparison<-NA
        
        indx.diplotypes<-grep("__",df$SNP)
        
        if(DEBUG ==1)
        {
          cat("indx.diplotypes_0\n")
          cat(str(indx.diplotypes))
          cat("\n")
        }
        
        if(length(indx.diplotypes) >0){
          
          df$comparison[indx.diplotypes]<-comparison_sel
          df$comparison[-indx.diplotypes]<-'ref_vs_alt'
          
        }else{
          
          df$comparison<-'ref_vs_alt'
        }#length(indx.diplotypes) >0
        
        
        
        if(DEBUG ==1)
        {
          cat("df_2\n")
          cat(str(df))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(df$comparison))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(df$comparison)))))
          cat("\n")
        }
        
        check<-df[which(df$SNP == 'chr10_103835537_C_T'),]
        
        if(DEBUG ==1)
        {
          
          cat("check_2\n")
          cat(str(check))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(check$comparison))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(check$comparison)))))
          cat("\n")
        }
        
        FLAG_NA<-dim(df[is.na(df$comparison),])[1]
        
        if(DEBUG ==1)
        {
          
          cat("FLAG_NA_0\n")
          cat(str(FLAG_NA))
          cat("\n")
        }
        
        if(FLAG_NA>0){
          
          stop("NAs_in_comparison")
        }#length(FLAG_NA)>0
        
        
        List_RESULTS[[i]]<-df
        
       
        
        
        
      }#LINE_gate
      else{
        
        cat(sprintf("empty_file\n"))
        cat(sprintf(as.character(read_file_sel)))
        cat("\n")
        
        
      }
    }#SIZE_gates
    
  }#i in 1:dim(file_list_sel)[1]
  
  
  Results_redundant = unique(as.data.frame(data.table::rbindlist(List_RESULTS, fill = T)))
  
  cat("Results_redundant_0\n")
  cat(str(Results_redundant))
  cat("\n")
  
  Results_redundant.dt<-data.table(Results_redundant, key=c("SNP","window","Cell_Type","comparison"))
  
  Results_non_redundant<-data.frame(Results_redundant.dt[,.SD[which.max(A_log2FC)], by=key(Results_redundant.dt)], stringsAsFactors = F)
  
  cat("Results_non_redundant_0\n")
  cat(str(Results_non_redundant))
  cat("\n")
  
  Results_non_redundant$Cell_Type<-factor(Results_non_redundant$Cell_Type,
                                          levels=c('K562','CHRF','HL60','THP1'),
                                          ordered=T)
  
  Results_non_redundant$comparison<-factor(Results_non_redundant$comparison,
                                           levels=c('ref_vs_alt','ref_ref_vs_alt_ref','ref_ref_vs_ref_alt','ref_ref_vs_alt_alt'),
                                           ordered=T)
  
  Results_non_redundant<-Results_non_redundant[order(Results_non_redundant$Cell_Type,Results_non_redundant$Cell_Type,Results_non_redundant$comparison),]
  
  cat("Results_non_redundant_1\n")
  cat(str(Results_non_redundant))
  cat("\n")
  
  colnames(Results_non_redundant)[which(colnames(Results_non_redundant) == 'Skew_logFDR')]<-'Skew_logPadj'
  colnames(Results_non_redundant)[which(colnames(Results_non_redundant) == 'Log2Skew')]<-'log2Skew'
  
  cat("Results_non_redundant_2\n")
  cat(str(Results_non_redundant))
  cat("\n")
  
  
  
  
  

  #### collapse A_log2FC and B_log2FC into log2FC, A_log2FC_SE and B_log2FC_SE into log2FC_SE
  
  mutate_tiles<-Results_non_redundant %>%
    ungroup() %>%
    dplyr::select(ID, SNP, comparison, Cell_Type, Project, window, 
                  A_Ctrl_Mean, A_Exp_Mean, A_log2FC, A_log2FC_SE, A_logPadj_BF,
                  B_Ctrl_Mean, B_Exp_Mean, B_log2FC, B_log2FC_SE, B_logPadj_BF,
                  log2Skew, Skew_SE, Skew_logPadj) %>%
    group_by(ID, SNP, comparison, Cell_Type, window) %>%
    dplyr::mutate(log2FC = ifelse(max(abs(A_log2FC) > abs(B_log2FC), na.rm = T), A_log2FC, B_log2FC),
                  log2FC_SE = ifelse(max(abs(A_log2FC) > abs(B_log2FC), na.rm = T), A_log2FC_SE, B_log2FC_SE),
                  logPadj_BF = max(A_logPadj_BF, B_logPadj_BF, na.rm = T)) %>%
    ungroup()
  
  cat("mutate_tiles_0\n")
  cat(str(mutate_tiles))
  cat("\n")
  
  setwd(path_collapsed_results)
  
  write.table(mutate_tiles, 
              file=paste("MPRA_results_prior_to_meta_analysis",".tsv",sep=""), 
              row.names = F, quote=F, sep="\t")
  
  
  ######################## Meta-analysis by SNP, Cell_Type and comparison -----------------------------------------
  
  
  mpra_meta_df_SNP_CT <- mutate_tiles %>%
    ungroup() %>%
    dplyr::select(SNP, Cell_Type, comparison,Project, A_log2FC, A_log2FC_SE, log2FC, log2FC_SE, log2Skew, Skew_SE) %>%
    distinct() %>%
    na.omit() %>%
    dplyr::group_by(SNP,comparison, Cell_Type, Project) %>%
    dplyr::mutate(log2FC_A_meta = sum(A_log2FC * (1 / A_log2FC_SE^2)) / sum(1 / A_log2FC_SE^2),
                  log2FC_meta = sum(log2FC * (1 / log2FC_SE^2)) / sum(1 / log2FC_SE^2),
                  log2FC_meta_SE = sqrt(1 / sum(1 / log2FC_SE^2)),
                  log2Skew_meta =  sum(log2Skew * (1 / Skew_SE^2)) / sum(1 / Skew_SE^2),
                  log2Skew_meta_SE = sqrt(1 / sum(1 / Skew_SE^2)),
                  log2FC_meta_nlogp = -1 * pnorm(abs(log2FC_meta / log2FC_meta_SE), lower.tail = F, log.p = T),
                  log2Skew_meta_nlogp = -1 * pnorm(abs(log2Skew_meta / log2Skew_meta_SE), lower.tail = F, log.p = T)) %>%
    ungroup() %>%
    distinct(SNP,comparison, Cell_Type, Project, log2FC_A_meta, log2FC_meta, log2FC_meta_SE, log2Skew_meta, log2Skew_meta_SE, log2FC_meta_nlogp, log2Skew_meta_nlogp) %>%
    dplyr::mutate(lof2FC_meta_padj = -1 * log10(p.adjust(10^(-1 * log2FC_meta_nlogp), "bonferroni")),
                  lof2Skew_meta_padj = -1 * log10(p.adjust(10^(-1 * log2Skew_meta_nlogp), "fdr"))) %>%
    dplyr::mutate(active_meta = ifelse(lof2FC_meta_padj >= -log10(Threshold_FC_meta_padj) & abs(log2FC_meta) >= Threshold_log2FC_meta, T, F),
                  emVar_meta = ifelse(active_meta & lof2Skew_meta_padj >= -log10(Threshold_Skew_meta_padj) & !is.na(lof2Skew_meta_padj) & abs(log2Skew_meta) >= 0, T, F))

  cat("mpra_meta_df_SNP_CT_0\n")
  cat(str(mpra_meta_df_SNP_CT))
  cat("\n")

  mpra_meta_df_SNP_CT$Activity<-'NA'


  mpra_meta_df_SNP_CT$Activity[which(mpra_meta_df_SNP_CT$active_meta == 'TRUE')]<-'ACTIVE'
  mpra_meta_df_SNP_CT$Activity[which(mpra_meta_df_SNP_CT$active_meta == 'FALSE')]<-'INACTIVE'

  mpra_meta_df_SNP_CT$Activity<-factor(mpra_meta_df_SNP_CT$Activity,
                                         levels=c('INACTIVE','ACTIVE'),
                                         ordered=T)

  mpra_meta_df_SNP_CT$MPRA_CLASS<-'NA'


  mpra_meta_df_SNP_CT$MPRA_CLASS[which(mpra_meta_df_SNP_CT$active_meta == 'TRUE' & mpra_meta_df_SNP_CT$emVar_meta == 'TRUE')]<-'MPRA_positive'
  mpra_meta_df_SNP_CT$MPRA_CLASS[-which(mpra_meta_df_SNP_CT$active_meta == 'TRUE' & mpra_meta_df_SNP_CT$emVar_meta == 'TRUE')]<-'MPRA_negative'

  mpra_meta_df_SNP_CT$MPRA_CLASS<-factor(mpra_meta_df_SNP_CT$MPRA_CLASS,
                                                               levels=c('MPRA_negative','MPRA_positive'),
                                                               ordered=T)

  cat("mpra_meta_df_SNP_CT_1\n")
  cat(str(mpra_meta_df_SNP_CT))
  cat("\n")

  setwd(path_collapsed_results)

  write.table(mpra_meta_df_SNP_CT,
              file=paste("MPRA_results_meta_analysis_collapsed_by_SNP_and_Cell_Type_Threshold_log2FC_",Threshold_log2FC_meta,"_Threshold_FC_meta_padj_",Threshold_FC_meta_padj,".tsv",sep=""),
              row.names = F, quote=F, sep="\t")

  saveRDS(mpra_meta_df_SNP_CT,
              file=paste("MPRA_results_meta_analysis_collapsed_by_SNP_and_Cell_Type_Threshold_log2FC_",Threshold_log2FC_meta,"_Threshold_FC_meta_padj_",Threshold_FC_meta_padj,".rds",sep=""))


  
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
    make_option(c("--Threshold_log2FC_meta"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_FC_meta_padj"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_Skew_meta_padj"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--NCGR"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_CTRL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Enhancer_CTRL"), type="character", default=NULL, 
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
  
  meta_analysis_function(opt)
    
  
}

###########################################################################

system.time( main() )