########################################################################################
#                                    0.Function pre                                   #
#######################################################################################
rm(list = ls())
source("Function_Source.R")
#######################################################################################
#                                   1.PKG pre                                         #
#######################################################################################
library(readxl);library(ggplot2);library(biomaRt);library(tidyverse)
library(biomaRt)
#######################################################################################
#                                  2.Over represent                                   #
#######################################################################################
# data pre 1 : imput (read in)
FPM_CNTRL_raw <- read_excel("gene_exp_FMP_CNTRL_bam4_all.xlsx")
AR_CNTRL_raw <- read_excel("gene_exp_AR_vs_CNTRL_bam2_all.xlsx")
PRF_CNTRL_raw <- read_excel("gene_exp_PRF_vs_CNTRL_April4_UPDATES.xlsx")
SMP_CNTRL_raw <- read_excel("gene_exp_SMP_CNTRL_bam3_all.xlsx")
SMP_FMP_raw <- read_excel("gene_exp_SMP_vs_FMP_bam5_all.xlsx")
raw_data_all_index = c("FPM_CNTRL_raw","AR_CNTRL_raw","PRF_CNTRL_raw","SMP_CNTRL_raw","SMP_FMP_raw")
# data pre 2 : selection
# massage names for continuty
data_all_index = character() # will be containing data - all "OK" and "qval <= .1"
for (i in seq_along(raw_data_all_index)){
  data_all_index[i] = (paste(substring(raw_data_all_index[i],1,nchar(raw_data_all_index[i])-4)))
}

###
# create two containers for double looping, exhaust combanation of dataset - GO repo
total_genes_all = list(FPM_CNTRL = c(),AR_CNTRL = c(),PRF_CNTRL = c(),SMP_CNTRL = c(),SMP_FMP = c())
sig_genes_all   = list(FPM_CNTRL = c(),AR_CNTRL = c(),PRF_CNTRL = c(),SMP_CNTRL = c(),SMP_FMP = c())
for (i in seq_along(raw_data_all_index)){
  tmp1 <- dplyr::select(get(raw_data_all_index[i]),gene,status,q_value) %>% 
    dplyr::filter(status =="OK") %>% group_by(gene) %>% 
    distinct()
  total_genes_all[[i]] = dplyr::select(tmp1,gene)
  
  tmp2 <- dplyr::filter(tmp1,q_value <= 0.1) %>% group_by(gene) %>% 
    distinct()
  sig_genes_all[[i]] = dplyr::select(tmp2,gene)
  assign(data_all_index[i],tmp2)
}

# sig_genes_all[1];total_genes_all[1]
# look into the compilaton; all good!
#######################################################################################
#                                   3. Run test                                       #
#######################################################################################
# double looping
TestingSubsetNames = names(total_genes_all)

Enrich_Results_thres005 = Go_Enrich_Plot(total_genes_all,
                                         sig_genes_all,
                                         TestingSubsetNames,
                                         GOthres = 0.05,
                                         keyword = "GO_Enrichment_qval01_pval005_test3")

Enrich_Results_thres001 = Go_Enrich_Plot(total_genes_all,
                                         sig_genes_all,
                                         TestingSubsetNames,
                                         GOthres = 0.01,
                                         keyword = "GO_Enrichment_qval01_pval001")

