########################################################################################
#                                    0.Function pre                                   #
#######################################################################################
#source("Function_Source.R")
#######################################################################################
#                                   1.PKG pre                                         #
#######################################################################################
# Included in Function Source
#######################################################################################
#                                  2.Data represent                                   #
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
#
TestingSubsetNames = names(total_genes_all)
#######################################################################################
#                             2. Conversion to EntrezID                               #
#######################################################################################
# Here we have the data compilation
# List with lenght 5 (one for each, easier for looping)
#    sig_genes_all
#    total_genes_all
#
# Now we convert
# two steps - 
# 1. match the records in org.Bt.eg.db
# 2, if not exist, try alias2Symbol() to get new name
# 3. if still fail to match, we have to ignore
     # Potentially lose some, but that's life #

library(org.Bt.eg.db)
# keys return the keys for the database contained in the MeSHdb object
key.symbol = keys(org.Bt.eg.db,  keytype = c("SYMBOL"))
entrezUniverse = select(org.Bt.eg.db, as.character(key.symbol), 
                        columns = c("ENTREZID"),keytype = "SYMBOL") %>% 
  dplyr::distinct(SYMBOL,.keep_all= TRUE)
# dim(entrezUniverse)

#
Sig_list_out = list()
Total_list_out = list()
library(limma)
for (i in c(1:5)){
  ## for sig
  tmp1 = data.frame(SYMBOL = unlist(sig_genes_all[[i]]))
  tmp1 = dplyr::left_join(tmp1 ,entrezUniverse, by = c("SYMBOL" = "SYMBOL"))
  # find the gap
  gather1 = tmp1$SYMBOL
  for (n in seq_along(tmp1$ENTREZID)){
    if (is.na(tmp1$ENTREZID[n])){
      trans = alias2Symbol(tmp1$SYMBOL[n],species = "Bt", expand.symbols = F)[1]
      gather1[n] = trans }
  }
  tmp1 = dplyr::mutate(tmp1,limma_cvt = gather1) %>% dplyr::left_join(entrezUniverse, by = c("limma_cvt" = "SYMBOL"),suffix = c("_orig", "_final"))
  Sig_list_out[[i]] = tmp1
  
  ## for total 
  tmp2 = data.frame(SYMBOL = unlist(total_genes_all[[i]]))
  tmp2 = dplyr::left_join(tmp2,entrezUniverse,by = c("SYMBOL" = "SYMBOL"))
  # find the gap
  gather2 = tmp2$SYMBOL
  for (n in seq_along(tmp2$ENTREZID)){
    if (is.na(tmp2$ENTREZID[n])){
      trans = alias2Symbol(tmp2$SYMBOL[n],species = "Bt", expand.symbols = F)[1]
      gather2[n] = trans }
  }
  tmp2 = dplyr::mutate(tmp2,limma_cvt = gather2) %>% dplyr::left_join(entrezUniverse, by = c("limma_cvt" = "SYMBOL"),suffix = c("_orig", "_final"))
  Total_list_out[[i]] = tmp2
}

##
# Keep only the entrez ID: then we have one vector for each element of the list (some format as always)
Sig_list_out_entrez = list()
Total_list_out_entrez = list()
for (i in c(1:5)){
  Sig_list_out_entrez[[i]] = unique(na.omit(data.frame(Sig_list_out[[i]])$ENTREZID_final))
  names(Sig_list_out_entrez)[i] = names(Sig_list_out)[i]
  Total_list_out_entrez[[i]] = unique(na.omit(data.frame(Total_list_out[[i]])$ENTREZID_final))
  names(Total_list_out_entrez)[i] = names(Total_list_out)[i]
}
## save the convert for next time, it will take some time if run again
save(Total_list_out,Sig_list_out,Sig_list_out_entrez,Total_list_out_entrez,file = "ConvertName2Entrez.RData")
load("ConvertName2Entrez.RData")

# print out the converting details
require(openxlsx)
#write.xlsx(Sig_list_out,file = "test_convert_sig.xlsx")
#write.xlsx(Total_list_out,file = "test_convert_total.xlsx")
#########################################################################################################################

