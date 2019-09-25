########################################################################################
#                                    0.Function pre                                   #
#######################################################################################
#rm(list = ls())
source("Function_Source.R")
#options("scipen"= -100, "digits"=4)
#######################################################################################
#                                   1.PKG pre                                         #
#######################################################################################
# Included in Function Source
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
#                                   3. Run analysis - double loop                     #
#######################################################################################
# double looping
TestingSubsetNames = names(total_genes_all)
Enrich_Results_thres005 = Go_Enrich_Plot(total_genes_all,
                                         sig_genes_all,
                                         TestingSubsetNames,
                                         GOthres = 0.05,
                                         keyword = "GO_Enrichment_qval01_pval005_test0924")

#######################################################################################
#                                   4. Join datatest                                 #
#######################################################################################
# load results 
#list.files()
load("GO_Enrichment_qval01_pval005_0921.RData")

# get the spavcename index 
biomart="ensembl";dataset="btaurus_gene_ensembl";attributes = c("go_id","namespace_1003")
database = useMart(biomart);genome = useDataset(dataset, mart = database);gene = getBM(attributes,mart = genome)
namespace_index = dplyr::filter(gene,go_id != "",namespace_1003 != "")

# parse 1 by 1, and attach the space name in the end
compile_select_index = c("go_id","GO_Name","Total_Genes","Significant_Genes","pvalue","hitsPerc")
### group 1
#
AR_CNTRL_enrich = Parse_GO_Results(GO_results_b[2])
names(AR_CNTRL_enrich) = c("go_id","GO_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
AR_CNTRL_enrich = dplyr::select(AR_CNTRL_enrich,compile_select_index) %>% dplyr::inner_join(namespace_index,by = "go_id")
#
PRF_CNTRL_enrich = Parse_GO_Results(GO_results_b[3])
names(PRF_CNTRL_enrich) = c("go_id","GO_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
PRF_CNTRL_enrich   = dplyr::select(PRF_CNTRL_enrich,compile_select_index) %>% dplyr::inner_join(namespace_index,by = "go_id")
### group 2
#
FPM_CNTRL_enrich = Parse_GO_Results(GO_results_b[1])
names(FPM_CNTRL_enrich) = c("go_id","GO_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
FPM_CNTRL_enrich = dplyr::select(FPM_CNTRL_enrich,compile_select_index) %>% dplyr::inner_join(namespace_index,by = "go_id")
#
SMP_CNTRL_enrich= Parse_GO_Results(GO_results_b[4])
names(SMP_CNTRL_enrich) = c("go_id","GO_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
SMP_CNTRL_enrich = dplyr::select(SMP_CNTRL_enrich,compile_select_index) %>% dplyr::inner_join(namespace_index,by = "go_id")
#
SMP_FMP_enrich = Parse_GO_Results(GO_results_b[5])
names(SMP_FMP_enrich) = c("go_id","GO_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
SMP_FMP_enrich =  dplyr::select(SMP_FMP_enrich,compile_select_index) %>% dplyr::inner_join(namespace_index,by = "go_id")

#######################################################################################
#                                   5. Compilation                                    #
#######################################################################################
#### group 1 - regression
GO_Results_full_005_reg <-
  dplyr::full_join(AR_CNTRL_enrich,PRF_CNTRL_enrich,by = c("go_id" = "go_id")) %>% 
  tidyr::replace_na(list(go_id.x = "Not Found",go_id.y = "Not Found")) 
GO_Results_inner_005_reg <-  
  dplyr::inner_join(AR_CNTRL_enrich,PRF_CNTRL_enrich,by = c("go_id" = "go_id"))
#### group 2 - pregnancy
GO_Results_full_005_preg <-  
  dplyr::full_join(FPM_CNTRL_enrich,SMP_CNTRL_enrich,by = c("go_id" = "go_id")) %>% 
  dplyr::full_join(SMP_FMP_enrich,by = c("go_id" = "go_id"))%>% 
  tidyr::replace_na(list(go_id.x = "Not Found",go_id.y = "Not Found",go_id = "Not Found")) 
GO_Results_inner_005_preg <-  
  dplyr::inner_join(FPM_CNTRL_enrich,SMP_CNTRL_enrich,by = c("go_id" = "go_id")) %>% 
  dplyr::inner_join(SMP_FMP_enrich,by = c("go_id" = "go_id"))
require(openxlsx)
GO_Enrich_Regression_005 <- list("Full_join" = GO_Results_full_005_reg, "Inner_join" = GO_Results_inner_005_reg)
write.xlsx(GO_Enrich_Regression_005,file = "GO_Enrich_Regression_005_0921.xlsx")
GO_Enrich_Pregnancy_005 <- list("Full_join" = GO_Results_full_005_preg, "Inner_join" = GO_Results_inner_005_preg)
write.xlsx(GO_Enrich_Pregnancy_005,file = "GO_Enrich_Pregnancy_005_0921.xlsx")

#######################################################################################
#                           6.Take 0.01 - be more strengent                           #
#######################################################################################
str(AR_CNTRL_enrich)
# Function takes Enrichment table as imput (same format as "AR_CNTRL_enrich", for example, like the one above)
ReduceDim_GO_Plot = function(Enrich_Out,
                             GOthres = 0.001,
                             label_size1 = 0.4,
                             label_size2 = 0.4,
                             label_size3 = 0.4,
                             Database = "org.Bt.eg.db",
                             measure="Jiang",combine=NULL,
                             Dataset_Name){
  # load libraries + download ref database
  library(GOSemSim);library(corrplot);library(tidyverse)
  do.call(library,list(Database))
  semData_BP <- godata(paste(Database), ont="BP", computeIC=T)
  semData_MF <- godata(paste(Database), ont="MF", computeIC=T)
  semData_CC <- godata(paste(Database), ont="CC", computeIC=T)
  # selection + formating: for each category we have one vector containing all the sig GO terms
  BP_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "biological_process") %>% 
    dplyr::select(go_id) %>% unlist();attributes(BP_List) = NULL # name is an attribute and we dont them, so set null
  CC_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "cellular_component") %>% 
    dplyr::select(go_id) %>% unlist();attributes(CC_List) = NULL
  MF_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "molecular_function") %>% 
    dplyr::select(go_id) %>% unlist();attributes(MF_List) = NULL
  ### Now we are trying to get all similarity matrix ready. N x N, symetric, diag = 1
  # For BP
  goSimMatrix_BP = GOSemSim::mgoSim(BP_List,
                                    BP_List,
                                    semData=semData_BP,measure=measure,combine = combine)
  suspectID_BP = rownames(goSimMatrix_BP)[is.na(goSimMatrix_BP[,1])]
  if (length(suspectID_BP) != 0){BP_List_new = BP_List[-which(BP_List == suspectID)]
  message(length(suspectID_BP)," invalid ID captured in BP: ",suspectID_BP)
  } else {BP_List_new = BP_List;message("Nice! All IDs are valid in BP!")}
  goSimMatrix_BP_new = GOSemSim::mgoSim(BP_List_new,
                                        BP_List_new,
                                        semData=semData_BP,measure=measure,combine = combine)
  colnames(goSimMatrix_BP_new) = paste(BP_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)])
  rownames(goSimMatrix_BP_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)],BP_List_new)
  # For CC
  goSimMatrix_CC = GOSemSim::mgoSim(CC_List,
                                    CC_List,
                                    semData=semData_CC,measure=measure,combine = combine)
  suspectID_CC = rownames(goSimMatrix_CC)[is.na(goSimMatrix_CC[,1])]
  if (length(suspectID_CC) != 0){CC_List_new = CC_List[-which(CC_List == suspectID_CC)]
  message(length(suspectID_CC)," invalid ID captured in CC: ",suspectID_CC)
  } else {CC_List_new = CC_List;message("Nice! All IDs are valid in CC!")}
  goSimMatrix_CC_new = GOSemSim::mgoSim(CC_List_new,
                                        CC_List_new,
                                        semData=semData_CC,measure=measure,combine =combine)
  colnames(goSimMatrix_CC_new) = paste(CC_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)])
  rownames(goSimMatrix_CC_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)],CC_List_new)
  # For MF
  goSimMatrix_MF = GOSemSim::mgoSim(MF_List,
                                    MF_List,
                                    semData=semData_MF,measure=measure,combine = combine)
  suspectID_MF = rownames(goSimMatrix_MF)[is.na(goSimMatrix_MF[,1])]
  if (length(suspectID_MF) != 0){MF_List_new = MF_List[-which(MF_List == suspectID_MF)]
  message(length(suspectID_MF)," invalid ID captured in MF: ",suspectID_MF)
  } else {MF_List_new = MF_List;message("Nice! All IDs are valid in MF!")}
  goSimMatrix_MF_new = GOSemSim::mgoSim(MF_List_new,
                                        MF_List_new,
                                        semData=semData_MF,measure=measure,combine = combine)
  colnames(goSimMatrix_MF_new) = paste(MF_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)])
  rownames(goSimMatrix_MF_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)],MF_List_new)
  # Now we take the results and plot
  pdf(paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".pdf",sep = ""))
  corrplot(goSimMatrix_CC_new,title = "Semantic_Similarity_Measure_CC",
           tl.col = "black", tl.cex = label_size1, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  corrplot(goSimMatrix_BP_new,title = "Semantic_Similarity_Measure_BP",
           tl.col = "black", tl.cex = label_size2, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  corrplot(goSimMatrix_MF_new,title = "Semantic_Similarity_Measure_MF",
           tl.col = "black", tl.cex = label_size3, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  dev.off()
  message(dim(goSimMatrix_CC_new)[1],",",
          dim(goSimMatrix_BP_new)[1],",",
          dim(goSimMatrix_MF_new)[1]," GOs ploted in CC, BP and MF, respectively")
  save(goSimMatrix_CC,goSimMatrix_BP,goSimMatrix_MF,
       file = paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".RData",sep = ""))
  message("Nice! plot exported and RData saved!")
}

# now plot
ReduceDim_GO_Plot(AR_CNTRL_enrich,Dataset_Name = "AR_CNTRL_enrich")
ReduceDim_GO_Plot(PRF_CNTRL_enrich,Dataset_Name = "PRF_CNTRL_enrich")
ReduceDim_GO_Plot(FPM_CNTRL_enrich,Dataset_Name = "FPM_CNTRL_enrich")
ReduceDim_GO_Plot(SMP_CNTRL_enrichh,Dataset_Name = "SMP_CNTRL_enrich")
ReduceDim_GO_Plot(SMP_FMP_enrich,Dataset_Name = "SMP_FMP_enrich")
