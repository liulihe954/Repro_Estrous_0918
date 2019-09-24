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
library("org.Bt.eg.db")
semData_BP <- godata('org.Bt.eg.db', ont="BP", computeIC=T) #ont="BP"
semData_MF <- godata('org.Bt.eg.db', ont="MF", computeIC=T) #ont="BP"
semData_CC <- godata('org.Bt.eg.db', ont="CC", computeIC=T) #ont="BP"


#AR_CNTRL - 0.01
AR_CNTRL_enrich_BP_List = dplyr::filter(AR_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "biological_process") %>% 
  dplyr::select(go_id) %>% unlist();attributes(AR_CNTRL_enrich_BP_List) = NULL
AR_CNTRL_enrich_CC_List = dplyr::filter(AR_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "cellular_component") %>% 
  dplyr::select(go_id) %>% unlist();attributes(AR_CNTRL_enrich_CC_List) = NULL
AR_CNTRL_enrich_MF_List = dplyr::filter(AR_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "molecular_function") %>% 
  dplyr::select(go_id) %>% unlist();attributes(AR_CNTRL_enrich_MF_List) = NULL

# PRF_CNTRL - 0.01
PRF_CNTRL_enrich_BP_List = dplyr::filter(PRF_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "biological_process") %>% 
  dplyr::select(go_id) %>% unlist();attributes(PRF_CNTRL_enrich_BP_List) = NULL
PRF_CNTRL_enrich_CC_List = dplyr::filter(PRF_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "cellular_component") %>% 
  dplyr::select(go_id) %>% unlist();attributes(PRF_CNTRL_enrich_CC_List) = NULL
PRF_CNTRL_enrich_MF_List = dplyr::filter(PRF_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "molecular_function") %>% 
  dplyr::select(go_id) %>% unlist();attributes(PRF_CNTRL_enrich_MF_List) = NULL

# - 0.05 / 0.01
FPM_CNTRL_enrich_BP_List = dplyr::filter(FPM_CNTRL_enrich,pvalue<=0.05 & namespace_1003 == "biological_process") %>% 
  dplyr::select(go_id) %>% unlist();attributes(FPM_CNTRL_enrich_BP_List) = NULL
FPM_CNTRL_enrich_CC_List = dplyr::filter(FPM_CNTRL_enrich,pvalue<=0.05 & namespace_1003 == "cellular_component") %>% 
  dplyr::select(go_id) %>% unlist();attributes(FPM_CNTRL_enrich_CC_List) = NULL
FPM_CNTRL_enrich_MF_List = dplyr::filter(FPM_CNTRL_enrich,pvalue<=0.05 & namespace_1003 == "molecular_function") %>% 
  dplyr::select(go_id) %>% unlist();attributes(FPM_CNTRL_enrich_MF_List) = NULL

#
SMP_CNTRL_enrich_BP_List = dplyr::filter(SMP_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "biological_process") %>% 
  dplyr::select(go_id) %>% unlist();attributes(SMP_CNTRL_enrich_BP_List) = NULL
SMP_CNTRL_enrich_CC_List = dplyr::filter(SMP_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "cellular_component") %>% 
  dplyr::select(go_id) %>% unlist();attributes(SMP_CNTRL_enrich_CC_List) = NULL
SMP_CNTRL_enrich_MF_List = dplyr::filter(SMP_CNTRL_enrich,pvalue<=0.001 & namespace_1003 == "molecular_function") %>% 
  dplyr::select(go_id) %>% unlist();attributes(SMP_CNTRL_enrich_MF_List) = NULL

#
SMP_FMP_enrich_BP_List = dplyr::filter(SMP_FMP_enrich,pvalue<=0.001 & namespace_1003 == "biological_process") %>% 
  dplyr::select(go_id) %>% unlist();attributes(SMP_FMP_enrich_BP_List) = NULL
SMP_FMP_enrich_CC_List = dplyr::filter(SMP_FMP_enrich,pvalue<=0.001 & namespace_1003 == "cellular_component") %>% 
  dplyr::select(go_id) %>% unlist();attributes(SMP_FMP_enrich_CC_List) = NULL
SMP_FMP_enrich_MF_List = dplyr::filter(SMP_FMP_enrich,pvalue<=0.001 & namespace_1003 == "molecular_function") %>% 
  dplyr::select(go_id) %>% unlist();attributes(SMP_FMP_enrich_MF_List) = NULL




# Function takes
ReduceDim_Go_Plot = function(Enrich_Out,GOthres,Database = "org.Bt.eg.db",measure="Jiang",combine=NULL){
  library(GOSemSim);library(corrplot);library(Database);library(tidyverse)
  semData_BP <- godata(Database, ont="BP", computeIC=T)
  semData_MF <- godata(Database, ont="MF", computeIC=T)
  semData_CC <- godata(Database, ont="CC", computeIC=T)
  BP_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "biological_process") %>% 
    dplyr::select(go_id) %>% unlist();attributes(BP_List) = NULL
  CC_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "cellular_component") %>% 
    dplyr::select(go_id) %>% unlist();attributes(CC_Listt) = NULL
  MF_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "molecular_function") %>% 
    dplyr::select(go_id) %>% unlist();attributes(MF_List) = NULL
  
  goSimMatrix_BP = GOSemSim::mgoSim(BP_List,BP_List,semData=semData_BP,measure=measure,combine = combine)
  goSimMatrix_CC = GOSemSim::mgoSim(CC_List,CC_List,semData=semData_MF,measure=measure,combine = combine)
  goSimMatrix_MF = GOSemSim::mgoSim(BP_List,BP_List,semData=semData_CC,measure=measure,combine = combine)
  
  
  
}



paste(AR_CNTRL_enrich_BP_List,AR_CNTRL_enrich$GO_Name[(AR_CNTRL_enrich$go_id %in% AR_CNTRL_enrich_BP_List)])

library(GOSemSim);library(corrplot)
goSimMat_AR_CNTRL_BP = GOSemSim::mgoSim(AR_CNTRL_enrich_BP_List,
                              AR_CNTRL_enrich_BP_List,
                              semData=semData_BP,
                              measure="Jiang", combine=NULL) # combind why use null?
suspectID = rownames(goSimMat_AR_CNTRL_BP)[is.na(goSimMat_AR_CNTRL_BP[,1])]
suspectID
  if (length(suspectID) != 0){AR_CNTRL_enrich_BP_List_new = AR_CNTRL_enrich_BP_List[-which(AR_CNTRL_enrich_BP_List == suspectID)]
  } else {AR_CNTRL_enrich_BP_List_new = AR_CNTRL_enrich_BP_List}
goSimMat_test_new = mgoSim(AR_CNTRL_enrich_BP_List_new,
                           AR_CNTRL_enrich_BP_List_new,
                           semData=semData_BP,
                           measure="Jiang", combine=NULL) # combind why use null?
colnames(goSimMat_test_new) = paste(AR_CNTRL_enrich_BP_List,AR_CNTRL_enrich$GO_Name[(AR_CNTRL_enrich$go_id %in% AR_CNTRL_enrich_BP_List)])
rownames(goSimMat_test_new) = paste(AR_CNTRL_enrich$GO_Name[(AR_CNTRL_enrich$go_id %in% AR_CNTRL_enrich_BP_List)],AR_CNTRL_enrich_BP_List)
x = corrplot(goSimMat_test_new,tl.col = "black", tl.cex = 0.4, 
             method = "shade", order = "hclust", 
             hclust.method = "centroid", is.corr = FALSE)

x



dev.off()
write(labels,"labels_AR_CONTROL_BP_0.01.csv")



colnames(M) <- c("a", "set", "of", "x", "labels", 1:6)
corrplot(M, method = "color")

?corrplot()


goSimMat_test[rownames(goSimMat_test) == "GO:0120162",]

1] "GO:0033089" "GO:0032534"
[3] "GO:0070120" "GO:0048843"
[5] "GO:1905606" "GO:0120162"





select(org.Bt.eg.db,keys ="GO:1905606",keytype = "GO","ENTREZID" )


view(goSimMat_test)

view(goSimMat_test)
goSimMat_test2 = mgoSim(list,
                        list,
                        semData=semData_BP,
                        measure="Jiang", combine=NULL) # combind why use null?
dev.off()
x = corrplot(goSimMat_test2, tl.col = "black", tl.cex = 0.8, 
             method = "shade", order = "hclust", 
             hclust.method = "centroid", is.corr = FALSE)






