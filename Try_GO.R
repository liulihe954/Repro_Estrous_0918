source("Function_Source.R")
source("Enrich_pipeline.R")
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
#                                  4. Enrichment of GO                                #
#                               Join datatest and  Compilation                        #
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

########

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
#write.xlsx(GO_Enrich_Regression_005,file = "GO_Enrich_Regression_005_0921.xlsx")
GO_Enrich_Pregnancy_005 <- list("Full_join" = GO_Results_full_005_preg, "Inner_join" = GO_Results_inner_005_preg)
#write.xlsx(GO_Enrich_Pregnancy_005,file = "GO_Enrich_Pregnancy_005_0921.xlsx")

########################################################################################
#                            5.Take 0.01 - be more strengent                           #
########################################################################################
# This is for test pval, to see how many you'll get at which pval; change groups/thres first
Enrich_Out = FPM_CNTRL_enrich# choices: "FPM_CNTRL" "AR_CNTRL"  "PRF_CNTRL" "SMP_CNTRL" "SMP_FMP" 
GOthres_testp = 0.005
BP_List = dplyr::filter(Enrich_Out,pvalue<= GOthres_testp & namespace_1003 == "biological_process") %>% 
  dplyr::select(go_id) %>% unlist();attributes(BP_List) = NULL 
CC_List = dplyr::filter(Enrich_Out,pvalue<=GOthres_testp & namespace_1003 == "cellular_component") %>% 
  dplyr::select(go_id) %>% unlist();attributes(CC_List) = NULL
MF_List = dplyr::filter(Enrich_Out,pvalue<=GOthres_testp & namespace_1003 == "molecular_function") %>% 
  dplyr::select(go_id) %>% unlist();attributes(MF_List) = NULL
message(" BP: ",length(BP_List)," CC: ",length(CC_List)," MF: ",length(MF_List))

# Function takes Enrichment table as imput (same format as "AR_CNTRL_enrich")
#####
# ReduceDim_GO_Plot(AR_CNTRL_enrich,GOthres = 0.001, Dataset_Name = "AR_CNTRL_enrich")
# ReduceDim_GO_Plot(PRF_CNTRL_enrich,GOthres = 0.001,Dataset_Name = "PRF_CNTRL_enrich")
# ReduceDim_GO_Plot(FPM_CNTRL_enrich,GOthres = 0.01,Dataset_Name = "FPM_CNTRL_enrich")
# ReduceDim_GO_Plot(SMP_CNTRL_enrich,GOthres = 0.001,Dataset_Name = "SMP_CNTRL_enrich")
# ReduceDim_GO_Plot(SMP_FMP_enrich,GOthres = 0.001,Dataset_Name = "SMP_FMP_enrich")
