source("Function_Source.R")
source("Enrich_pipeline.R")
#######################################################################################
# input pre
str(Sig_list_out_entrez)
str(Total_list_out_entrez)
TestingSubsetNames
#######################################################################################
#                                   8.     kegg  Enrich                              #
#######################################################################################
KEGG_Enrichment_thres005_1008 = 
  Kegg_Enrich_Plot(Sig_list_out_entrez,
                   Total_list_out_entrez,
                   TestingSubsetNames, 
                   KEGGthres = 0.05,
                   species = "bta", 
                   id.type = "kegg",
                   keyword = "KEGG_Enrichment_thres005_1008")
#######################################################################################
load("KEGG_Enrichment_thres005_1008.RData")
# parse 1 by 1, and attach the space name in the end
compile_select_index = c("KEGGID","KEGGTERM","Total_Genes","Significant_Genes","pvalue","hitsPerc")
### group 1
AR_CNTRL_enrich_KEGG = Parse_Results(KEGG_results_b[2]) 
names(AR_CNTRL_enrich_KEGG) = c("KEGGID","KEGGTERM","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
AR_CNTRL_enrich_KEGG = dplyr::select(AR_CNTRL_enrich_KEGG,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))

#
PRF_CNTRL_enrich_KEGG = Parse_Results(KEGG_results_b[3])
names(PRF_CNTRL_enrich_KEGG) = c("KEGGID","KEGGTERM","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
PRF_CNTRL_enrich_KEGG  = dplyr::select(PRF_CNTRL_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))

### group 2
FPM_CNTRL_enrich_KEGG = Parse_Results(KEGG_results_b[1])
names(FPM_CNTRL_enrich_KEGG) = c("KEGGID","KEGGTERM","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
FPM_CNTRL_enrich_KEGG = dplyr::select(FPM_CNTRL_enrich_KEGG,compile_select_index) # %>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
#
SMP_CNTRL_enrich_KEGG= Parse_Results(KEGG_results_b[4])
names(SMP_CNTRL_enrich_KEGG) = c("KEGGID","KEGGTERM","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
SMP_CNTRL_enrich_KEGG = dplyr::select(SMP_CNTRL_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
#
SMP_FMP_enrich_KEGG = Parse_Results(KEGG_results_b[5])
names(SMP_FMP_enrich_KEGG) = c("KEGGID","KEGGTERM","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
SMP_FMP_enrich_KEGG =  dplyr::select(SMP_FMP_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))


#### group 1 - regression
KEGG_Results_full_005_reg <-
  dplyr::full_join(AR_CNTRL_enrich_KEGG,PRF_CNTRL_enrich_KEGG,
                   by = c("KEGGID" = "KEGGID")) 
#%>% tidyr::replace_na(list(InterproID.x = "Not Found",InterproID.y = "Not Found")) 
KEGG_Results_inner_005_reg <-  
  dplyr::inner_join(AR_CNTRL_enrich_KEGG,PRF_CNTRL_enrich_KEGG,
                    by = c("KEGGID" = "KEGGID")) 

#### group 2 - pregnancy
KEGG_Results_full_005_preg <-  
  dplyr::full_join(FPM_CNTRL_enrich_KEGG,SMP_CNTRL_enrich_KEGG,
                   by = c("KEGGID" = "KEGGID"))  %>% 
  dplyr::full_join(SMP_FMP_enrich_KEGG,
                   by = c("KEGGID" = "KEGGID"))  
# %>%  tidyr::replace_na(list(InterproID.x = "Not Found",InterproID.y = "Not Found",InterproID = "Not Found"))
KEGG_Results_inner_005_preg <-  
  dplyr::inner_join(FPM_CNTRL_enrich_KEGG,SMP_CNTRL_enrich_KEGG, 
                    by = c("KEGGID" = "KEGGID"))  %>% 
  dplyr::inner_join(SMP_FMP_enrich_KEGG,by = c("KEGGID" = "KEGGID")) 

require(openxlsx)
KEGG_Enrich_Regression_005 <- list("Full_join" = KEGG_Results_full_005_reg, "Inner_join" = KEGG_Results_inner_005_reg)
write.xlsx(KEGG_Enrich_Regression_005,file = "KEGG_Enrich_Regression_005_1007.xlsx")
KEGG_Enrich_Pregnancy_005 <- list("Full_join" = KEGG_Results_full_005_preg, "Inner_join" = KEGG_Results_inner_005_preg)
write.xlsx(KEGG_Enrich_Pregnancy_005,file = "KEGG_Enrich_Pregnancy_005_1007.xlsx")


