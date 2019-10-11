source("Function_Source.R")
source("Enrich_pipeline.R")
#######################################################################################
#                             6. Enrichment of Interpro                               #
#                            Join datatest and Compilation                            #
#######################################################################################
# main body of the function 
Interpro_Enrich_Results_thres005 = 
  InterPro_Enrich(total_genes_all,
                  sig_genes_all,
                  TestingSubsetNames,
                  IPthres = 0.05,
                  biomart="ensembl",
                  dataset="btaurus_gene_ensembl",
                  Identifier = "external_gene_name",
                  attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                  keyword = "Interpro_Enrichment_thres005_1011")

## load in the results just created. Tip: use the key word above
load("Interpro_Enrichment_thres005_1011.RData")

########
library(tidyverse)
match_family = read.table("entry.list",sep = "\t",header = T)
match_family = as_tibble(match_family) %>% dplyr::select(ENTRY_AC,ENTRY_TYPE) %>% 
  rename(InterproID = ENTRY_AC,Type =ENTRY_TYPE)

# parse 1 by 1, and attach the space name in the end
compile_select_index = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","findG","hitsPerc")
### group 1
AR_CNTRL_enrich_IP = Parse_Interpro_Results(Interpro_results_b[2]) 
names(AR_CNTRL_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
AR_CNTRL_enrich_IP = dplyr::select(AR_CNTRL_enrich_IP,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
AR_CNTRL_enrich_IP = merge(AR_CNTRL_enrich_IP,match_family,by = "InterproID")

#
PRF_CNTRL_enrich_IP = Parse_Interpro_Results(Interpro_results_b[3])
names(PRF_CNTRL_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
PRF_CNTRL_enrich_IP  = dplyr::select(PRF_CNTRL_enrich_IP,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
PRF_CNTRL_enrich_IP = merge(PRF_CNTRL_enrich_IP,match_family,by = "InterproID")

### group 2
FPM_CNTRL_enrich_IP = Parse_Interpro_Results(Interpro_results_b[1])
names(FPM_CNTRL_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
FPM_CNTRL_enrich_IP = merge(FPM_CNTRL_enrich_IP,match_family,by = "InterproID")
#
SMP_CNTRL_enrich_IP= Parse_Interpro_Results(Interpro_results_b[4])
names(SMP_CNTRL_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
SMP_CNTRL_enrich_IP = dplyr::select(SMP_CNTRL_enrich_IP,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
SMP_CNTRL_enrich_IP = merge(SMP_CNTRL_enrich_IP,match_family,by = "InterproID")

#
SMP_FMP_enrich_IP = Parse_Interpro_Results(Interpro_results_b[5])
names(SMP_FMP_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
SMP_FMP_enrich_IP =  dplyr::select(SMP_FMP_enrich_IP,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
SMP_FMP_enrich_IP = merge(SMP_FMP_enrich_IP,match_family,by = "InterproID")

#### group 1 - regression
Interpro_Results_full_005_reg <-
  dplyr::full_join(AR_CNTRL_enrich_IP,PRF_CNTRL_enrich_IP,
                   by = c("InterproID" = "InterproID")) 
#%>% tidyr::replace_na(list(InterproID.x = "Not Found",InterproID.y = "Not Found")) 
Interpro_Results_inner_005_reg <-  
  dplyr::inner_join(AR_CNTRL_enrich_IP,PRF_CNTRL_enrich_IP,
                    by = c("InterproID" = "InterproID"))
#### group 2 - pregnancy
Interpro_Results_full_005_preg <-  
  dplyr::full_join(FPM_CNTRL_enrich_IP,SMP_CNTRL_enrich_IP,
                   by = c("InterproID" = "InterproID")) %>% 
  dplyr::full_join(SMP_FMP_enrich_IP,
                   by = c("InterproID" = "InterproID")) 
# %>%  tidyr::replace_na(list(InterproID.x = "Not Found",InterproID.y = "Not Found",InterproID = "Not Found"))
Interpro_Results_inner_005_preg <-  
  dplyr::inner_join(FPM_CNTRL_enrich_IP,SMP_CNTRL_enrich_IP, 
                    by = c("InterproID" = "InterproID")) %>% 
  dplyr::inner_join(SMP_FMP_enrich_IP,by = c("InterproID" = "InterproID"))


require(openxlsx)
Interpro_Enrich_Regression_005 <- list("Full_join" = Interpro_Results_full_005_reg, "Inner_join" = Interpro_Results_inner_005_reg)
write.xlsx(Interpro_Enrich_Regression_005,file = "Interpro_Enrich_Regression_005_1011_withFamily2.xlsx")
Interpro_Enrich_Pregnancy_005 <- list("Full_join" = Interpro_Results_full_005_preg, "Inner_join" = Interpro_Results_inner_005_preg)
write.xlsx(Interpro_Enrich_Pregnancy_005,file = "Interpro_Enrich_Pregnancy_005_1011_withFamily2.xlsx")
