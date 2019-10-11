source("Function_Source.R")
source("Enrich_pipeline.R")
#######################################################
total_genes_all = Total_list_out_entrez
sig_genes_all = Sig_list_out_entrez
#TestingSubsetNames
MESH_Enrichment_1010 = MESH_Enrich(total_genes_all,
                   sig_genes_all,
                   TestingSubsetNames,
                   Meshthres = 0.05,
                   Sig_list_out = Sig_list_out,
                   MeshCate = c("G"),
                   dataset="MeSH.Bta.eg.db",
                   keyword = "MESH_Enrichment_1010_G")
                       
##############################################################################################################
load("MESH_Enrichment_1010_G.RData")
keyword1 = "Mesh_Enrich_Regression_005_1010_G.xlsx"
keyword2 = "Mesh_Enrich_Pregnancy_005_1010_G.xlsx"

# parse 1 by 1, and attach the space name in the end
compile_select_index = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","findG","hitsPerc")
### group 1
AR_CNTRL_enrich_KEGG = Parse_Results(Mesh_results_b[2]) 
names(AR_CNTRL_enrich_KEGG) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
AR_CNTRL_enrich_KEGG = dplyr::select(AR_CNTRL_enrich_KEGG,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
#
PRF_CNTRL_enrich_KEGG = Parse_Results(Mesh_results_b[3])
names(PRF_CNTRL_enrich_KEGG) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
PRF_CNTRL_enrich_KEGG  = dplyr::select(PRF_CNTRL_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
### group 2
FPM_CNTRL_enrich_KEGG = Parse_Results(Mesh_results_b[1])
names(FPM_CNTRL_enrich_KEGG) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
FPM_CNTRL_enrich_KEGG = dplyr::select(FPM_CNTRL_enrich_KEGG,compile_select_index) # %>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
#
SMP_CNTRL_enrich_KEGG= Parse_Results(Mesh_results_b[4])
names(SMP_CNTRL_enrich_KEGG) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
SMP_CNTRL_enrich_KEGG = dplyr::select(SMP_CNTRL_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
#
SMP_FMP_enrich_KEGG = Parse_Results(Mesh_results_b[5])
names(SMP_FMP_enrich_KEGG) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
SMP_FMP_enrich_KEGG =  dplyr::select(SMP_FMP_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))

#### group 1 - regression
KEGG_Results_full_005_reg <-
  dplyr::full_join(AR_CNTRL_enrich_KEGG,PRF_CNTRL_enrich_KEGG,
                   by = c("MeshID" = "MeshID")) 
#%>% tidyr::replace_na(list(InterproID.x = "Not Found",InterproID.y = "Not Found")) 
KEGG_Results_inner_005_reg <-  
  dplyr::inner_join(AR_CNTRL_enrich_KEGG,PRF_CNTRL_enrich_KEGG,
                    by = c("MeshID" = "MeshID")) 

#### group 2 - pregnancy
KEGG_Results_full_005_preg <-  
  dplyr::full_join(FPM_CNTRL_enrich_KEGG,SMP_CNTRL_enrich_KEGG,
                   by = c("MeshID" = "MeshID"))  %>% 
  dplyr::full_join(SMP_FMP_enrich_KEGG,
                   by = c("MeshID" = "MeshID"))  
# %>%  tidyr::replace_na(list(InterproID.x = "Not Found",InterproID.y = "Not Found",InterproID = "Not Found"))
KEGG_Results_inner_005_preg <-  
  dplyr::inner_join(FPM_CNTRL_enrich_KEGG,SMP_CNTRL_enrich_KEGG, 
                    by =c("MeshID" = "MeshID"))  %>% 
  dplyr::inner_join(SMP_FMP_enrich_KEGG,by = c("MeshID" = "MeshID"))

require(openxlsx)
KEGG_Enrich_Regression_005 <- list("Full_join" = KEGG_Results_full_005_reg, "Inner_join" = KEGG_Results_inner_005_reg)
write.xlsx(KEGG_Enrich_Regression_005,file = keyword1)
KEGG_Enrich_Pregnancy_005 <- list("Full_join" = KEGG_Results_full_005_preg, "Inner_join" = KEGG_Results_inner_005_preg)
write.xlsx(KEGG_Enrich_Pregnancy_005,file = keyword2)


