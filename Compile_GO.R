source("Function_Source.R")
#getwd()
#rm(list = ls())

library(GOSemSim)
library(corrplot)


load("GO_Enrichment_qval01_pval005.RData")



# group 1
AR_CNTRL_enrich = Parse_GO_Results(GO_results_b_raw[])
PRF_CNTRL_enrich = Parse_GO_Results(GO_results_b_raw[])
# group 2
FPM_CNTRL_enrich = Parse_GO_Results(GO_results_b_raw[])
SMP_CNTRL_enrich= Parse_GO_Results(GO_results_b_raw[])
SMP_FMP_enrich = Parse_GO_Results(GO_results_b[])
#
# thres at 0.05
AR_CNTRL_enrich <- AR_CNTRL_enrich %>% 
  

KEGG_Results_Day14_Z5 <-  
  KEGG_Results_Day14 %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::select(KEGG.ID,KEGG_Name,pvalue) %>% 
  dplyr::group_by(KEGG.ID) %>% 
  dplyr::arrange((KEGG.ID))%>% 
  dplyr::rename(KEGG_Name_14 = KEGG_Name)
#KEGG_Results_Day14_Z1

KEGG_Results_Day42 = Parse_KEGG_Results(KEGG_results_b_42)
KEGG_Results_Day42_Z5 <-  
  KEGG_Results_Day42 %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::select(KEGG.ID,KEGG_Name,pvalue) %>% 
  dplyr::group_by(KEGG.ID) %>% 
  dplyr::arrange((KEGG.ID)) %>% 
  dplyr::rename(KEGG_Name_42 = KEGG_Name)

KEGG_Results_Day84 = Parse_KEGG_Results(KEGG_results_b_84)
KEGG_Results_Day84_Z5 <-  
  KEGG_Results_Day84 %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::select(KEGG.ID,KEGG_Name,pvalue) %>% 
  dplyr::group_by(KEGG.ID) %>% 
  dplyr::arrange((KEGG.ID))%>% 
  dplyr::rename(KEGG_Name_84 = KEGG_Name)
#KEGG_Results_Day84_Z1

####
KEGG_Results_Top50_full_Z5 <-  
  dplyr::full_join(KEGG_Results_Day14_Z5,KEGG_Results_Day42_Z5,by = c("KEGG.ID" = "KEGG.ID")) %>% 
  dplyr::full_join(KEGG_Results_Day84_Z5,by = c("KEGG.ID" = "KEGG.ID"))%>% 
  tidyr::replace_na(list(KEGG_Name.x = "Not Found",KEGG_Name.y = "Not Found")) 

KEGG_Results_Top50_inner_Z5 <- 
  dplyr::inner_join(KEGG_Results_Day14_Z5,KEGG_Results_Day42_Z5,by = c("KEGG.ID" = "KEGG.ID")) %>% 
  dplyr::inner_join(KEGG_Results_Day84_Z5,by = c("KEGG.ID" = "KEGG.ID")) 


require(openxlsx)
kegg_out_list_Z5 <- list("Full_join" = KEGG_Results_Top50_full_Z5, "Inner_join" = KEGG_Results_Top50_inner_Z5)
setwd("/Users/liulihe95/Desktop/CoolHeat_Results_Top20")
write.xlsx(kegg_out_list_Z1, file = "KEGG_Results_Top20_Z5.xlsx")







load("GO_Enrichment_qval01_pval001.RData")











# plotting 

# BiocManager::install("org.Bt.eg.db")
library("org.Bt.eg.db")
parse_test = Parse_GO_Results(GO_results_b)
dim(parse_test)
goList = parse_test$AR_CNTRL.with.134.enriched.GO.GOID # vector of significant GO terms


semData <- godata('org.Bt.eg.db', ont="BP", computeIC=T) #ont="BP"
goSimMat = mgoSim(goList,goList,semData,measure="Jiang", combine=NULL) #
corrplot(goSimMat, tl.col = "black", tl.cex = 0.8, 
         method = "shade", order = "hclust", 
         hclust.method = "centroid", is.corr = FALSE)

