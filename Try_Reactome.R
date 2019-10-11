source("Function_Source.R")
source("Enrich_pipeline.R")
#######################################################################################
#                               Prepare Reactome Database                             #
#######################################################################################
# Read in database
# lowest_path
NCBI2Reactome_lowest_path = read.csv("NCBI2Reactome.txt",sep = "\t",header = F)
NCBI2Reactome_lowest_path_bt = dplyr::filter(NCBI2Reactome_lowest_path, V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)

#head(NCBI2Reactome_lowest_path_bt,10)
# all_path
NCBI2Reactome_all_path = read.csv("NCBI2Reactome_All_Levels.txt",sep = "\t",header = F)
NCBI2Reactome_all_path_bt = 
  dplyr::filter(NCBI2Reactome_all_path,V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
dplyr::rename(EntrezID = V1,
              ReactomeID = V2,
              Reactome_Description = V4, 
              Source = V5, 
              Species = V6)

#head(NCBI2Reactome_all_path_bt,30)
#rm(list = c("test1","test2","target","j","InputSource","ReactomeID"))
#str(NCBI2Reactome_all_react_bt)
#InputSource = NCBI2Reactome_all_react_bt
#ReactomeID =  unique(InputSource[,2])
#j = 24
#ReactomeID[j] == "R-BTA-72127"
#test1 = dplyr::filter(NCBI2Reactome_all_react_bt,ReactomeID == ReactomeID[j])
#table(test1$ReactomeID == ReactomeID[j])
#
#target = ReactomeID[j]
#test2 = dplyr::filter(NCBI2Reactome_all_react_bt,ReactomeID == target)
#table(test2$ReactomeID == ReactomeID[j])
#test

#head(NCBI2Reactome_all_path_bt)
# all_react
NCBI2Reactome_all_react = read.csv("NCBI2Reactome_PE_Reactions.txt",sep = "\t",header = F)
NCBI2Reactome_all_react_bt = 
  dplyr::filter(NCBI2Reactome_all_react,V8 == "Bos taurus") %>% 
  dplyr::select(V1,V4,V6,V2,V3,V7,V8) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V4, 
                Reaction_Description = V6,
                ProteinID = V2,
                Protein_Description = V3,
                Source = V7, Species = V8)


#head(NCBI2Reactome_all_react_bt,50)
# turn data input as charactor
NCBI2Reactome_all_react_bt[] <-   lapply(NCBI2Reactome_all_react_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_lowest_path_bt[] <- lapply(NCBI2Reactome_lowest_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_all_path_bt[] <-   lapply(NCBI2Reactome_all_path_bt, function(x) if(is.factor(x)) as.character(x) else x)


# alternative: get from the web; dont forget to massage them like above (selecting, filtering, renaming)
library(data.table)
Reactome_lowest_path <- fread('https://reactome.org/download/current/NCBI2Reactome.txt')
Reactome_all_path <- fread('https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt')
Reactome_all_react <- fread('https://reactome.org/download/current/NCBI2ReactomeReactions.txt')



#### run test
str(Sig_list_out_entrez)
str(Total_list_out_entrez)
TestingSubsetNames
## all react
Reactome_Enrich_all_react_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                 sig_genes_all=Sig_list_out_entrez,
                                                 TestingSubsetNames = TestingSubsetNames,
                                                 InputSource=  NCBI2Reactome_all_react_bt,
                                                 Sig_list_out = Sig_list_out,
                                                 Reacthres = 0.05,
                                                 keyword = "Reactome_Enrichment_all_react_1011")
## lowest path
Reactome_Enrich_lowest_path_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                   sig_genes_all=Sig_list_out_entrez,
                                                   TestingSubsetNames = TestingSubsetNames,
                                                   InputSource=  NCBI2Reactome_lowest_path_bt,
                                                   Sig_list_out = Sig_list_out,
                                                   Reacthres = 0.05,
                                                   keyword = "Reactome_Enrich_lowest_path_1011")



## all path
Reactome_Enrich_all_path_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                sig_genes_all=Sig_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames,
                                                InputSource=  NCBI2Reactome_all_path_bt,
                                                Sig_list_out = Sig_list_out,
                                                Reacthres = 0.05,
                                                keyword = "Reactome_Enrich_all_path_1011")
#########################################################################################################################
#########################################################################################################################
#total_genes_all=Total_list_out_entrez
#sig_genes_all=Sig_list_out_entrez
#InputSource=  NCBI2Reactome_all_react_bt

#(dplyr::filter(InputSource,ProteinID == "R-BTA-72185"))
#table(duplicated(InputSource$Protein_Description))
##############################
### formating the results  ##
##############################
All_dataset = c("Reactome_Enrich_lowest_path_1011.RData",
                "Reactome_Enrich_all_path_1011.RData",
                "Reactome_Enrichment_all_react_1011.RData")
Keyword1 = c("Reactome_Enrich_Regression_005_1011_lowest_path.xlsx",
             "Reactome_Enrich_Regression_005_1011_all_path.xlsx",
             "Reactome_Enrich_Regression_005_1011_all_react.xlsx")
Keyword2 = c("Reactome_Enrich_Pregnancy_005_1011_lowest_path.xlsx",
             "Reactome_Enrich_Pregnancy_005_1011_all_path.xlsx",
             "Reactome_Enrich_Pregnancy_005_1011_all_react.xlsx")
for ( i in seq_along(All_dataset)){
  compile_select_index = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","findG","hitsPerc")
  ### group 1
  AR_CNTRL_enrich_KEGG = Parse_Results(Reactome_results_b[2]) 
  names(AR_CNTRL_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  AR_CNTRL_enrich_KEGG = dplyr::select(AR_CNTRL_enrich_KEGG,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))

  PRF_CNTRL_enrich_KEGG = Parse_Results(Reactome_results_b[3])
  names(PRF_CNTRL_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  PRF_CNTRL_enrich_KEGG  = dplyr::select(PRF_CNTRL_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
  ### group 2
  FPM_CNTRL_enrich_KEGG = Parse_Results(Reactome_results_b[1])
  names(FPM_CNTRL_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  FPM_CNTRL_enrich_KEGG = dplyr::select(FPM_CNTRL_enrich_KEGG,compile_select_index) # %>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
  #
  SMP_CNTRL_enrich_KEGG= Parse_Results(Reactome_results_b[4])
  names(SMP_CNTRL_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  SMP_CNTRL_enrich_KEGG = dplyr::select(SMP_CNTRL_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
  #
  SMP_FMP_enrich_KEGG = Parse_Results(Reactome_results_b[5])
  names(SMP_FMP_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  SMP_FMP_enrich_KEGG =  dplyr::select(SMP_FMP_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
  #### group 1 - regression
  KEGG_Results_full_005_reg <-
    dplyr::full_join(AR_CNTRL_enrich_KEGG,PRF_CNTRL_enrich_KEGG,
                   by = c("ReactomeID" = "ReactomeID")) 
  #%>% tidyr::replace_na(list(InterproID.x = "Not Found",InterproID.y = "Not Found")) 
  KEGG_Results_inner_005_reg <-  
    dplyr::inner_join(AR_CNTRL_enrich_KEGG,PRF_CNTRL_enrich_KEGG,
                    by = c("ReactomeID" = "ReactomeID")) 
  #### group 2 - pregnancy
  KEGG_Results_full_005_preg <-  
    dplyr::full_join(FPM_CNTRL_enrich_KEGG,SMP_CNTRL_enrich_KEGG,
                     by = c("ReactomeID" = "ReactomeID"))  %>% 
    dplyr::full_join(SMP_FMP_enrich_KEGG,
                   by = c("ReactomeID" = "ReactomeID"))  
  # %>%  tidyr::replace_na(list(InterproID.x = "Not Found",InterproID.y = "Not Found",InterproID = "Not Found"))
  KEGG_Results_inner_005_preg <-  
    dplyr::inner_join(FPM_CNTRL_enrich_KEGG,SMP_CNTRL_enrich_KEGG, 
                    by = c("ReactomeID" = "ReactomeID"))  %>% 
    dplyr::inner_join(SMP_FMP_enrich_KEGG,by = c("ReactomeID" = "ReactomeID")) 
  require(openxlsx)
  KEGG_Enrich_Regression_005 <- list("Full_join" = KEGG_Results_full_005_reg, "Inner_join" = KEGG_Results_inner_005_reg)
  write.xlsx(KEGG_Enrich_Regression_005,file = Keyword1[i])
  KEGG_Enrich_Pregnancy_005 <- list("Full_join" = KEGG_Results_full_005_preg, "Inner_join" = KEGG_Results_inner_005_preg)
  write.xlsx(KEGG_Enrich_Pregnancy_005,file = Keyword2[i])
} 



  
