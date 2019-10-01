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
#Enrich_Results_thres005 = Go_Enrich_Plot(total_genes_all,
#                                         sig_genes_all,
#                                         TestingSubsetNames,
#                                         GOthres = 0.05,
#                                         keyword = "GO_Enrichment_qval01_pval005_test0924")

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

#######################################################################################
#                             6. Enrichment of Interpro                               #
#                            Join datatest and Compilation                            #
#######################################################################################
# main body of the function 
#Interpro_Enrich_Results_thres005 = 
#      InterPro_Enrich(total_genes_all,
#                      sig_genes_all,
#                      TestingSubsetNames,
#                      IPthres = 0.05,
#                      biomart="ensembl",
#                      dataset="btaurus_gene_ensembl",
#                      Identifier = "external_gene_name",
#                      attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
#                      keyword = "Interpro_Enrichment_thres005_0925")

## load in the results just created. Tip: use the key word above
load("Interpro_Enrichment_thres005_0925.RData")

########
library(tidyverse)
match_family = read.table("entry.list",sep = "\t",header = T)
match_family = as_tibble(match_family) %>% dplyr::select(ENTRY_AC,ENTRY_TYPE) %>% 
rename(InterproID = ENTRY_AC,Type =ENTRY_TYPE)

# parse 1 by 1, and attach the space name in the end
compile_select_index = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","hitsPerc")
### group 1
AR_CNTRL_enrich_IP = Parse_Interpro_Results(Interpro_results_b[2]) 
names(AR_CNTRL_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
AR_CNTRL_enrich_IP = dplyr::select(AR_CNTRL_enrich_IP,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
AR_CNTRL_enrich_IP = merge(AR_CNTRL_enrich_IP,match_family,by = "InterproID")

#
PRF_CNTRL_enrich_IP = Parse_Interpro_Results(Interpro_results_b[3])
names(PRF_CNTRL_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
PRF_CNTRL_enrich_IP  = dplyr::select(PRF_CNTRL_enrich_IP,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
PRF_CNTRL_enrich_IP = merge(PRF_CNTRL_enrich_IP,match_family,by = "InterproID")

### group 2
FPM_CNTRL_enrich_IP = Parse_Interpro_Results(Interpro_results_b[1])
names(FPM_CNTRL_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
FPM_CNTRL_enrich_IP = dplyr::select(FPM_CNTRL_enrich_IP,compile_select_index) # %>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
FPM_CNTRL_enrich_IP = merge(FPM_CNTRL_enrich_IP,match_family,by = "InterproID")
#
SMP_CNTRL_enrich_IP= Parse_Interpro_Results(Interpro_results_b[4])
names(SMP_CNTRL_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
SMP_CNTRL_enrich_IP = dplyr::select(SMP_CNTRL_enrich_IP,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
SMP_CNTRL_enrich_IP = merge(SMP_CNTRL_enrich_IP,match_family,by = "InterproID")

#
SMP_FMP_enrich_IP = Parse_Interpro_Results(Interpro_results_b[5])
names(SMP_FMP_enrich_IP) = c("InterproID","Interpro_Name","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
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
#write.xlsx(Interpro_Enrich_Regression_005,file = "Interpro_Enrich_Regression_005_0925_withFamily2.xlsx")
Interpro_Enrich_Pregnancy_005 <- list("Full_join" = Interpro_Results_full_005_preg, "Inner_join" = Interpro_Results_inner_005_preg)
#write.xlsx(Interpro_Enrich_Pregnancy_005,file = "Interpro_Enrich_Pregnancy_005_0925_withFamily2.xlsx")



#######################################################################################
#                             7. Conversion to EntrezID                               #
#                                  Mesh Enrich                                        #
#######################################################################################
# Here we have the data compilation
# List with lenght 5 (one for each, easier for looping)
#    sig_genes_all
#    total_genes_all
#
# Now we convert
# two steps - 
# 1. using alias2Symbol(): get the symbol
# 2, convert using biomart
# Potentially loose some, but that's life
## transform ensemble ID or External Gene name to EntrezID
ensembl_try=useMart("ENSEMBL_MART_ENSEMBL",dataset="btaurus_gene_ensembl")
attributes_try = c("external_gene_name","entrezgene_accession","entrezgene_id") #"entrezgene_accession"
#,"reactome","reactome_gene")
# "entrezgene_trans_name"
gene_try = getBM(attributes=attributes_try,mart = ensembl_try)
#
match_source = dplyr::select(gene_try,external_gene_name,entrezgene_id)

for (i in c(1:5)){
  tmp1 = unlist(sig_genes_all[[i]])
  for (n in seq_along(tmp1)){
    trans = alias2Symbol(tmp1[n],species = "Bt", expand.symbols = F)
    if (length(trans) == 0){tmp1[n] = tmp1[n]}
    else {tmp1[n] = trans}
  }
  tmp2 = unlist(total_genes_all[[i]])
  for (m in seq_along(tmp2)){
    trans = alias2Symbol(tmp2[m],species = "Bt", expand.symbols = F)
    if (length(trans) == 0){tmp2[n] = tmp2[n]}
    else {tmp2[n] = trans}
  }
  #
  origin1 = data.frame(sig_genes_all[[i]])
  origin2 = data.frame(total_genes_all[[i]])
  #names(tmp1) = names(tmp2) = "external_gene_name"
  test_inner1 = data.frame(origin = origin1,name_trans = tmp1)
  test_inner2 = data.frame(origin = origin2,name_trans = tmp2)
  test_inner1 = dplyr::left_join(test_inner1,match_source,by=c("name_trans" = "external_gene_name"))
  test_inner2 = dplyr::left_join(test_inner2,match_source,by=c("name_trans" = "external_gene_name"))
  # substitue entrez id with ### following "LOC"
  index1 = is.na(test_inner1$entrezgene_id)
  test_inner1$entrezgene_id[index1] = sub("^LOC*","",test_inner1$gene)[index1]
  index2 = is.na(test_inner2$entrezgene_id)
  test_inner2$entrezgene_id[index2] = sub("^LOC*","",test_inner2$gene)[index2]
  #
  Sig_list_out[[i]] = test_inner1
  Total_list_out[[i]] = test_inner2
  names(Sig_list_out)[i] = names(sig_genes_all)[i]
  names(Total_list_out)[i] = names(total_genes_all)[i]
  t1 = table(is.na(test_inner1[,3]))
  message("this is for sig")
  print(t1)
  t2 = table(is.na(test_inner2[,3]))
  message("this is for total")
  print(t2)
}

# print out
require(openxlsx)
#write.xlsx(Sig_list_out,file = "test_convert_sig.xlsx")
#write.xlsx(Total_list_out,file = "test_convert_total.xlsx")

# Keep only the entrez ID: then we have one vector for each element of the list (some format as always)
Sig_list_out_entrez = list()
Total_list_out_entrez = list()
for (i in c(1:5)){
  Sig_list_out_entrez[[i]] = data.frame(Sig_list_out[[i]])$entrezgene_id
  names(Sig_list_out_entrez)[i] = names(Sig_list_out)[i]
  Total_list_out_entrez[[i]] = data.frame(Total_list_out[[i]])$entrezgene_id
  names(Total_list_out_entrez)[i] = names(Total_list_out)[i]
}
#str(Total_list_out_entrez)
#str(Sig_list_out_entrez)


#####################
##  run analysis   ##
####################
# just in case that does not work
keyword = "MESH_Enrichment_1001"
DB = paste(keyword,".RData",sep = "")
load(DB)
###
MESH_Enrich_Result1001 = MESH_Enrich(total_genes_all= Total_list_out_entrez,
                                     sig_genes_all = Sig_list_out_entrez,
                                     TestingSubsetNames = TestingSubsetNames,
                                     Meshthres = 0.05,
                                     dataset="MeSH.Bta.eg.db",
                                     keyword = "MESH_Enrichment_1001")
#######################################################################################
#                          7. Conversion to EntrezID　(use above)                     #
#                                Reactome Enrich                                      #
#######################################################################################
# Get data from web (to be specified)

# Read in database
# lowest_path
NCBI2Reactome_lowest_path = read.csv("NCBI2Reactome.txt",sep = "\t",header = F)
NCBI2Reactome_lowest_path_bt = dplyr::filter(NCBI2Reactome_lowest_path, V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_lowest_path_bt,10)
# all_path
NCBI2Reactome_all_path = read.csv("NCBI2Reactome_All_Levels.txt",sep = "\t",header = F)
NCBI2Reactome_all_path_bt = dplyr::filter(NCBI2Reactome_all_path,V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_all_path_bt)
# all_react
NCBI2Reactome_all_react = read.csv("NCBI2Reactome_PE_Reactions.txt",sep = "\t",header = F)
NCBI2Reactome_all_react_bt = dplyr::filter(NCBI2Reactome_all_react,V8 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V3,V4,V6,V7,V8) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2, 
                Reaction_Description = V3,
                ProteinID = V4,
                Protein_Description = V6,
                Source = V7, Species = V8)

# turn data input as charactor
NCBI2Reactome_all_react_bt[] <-   lapply(NCBI2Reactome_all_react_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_lowest_path_bt[] <- lapply(NCBI2Reactome_lowest_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_all_path_bt[] <-   lapply(NCBI2Reactome_all_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
# data massage done 

###
### just for testing
#Total_list_out_entrez_test = Total_list_out_entrez[1:2]
#Sig_list_out_entrez_test = Sig_list_out_entrez[1:2]
#TestingSubsetNames
#InputSource = NCBI2Reactome_all_react_bt
## testing ends

## all react
Reactome_Enrich_all_react_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                sig_genes_all=Sig_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames,
                                                InputSource=  NCBI2Reactome_all_react_bt,
                                                Reacthres = 0.05,
                                                keyword = "Reactome_Enrichment_all_react_1001")
## lowest path
Reactome_Enrich_lowest_path_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                 sig_genes_all=Sig_list_out_entrez,
                                                 TestingSubsetNames = TestingSubsetNames,
                                                 InputSource=  NCBI2Reactome_alowest_path_bt,
                                                 Reacthres = 0.05,
                                                 keyword = "Reactome_Enrich_lowest_path_1001")
## all path
Reactome_Enrich_all_path_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                   sig_genes_all=Sig_list_out_entrez,
                                                   TestingSubsetNames = TestingSubsetNames,
                                                   InputSource=  NCBI2Reactome_all_path_bt,
                                                   Reacthres = 0.05,
                                                   keyword = "Reactome_Enrich_all_path_1001")

#########################################################################################################################
