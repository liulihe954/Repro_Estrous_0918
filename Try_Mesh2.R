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
#
Sig_list_out = list()
Total_list_out = list()


#####################
##  run analysis   ##
####################
# just in case that does not work
keyword = "MeshDB"
DB = paste(keyword,".RData",sep = "")
load(DB)


######################################
################
## MeSH Analysis
library(org.Bt.eg.db)
library(meshr)
library(MeSH.db)
library(MeSH.Bta.eg.db)

# keys return the keys for the database contained in the MeSHdb object
key.symbol = keys(org.Bt.eg.db,  keytype = c("SYMBOL"))
entrezUniverse = select(org.Bt.eg.db, as.character(key.symbol), 
                        columns = c("ENTREZID"),keytype = "SYMBOL") %>% 
  dplyr::distinct(SYMBOL,.keep_all= TRUE)
# dim(entrezUniverse)
#####
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


# testing 
Sig_list_out_entrez_test2 = unlist(Sig_list_out_entrez[[4]]);attributes(Sig_list_out_entrez_test2) = NULL
Total_list_out_entrez_test2 = unlist(Total_list_out_entrez[[4]]);attributes(Total_list_out_entrez_test2) = NULL

Sig_list_out_entrez_test2 = Sig_list_out_entrez[4]#;attributes(Sig_list_out_entrez_test2) = NULL
Total_list_out_entrez_test2 = Total_list_out_entrez[4]#;attributes(Total_list_out_entrez_test2) = NULL
TestingSubsetNames = "test"

test = MESH_Enrich(Total_list_out_entrez_test2,
                   Sig_list_out_entrez_test2,
                   TestingSubsetNames,
                   Meshthres = 0.05,
                   MeshCate = c("D","G"),
                   #biomart="ensembl",
                   dataset="MeSH.Bta.eg.db",
                   #dataset= "btaurus_gene_ensembl",
                   #Identifier = "external_gene_name",
                   #attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                   keyword = "MESH_Enrichment_test1007")

load("MESH_Enrichment_test1007.RData")

##

## FULL GENES
genes.back = data.frame(total.genes)
colnames(genes.back) <- "ENSEMBL"
geneID.back <- merge(genes.back, entrezUniverse3, by ="ENSEMBL")
geneID2.back <- geneID.back[ !duplicated(geneID.back[,2]),]

## SIGNIFICANT GENES
genes.sig = data.frame(sig.genes)
colnames(genes.sig) <- "ENSEMBL"
geneID.sig <- merge(genes.sig, entrezUniverse3, by ="ENSEMBL")
geneID2.sig <- geneID.sig[ !duplicated(geneID.sig[,2]),]

## Total Genes
ns = length(geneID2.sig[,1])
nt = length(geneID2.back[,1])
cat(paste("Significant Genes:", ns, " and Backgroung Genes:", nt - ns), "\n")

### MeSH Phenomena and Processes
# Total_list_out_entrez_test = Total_list_out_entrez[2]
# Sig_list_out_entrez_test = Sig_list_out_entrez[2]
Sig_list_out_entrez_test2 = unlist(Sig_list_out_entrez[[4]]);attributes(Sig_list_out_entrez_test2) = NULL
Total_list_out_entrez_test2 = unlist(Total_list_out_entrez[[4]]);attributes(Total_list_out_entrez_test2) = NULL


head(Sig_list_out_entrez_test2);length(Sig_list_out_entrez_test2)
head(Total_list_out_entrez_test2);length(Total_list_out_entrez_test2)

meshParams <- new("MeSHHyperGParams", 
                  geneIds = Sig_list_out_entrez_test2, 
                  universeGeneIds = Total_list_out_entrez_test2, 
                  annotation = "MeSH.Bta.eg.db",
                  category = c("G"), database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")

meshR <- meshHyperGTest(meshParams)

out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH_Term_ID", "MeSH_Term_Name", 
                  "Total_Genes", "DE_Genes", "P-value")
print(unique(out), row.names = F)

#######
meshParams <- new("MeSHHyperGParams", 
                  geneIds = Sig_list_out_entrez_test2, 
                  universeGeneIds = Total_list_out_entrez_test2, 
                  annotation = "MeSH.Bta.eg.db",
                  category = c("G"), database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")

str(meshParams)
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH_Term_ID", "MeSH_Term_Name", 
                  "Total_Genes", "DE_Genes", "P-value")
print(unique(out), row.names = F)

#write.table(unique(out),file="MESH_Phenomena_process.txt",append = F, quote = F,sep ="\t",na = "NA",row.names=F,col.names=T)
