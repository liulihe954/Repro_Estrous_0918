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

#######################################################################################
#                             7. Conversion to EntrezID                               #
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
# Potentially lose some, but that's life
## MeSH Analysis
library(org.Bt.eg.db)
library(meshr)
library(MeSH.db)
library(MeSH.Bta.eg.db)
# keys return the keys for the database contained in the MeSHdb object
key.symbol = AnnotationDbi::keys(org.Bt.eg.db,  keytype = c("SYMBOL"))
entrezUniverse = AnnotationDbi::select(org.Bt.eg.db, as.character(key.symbol), 
                        columns = c("ENTREZID"),keytype = "SYMBOL") %>% 
  dplyr::distinct(SYMBOL,.keep_all= TRUE)
# dim(entrezUniverse)
#
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
  Sig_list_out[[i]] = tmp1;names(Sig_list_out)[i] = names(sig_genes_all)[i]
  
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
  Total_list_out[[i]] = tmp2;names(Total_list_out)[i] = names(total_genes_all)[i]
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


#######################################################################################
#                          7. Conversion to EntrezID　(use above)                     #
#                                Reactome Enrich                                      #
#######################################################################################
# Get data from web (to be specified)

#rm(NCBI2Reactome_all_react)
#rm(NCBI2Reactome_all_react_bt)

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
#head(NCBI2Reactome_all_path_bt)
# all_react
NCBI2Reactome_all_react = read.csv("NCBI2Reactome_PE_Reactions.txt",sep = "\t",header = F)
NCBI2Reactome_all_react_bt = 
  dplyr::filter(NCBI2Reactome_all_react,V8 == "Bos taurus") %>% 
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
#TestingSubsetNames_test
#InputSource = NCBI2Reactome_all_react_bt
### testing ends - results look good
#Reactome_Enrich_all_react_1001 = Reactome_Enrich(Total_list_out_entrez_test,
#                                                 Sig_list_out_entrez_test,
#                                                 TestingSubsetNames_test,
#                                                 InputSource,
#                                                 Reacthres = 0.05,
#                                                keyword = "Reactome_Enrichment_localtest")



# str(Total_list_out_entrez)
# str(Sig_list_out_entrez)
TestingSubsetNames

## all react
Reactome_Enrich_all_react_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                 sig_genes_all=Sig_list_out_entrez,
                                                 TestingSubsetNames = TestingSubsetNames,
                                                 InputSource=  NCBI2Reactome_all_react_bt,
                                                 Reacthres = 0.05,
                                                 keyword = "Reactome_Enrichment_all_react_1008")
## lowest path
Reactome_Enrich_lowest_path_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                   sig_genes_all=Sig_list_out_entrez,
                                                   TestingSubsetNames = TestingSubsetNames,
                                                   InputSource=  NCBI2Reactome_lowest_path_bt,
                                                   Reacthres = 0.05,
                                                   keyword = "Reactome_Enrich_lowest_path_1008")
## all path
Reactome_Enrich_all_path_1001 = Reactome_Enrich(total_genes_all=Total_list_out_entrez,
                                                sig_genes_all=Sig_list_out_entrez,
                                                TestingSubsetNames = TestingSubsetNames,
                                                InputSource=  NCBI2Reactome_all_path_bt,
                                                Reacthres = 0.05,
                                                keyword = "Reactome_Enrich_all_path_1008")

#########################################################################################################################
