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

#
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
#str(Sig_list_out)
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


#head(NCBI2Reactome_all_path_bt,30)
#test = dplyr::filter(NCBI2Reactome_all_path_bt,ReactomeID == "R-BTA-4090294")
#dim(test)
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
Reactome_Enrich = function(total_genes_all,
                           sig_genes_all,
                           TestingSubsetNames,
                           InputSource,
                           Reacthres = 0.05,
                           #biomart="ensembl",
                           #dataset="btaurus_gene_ensembl",
                           #host="http://www.ensembl.org",
                           #attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                           keyword = "Reactome_Enrichment"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  Reactome_results_b = list()
  Reactome_results_b_raw = list()
  library(ggplot2);library(biomaRt);library(gage);library(magrittr);library(tidyverse)# load pkg\
  InputSource = NCBI2Reactome_all_path_bt # delete
  Reactome_gene =   unique(InputSource[,1])
  ReactomeID =      unique(InputSource[,2])
  ReactomeName =    unique(InputSource[,3])
  #ReactomeID = ReactomeID[1:300]
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Reactome to check: ",length(ReactomeID)," with total number of names: ",length(ReactomeName))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    #i = 1
    
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    #length(sig.genes)
    #length(total.genes)
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% Reactome_gene])
    S = length(sig.genes[sig.genes %in% Reactome_gene]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(ReactomeID=character(),
                     ReactomeTerm=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(ReactomeID)){
      if (j%%100 == 0) {message("tryingd on Reactome ",j," - ",ReactomeID[j]," - ",ReactomeName[j])}
      
      
      NCBI2Reactome_all_path_bt = 
        dplyr::filter(NCBI2Reactome_all_path,V6 == "Bos taurus") %>% 
        dplyr::select(V1,V2,V4,V5,V6) %>% 
        dplyr::rename(EntrezID = V1,
                      ReactomeID = V2,
                      Reactome_Description = V4, 
                      Source = V5, 
                      Species = V6)
      
      
      
    
      #head(NCBI2Reactome_all_path_bt,30)
      InputSource =  NCBI2Reactome_all_path_bt # prepare data
      ReactomeID = unique(InputSource[,2]) # extract the index
      str(InputSource)
      str(ReactomeID) # just take one column of the data.frame above
      j = 24 # pick up a random target
      ReactomeID[j] == "R-BTA-190236" # they are same 
      # now we check if filter can find the corresponding rows - test1
      test1 = dplyr::filter(InputSource, ReactomeID == ReactomeID[j])
      table(test1$ReactomeID == ReactomeID[j])
      # now we check if filter can find the corresponding rows - test2
      target = ReactomeID[j] # re-assign this value to something in the middle
      test2 = dplyr::filter(InputSource, ReactomeID == target)
      table(test2$ReactomeID == ReactomeID[j])
      # I wonder what is going on here???
      # Thank you!
      
      
      
      test1
      test2 = dplyr::filter(InputSource, ReactomeID == "R-BTA-190236")
      #test2
      dim(test1)[1] == dim(test2)[1]
      
      
      test = InputSource[which(InputSource$ReactomeID == ReactomeID[j]),]$EntrezID
      test
      dim(test)
   
      
      dim(subset(InputSource, ReactomeID == ReactomeID[j]))
      gENEs = unique(subset(InputSource, ReactomeID == ReactomeID[j])$EntrezID)
      gENEs
      
      m = length(total.genes[total.genes %in% gENEs]) 
      
      s = length(sig.genes[sig.genes %in% gENEs]) # 
      m
      s
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(ReactomeID = ReactomeID[j], 
                       ReactomeName = ReactomeName[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Reactome_results_b_raw[[i]] = final_raw; names(Reactome_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Reactomeid raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Reacthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Reactome_results_b[[i]] = final;names(Reactome_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched ReactomeID")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>% 
    #           top_n(dim(final)[1], wt= -pvalue)%>% 
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])), 
    #                              x="Hits (%)", y="GO term", 
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Reactome = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Reactome_results_b, Reactome_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant Reactome domains found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",Reacthres)
  message("Nice! - Reactome enrichment finished and data saved")}

#########################################################################################################################
#total_genes_all=Total_list_out_entrez
#sig_genes_all=Sig_list_out_entrez
#InputSource=  NCBI2Reactome_all_react_bt

#(dplyr::filter(InputSource,ProteinID == "R-BTA-72185"))
#table(duplicated(InputSource$Protein_Description))

##############################
### formating the results  ##
##############################
load("Reactome_Enrich_lowest_path_1008.RData")
keyword1 = "Reactome_Enrich_Regression_005_1008_lowest_path.xlsx"
keyword2 = "Reactome_Enrich_Pregnancy_005_1008_lowest_path.xlsx"

load("Reactome_Enrich_all_path_1008.RData")
keyword1 = "Reactome_Enrich_Regression_005_1008_all_path.xlsx"
keyword2 = "Reactome_Enrich_Pregnancy_005_1008_all_path.xlsx"

load("Reactome_Enrichment_all_react_1008.RData")
keyword1 = "Reactome_Enrich_Regression_005_1008_all_react.xlsx"
keyword2 = "Reactome_Enrich_Pregnancy_005_1008_all_react.xlsx"


# parse 1 by 1, and attach the space name in the end
compile_select_index = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","hitsPerc")
### group 1
AR_CNTRL_enrich_KEGG = Parse_Results(Reactome_results_b[2]) 
names(AR_CNTRL_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
AR_CNTRL_enrich_KEGG = dplyr::select(AR_CNTRL_enrich_KEGG,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
#
PRF_CNTRL_enrich_KEGG = Parse_Results(Reactome_results_b[3])
names(PRF_CNTRL_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
PRF_CNTRL_enrich_KEGG  = dplyr::select(PRF_CNTRL_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
### group 2
FPM_CNTRL_enrich_KEGG = Parse_Results(Reactome_results_b[1])
names(FPM_CNTRL_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
FPM_CNTRL_enrich_KEGG = dplyr::select(FPM_CNTRL_enrich_KEGG,compile_select_index) # %>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
#
SMP_CNTRL_enrich_KEGG= Parse_Results(Reactome_results_b[4])
names(SMP_CNTRL_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
SMP_CNTRL_enrich_KEGG = dplyr::select(SMP_CNTRL_enrich_KEGG,compile_select_index) #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
#
SMP_FMP_enrich_KEGG = Parse_Results(Reactome_results_b[5])
names(SMP_FMP_enrich_KEGG) = c("ReactomeID","ReactomeTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","hitsPerc")
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
write.xlsx(KEGG_Enrich_Regression_005,file = keyword1)
KEGG_Enrich_Pregnancy_005 <- list("Full_join" = KEGG_Results_full_005_preg, "Inner_join" = KEGG_Results_inner_005_preg)
write.xlsx(KEGG_Enrich_Pregnancy_005,file = keyword2)




