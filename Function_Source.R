##############################################################################################################
library(readxl);library(ggplot2);library(biomaRt);library(tidyverse)
library(biomaRt);library(GOSemSim);library(corrplot)
##############################################################################################################
Go_Enrich_Plot = function(total.genes,
                          sig.genes,
                          TestingSubsetNames,
                          GOthres = 0.05,
                          biomart="ensembl",
                          dataset="btaurus_gene_ensembl",
                          host="http://www.ensembl.org",
                          attributes = c("external_gene_name","go_id","name_1006"),
                          keyword = "GO_Enrichment_thres_point1_5sets"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  GO_results_b = list()
  GO_results_b_raw = list()
  library(ggplot2);library(biomaRt);library(gage);library(magrittr)# load pkg
  ## Analysis bosTau annotation: GO
  database = useMart(biomart)
  genome = useDataset(dataset, mart = database)
  gene = getBM(attributes,mart = genome)
  goName = unique(gene[,c(2,3)]); goName = goName[order(goName$go_id),];goName = goName[-1,]
  GO = goName$go_id
  Name = goName$name_1006
  genesGO = unique(subset(gene,go_id != "")$external_gene_name)[-1]
  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of GO sets to check: ",length(GO)," with total number of names: ",length(Name))
  # plot
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesGO])
    S = length(sig.genes[sig.genes %in% genesGO]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(GO=character(),
                     Name=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(GO)){
      if (j%%100 == 0) {message("tryingd on GO ",j," - ",GO[j]," - ",Name[j])}
      gENEs = subset(gene, go_id == GO[j])$external_gene_name # all gene in target GO
      m = length(total.genes[total.genes %in% gENEs]) # genes from target GO and in our dataset
      s = length(sig.genes[sig.genes %in% gENEs]) # # genes from target GO also in the non-preserved module
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(GO = GO[j], 
                       Name = Name[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("GOID","GO_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    GO_results_b_raw[[i]] = final_raw; names(GO_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched GO raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= GOthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("GOID","GO_Name", "Total_Genes", "Significant_Genes", "pvalue","ExternalLoss_total","InternalLoss_sig")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    GO_results_b[[i]] = final;names(GO_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched GO")
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
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_GO = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(GO_results_b, GO_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant GO terms found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",GOthres)
  message("Nice! - GO enrichment finished and data saved")}
##############################################################################################################
Parse__Results = function(KEGG_results_b){
  all_enrich_KEGG = data.frame(ID=character(),
                               Description=character(),
                               GeneRatio=character(),
                               BgRatio=character(),
                               pvalue=numeric(),
                               p.adjust=numeric(),
                               qvalue=numeric(),
                               geneID=character(),
                               Count=numeric(),
                               stringsAsFactors=FALSE)
  for (i in 1:length(KEGG_results_b)){
    len = dim(data.frame(KEGG_results_b[i]))[1]
    if (len> 0){
      all_enrich_KEGG = rbind(all_enrich_KEGG,data.frame(KEGG_results_b[i]))
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich_KEGG)[1]
  total_modules = length(KEGG_results_b)
  print(paste(total_hits,"hits found in",total_modules,"non-preserved modules"))
  return(ParseResults = all_enrich_KEGG)
}


Parse_GO_Results = function(GO_results_b){
  all_enrich_GO = data.frame(ID=character(),
                             Description=character(),
                             GeneRatio=character(),
                             BgRatio=character(),
                             pvalue=numeric(),
                             p.adjust=numeric(),
                             qvalue=numeric(),
                             geneID=character(),
                             Count=numeric(),
                             stringsAsFactors=FALSE)
  for (i in 1:length(GO_results_b)){
    len = dim(data.frame(GO_results_b[i]))[1]
    if (len> 0){
      all_enrich_GO = rbind(all_enrich_GO,data.frame(GO_results_b[i]))
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich_GO)[1]
  total_modules = length(GO_results_b)
  print(paste(total_hits,"hits found in",total_modules,"non-preserved modules"))
  return(ParseResults = all_enrich_GO)
}



