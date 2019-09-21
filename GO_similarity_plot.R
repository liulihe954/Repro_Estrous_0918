library(GOSemSim)
library(corrplot)

goList = # vector of significant GO terms
goSimMat = mgoSim(goList, goList, ont="BP", measure="Jiang", combine=NULL)
corrplot(goSimMat, tl.col = "black", tl.cex = 0.8, method = "shade", order = "hclust", hclust.method = "centroid", is.corr = FALSE)


#
##########################################################################################
#                         6. reduce hierarcy FOR Similarity                              #
#########################################################################################
# dataset pre
# BiocManager::install("org.Bt.eg.db")
library("org.Bt.eg.db")
semData_BP <- godata('org.Bt.eg.db', ont="BP", computeIC=T) #ont="BP"
semData_MF <- godata('org.Bt.eg.db', ont="MF", computeIC=T) #ont="BP"
semData_CC <- godata('org.Bt.eg.db', ont="CC", computeIC=T) #ont="BP"
class(semData_BP)
#AR_CNTRL
AR_CNTRL_enrich_BP_List = dplyr::filter(AR_CNTRL_enrich,namespace_1003 == "biological_process") %>% 
  dplyr::select(go_id) %>% unlist();attributes(AR_CNTRL_enrich_BP_List) = NULL
AR_CNTRL_enrich_CC_List = dplyr::filter(AR_CNTRL_enrich,namespace_1003 == "cellular_component") %>% 
  dplyr::select(go_id) %>% unlist();attributes(AR_CNTRL_enrich_CC_List) = NULL
AR_CNTRL_enrich_MF_List = dplyr::filter(AR_CNTRL_enrich,namespace_1003 == "molecular_function") %>% 
  dplyr::select(go_id) %>% unlist();attributes(AR_CNTRL_enrich_MF_List) = NULL

goSimMat_test = mgoSim(AR_CNTRL_enrich_MF_List,
                       AR_CNTRL_enrich_MF_List,
                       semData = semData_MF,
                       measure="Jiang", combine="BMA") # combind why use null?

corrplot(goSimMat_test, tl.col = "black", tl.cex = 0.8, 
         method = "shade", order = "hclust", 
         hclust.method = "centroid", is.corr = FALSE)
