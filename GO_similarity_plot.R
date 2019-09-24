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

goSim("GO:0055114","GO:0120162",semData=semData_BP,measure="Jiang")

goSimMat_test = mgoSim(AR_CNTRL_enrich_BP_List,
                       AR_CNTRL_enrich_BP_List,
                       semData=semData_BP,
                       measure="Jiang", combine=NULL) # combind why use null?
view(goSimMat_test)
goSimMat_test2 = mgoSim(list,
                       list,
                       semData=semData_BP,
                       measure="Jiang", combine=NULL) # combind why use null?
dev.off()
x = corrplot(goSimMat_test2, tl.col = "black", tl.cex = 0.8, 
         method = "shade", order = "hclust", 
         hclust.method = "centroid", is.corr = FALSE)



 