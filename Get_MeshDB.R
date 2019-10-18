# Save mesh for later use
library(MeSH.db);library(MeSH.Bta.eg.db)
library(tidyverse);library(gage);library(magrittr)
#
#keyword_outer = "MeshDB"
#DB = paste(keyword_outer,".RData",sep = "")
#load(DB)
#
#library(ggplot2);library(biomaRt) # load pkg
# raw data for retrive MESHid and all details linked
key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
#length(key_Bta)
List = MeSHDbi::select(MeSH.db, keys = key_Bta, 
                       columns = c("MESHID","MESHTERM"), 
                       keytype = "MESHID")

library(MeSH.Bta.eg.db)
key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
list_Bta = MeSHDbi::select(MeSH.Bta.eg.db, 
                           keys = key_Bta, 
                           columns = columns(MeSH.Bta.eg.db)[1:3], 
                           keytype = "MESHID") %>% 
  dplyr::filter(MESHCATEGORY %in% c("D","G")) %>% 
  dplyr::left_join(List,by= c("MESHID" = "MESHID"))

#List = MeSHDbi::select(MeSH.db, keys = key_Bta[3], columns = columns(MeSH.db), keytype = "MESHID")

#str(List)


#dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>%

#%>% 
 # dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY %in% c("D","G")) %>% 
#  dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))
save(key_Bta,List,list_Bta,file = "MeshDB.RData")


