BiocManager::install("MeSH.db")
BiocManager::install("MeSH.Bta.eg.db")
library(MeSH.db)
library(MeSH.Bta.eg.db)

# raw data for retrive MESHid and all details linked
KEY = keys(MeSH.db, keytype = "MESHID")
List = select(MeSH.db, keys = KEY, columns = columns(MeSH.db), keytype = "MESHID")
Match_List = dplyr::select(List, MESHID, MESHTERM)
# head(Match_List) 

# Prepare Bta database
key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
list_Bta = select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>% 
  dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY == c("D","G")) %>% 
  dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))
# head(list_Bta,30)




