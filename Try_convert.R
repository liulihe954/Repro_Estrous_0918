library(IRanges)
library(AnnotationDbi)

library(org.Bt.eg.db)
library(limma)

?alias2Symbol()
gene_try[gene_try$external_gene_name %in% test_assesion3,]

test_assesion = c("LOC100139916","LOC100297594","LOC100298356","LOC100336786",
                  "LOC100337076","LOC100847245","LOC100848491","LOC100848655")

test_assesion2 = c("ANKRD57","CAGE1","DMBT1","GALR3","ULBP27","ZNF791")
test_assesion3 = c("SUGCT")

sessionInfo()

alias2Symbol(alias, species = "Hs", expand.symbols = FALSE)
alias2SymbolTable(alias, species = "Hs")
alias2SymbolUsingNCBI(alias, gene.info.file,
                      required.columns = c("GeneID","Symbol","description"))

find.package("IRanges")
find.package("AnnotationDbi")
find.package("limma")
library(limma)


alias =c("LOC100297594")

alias2SymbolTable(alias, species = "Bt",expand.symbols = T)



library(IRanges)

detach(package:limma)
browseVignettes("IRanges")
if (!require("BiocManager"))
install.packages("BiocManager")
BiocManager::install("IRanges")



library(meshr)
?meshr::annotation()


data(sig.geneid.cummeRbund)
names(sig.geneid.cummeRbund)

data(geneid.cummeRbund)


dim(sig.geneid.cummeRbund)


