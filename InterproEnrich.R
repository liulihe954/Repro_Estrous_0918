
test= InterPro_Enrich(total_genes_all,
                      sig_genes_all,
                      TestingSubsetNames,
                      IPthres = 0.05,
                      biomart="ensembl",
                      dataset="btaurus_gene_ensembl",
                      Identifier = "external_gene_name",
                      attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                      keyword = "Interpro_Enrichment_test")



biomart= "ensembl"
dataset="btaurus_gene_ensembl"
host="http://www.ensembl.org"
attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description")
database = useMart(biomart)
genome = useDataset(dataset, mart = database)
gene = getBM(attributes,mart = genome)




