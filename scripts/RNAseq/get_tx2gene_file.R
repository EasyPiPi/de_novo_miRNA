library(GenomicFeatures)
library(biomaRt)

################################################################################
# download data from ensembl
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
tmp <- datasets$description

ensembl <- useDataset("dmelanogaster_gene_ensembl",mart=ensembl)

# Dmel_TxDB <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
#                                  dataset="dmelanogaster_gene_ensembl")
# saveDb(Dmel_TxDB, "Dmel_ensembl_release_97")

attributes_ls <- listAttributes(ensembl)

# Drosophila ID
tx2gene <- getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id"), mart = ensembl)
write.csv(tx2gene, "dmel_tx2gene.csv")

################################################################################
# download from flybase
# orthologs
download.file("ftp://ftp.flybase.org/releases/FB2019_04/precomputed_files/orthologs/dmel_orthologs_in_drosophila_species_fb_2019_04.tsv.gz", "/home/yixin/NutstoreFiles/Nutstore/Wu_lab/Snakemake_projects/miR983_975/output/expression/dmel_orthologs_in_drosophila_species_fb_2019_04.tsv.gz")

# tx2gene
download.file("ftp://ftp.flybase.org/releases/FB2019_04/precomputed_files/genes/fbgn_fbtr_fbpp_fb_2019_04.tsv.gz", "/home/yixin/NutstoreFiles/Nutstore/Wu_lab/Snakemake_projects/miR983_975/output/expression/fbgn_fbtr_fbpp_fb_2019_04.tsv.gz")

# tx2gene expanded
download.file("ftp://ftp.flybase.org/releases/FB2019_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2019_04.tsv.gz", "/home/yixin/NutstoreFiles/Nutstore/Wu_lab/Snakemake_projects/miR983_975/output/expression/fbgn_fbtr_fbpp_expanded_fb_2019_04.tsv.gz")


