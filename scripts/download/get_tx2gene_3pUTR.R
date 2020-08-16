library(biomaRt)
library(tidyverse)

root_dir <- '/home/yixin/Desktop/github_repo/de_novo_miRNA'
utr_output_dir <- file.path(root_dir, "Snakemake_projects/miR983_975/output/resources/3utr")
tx2gene_output_dir <- file.path(root_dir, "Snakemake_projects/miR983_975/output/resources/tx2gene")
ortholog_output_dir <- file.path(root_dir, "Snakemake_projects/miR983_975/output/resources/ortholog")

# connect to metazoa mart
mart = useEnsembl("metazoa_mart", host="metazoa.ensembl.org")
marts <- listDatasets(mart)

# genome version "Drosophila melanogaster genes (BDGP6.22)" and "Drosophila simulans genes (ASM75419v3)"
dataset_names <- c("dmelanogaster_eg_gene", "dsimulans_eg_gene")

# get tx/gene name and other info
get_tx2gene <- function(dataset) {
    tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "external_gene_name", "transcript_biotype"), mart = dataset)
}

# get tx2gene dfs
tx2gene_ls <- list()
for (dataset_name in dataset_names){
    dataset <- useDataset(dataset_name, mart = mart)
    tx2gene_ls[[str_replace(dataset_name, "_eg_gene", "_tx2gene")]] <- get_tx2gene(dataset)
}

# save files
mapply(write_csv, tx2gene_ls, file.path(tx2gene_output_dir, str_c(names(tx2gene_ls), ".csv")))

# get sequence from biomart
get_sequence <- function(dataset, sequence_type) {
    getBM(attributes = c(sequence_type, "ensembl_transcript_id"), mart = dataset)
}

# get 3'UTR sequence
utr3_ls <- list()
for (dataset_name in dataset_names){
    dataset <- useDataset(dataset_name, mart = mart)
    # ensembl_attributes <- listAttributes(dataset)
    utr3_ls[[str_replace(dataset_name, "_eg_gene", "_3utr")]] <- get_sequence(dataset, "3utr")
}

# add species id to satisfy targetScan
utr3_ls$dmelanogaster_3utr$species_id <- "7227"
utr3_ls$dsimulans_3utr$species_id <- "7240"

utr3_ls <- lapply(utr3_ls, function(df) df[c("ensembl_transcript_id", "species_id", "3utr")])
# save utr files
# mapply(exportFASTA, utr3_ls, file.path(utr_output_dir, str_c(names(utr3_ls), ".fasta")))

mapply(write_delim, utr3_ls, file.path(utr_output_dir, str_c(names(utr3_ls), ".tab")), delim = "\t", col_names = F)

# get Othorlogy
dataset_name <- dataset_names[[2]]
ensembl_attributes <- listAttributes(dataset)
dataset <- useDataset(dataset_name, mart = mart)
dsim2dmel_ortholog <- getBM(attributes = c("ensembl_gene_id", "dmelanogaster_eg_homolog_ensembl_gene", "dmelanogaster_eg_homolog_orthology_type", "dmelanogaster_eg_homolog_orthology_confidence"), mart = dataset)
dsim2dmel_1to1_ortholog <- dsim2dmel_ortholog %>%
    filter(dmelanogaster_eg_homolog_orthology_type == "ortholog_one2one", dmelanogaster_eg_homolog_orthology_confidence == 1)

write_csv(dsim2dmel_ortholog, file.path(ortholog_output_dir, "dsim2dmel_ortholog.csv"))
write_csv(dsim2dmel_1to1_ortholog, file.path(ortholog_output_dir, "dsim2dmel_1to1_ortholog.csv"))

