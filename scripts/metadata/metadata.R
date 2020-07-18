suppressWarnings(suppressMessages(library(tidyverse)))

metadata_dir <- path.expand("~/Desktop/github_repo/de_novo_miRNA/metadata/small_RNAseq")
gse_all_out <- path.expand("~/Desktop/github_repo/de_novo_miRNA/metadata/small_RNAseq/metadata.csv")
metadata_out_dir <- path.expand("~/Desktop/github_repo/de_novo_miRNA/outputs/metadata/small_RNAseq")
dir.create(metadata_out_dir, showWarnings = FALSE, recursive = TRUE)
gse_sample_num_out <- file.path(metadata_out_dir, "sample_number.csv")

# read in meta data frames from multiple studies
gse11624 <-
    read_csv(file.path(metadata_dir, "GSE11624", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse15898 <-
    read_csv(file.path(metadata_dir, "GSE15898", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse22067 <-
    read_csv(file.path(metadata_dir, "GSE22067", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse24304 <-
    read_csv(file.path(metadata_dir, "GSE24304", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse37041 <-
    read_csv(file.path(metadata_dir, "GSE37041", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse48017 <-
    read_csv(file.path(metadata_dir, "GSE48017", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse56244 <-
    read_csv(file.path(metadata_dir, "GSE56244", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse97364 <-
    read_csv(file.path(metadata_dir, "GSE97364", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse98013 <-
    read_csv(file.path(metadata_dir, "GSE98013", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))
gse110086 <-
    read_csv(file.path(metadata_dir, "GSE110086", "SraRunTable.txt"),
             col_types = cols(Run = col_character()))

shared_cols <-
    reduce(list(colnames(gse11624), colnames(gse15898), colnames(gse22067),
                colnames(gse24304), colnames(gse37041), colnames(gse48017),
                colnames(gse56244), colnames(gse97364), colnames(gse98013),
                colnames(gse110086)),
           intersect)

# select subset cols
sel_cols <- c("Run", "GEO_Accession (exp)", "Organism", "source_name", "tissue")
gse_all <- bind_rows(gse11624, gse15898, gse22067, gse24304, gse37041,
                     gse48017, gse56244, gse97364, gse98013, gse110086)

# gse_all <-
gse_all <- gse_all %>%
    select(all_of(sel_cols)) %>%
    rename(run = "Run", sample = "GEO_Accession (exp)", organism = "Organism") %>%
    mutate(tissue = case_when(
        (tissue == "testes" | tissue == "testis") ~ "testes",
        TRUE ~ "male_body"
    )) %>%
    mutate(species = case_when(
        organism == "Drosophila melanogaster" ~ "dme",
        organism == "Drosophila simulans" ~ "dsi",
        organism == "Drosophila sechellia" ~ "dse",
        organism == "Drosophila erecta" ~ "der",
        organism == "Drosophila virilis" ~ "dvi"

    )) %>% base::unique()

gse_sample_num <- gse_all %>%
    count(organism, tissue) %>%
    spread(key = tissue, value = n)

write_csv(gse_all, gse_all_out)
write_csv(gse_sample_num, gse_sample_num_out)
