suppressMessages(suppressWarnings(library(tidyverse)))

# rtholog_in <- "/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/flybase/dmel_orthologs_in_drosophila_species_fb_2020_03.tsv"
ortholog_in <- snakemake@input[["dmel_orthologs"]]
ortholog_out <- snakemake@output[["dmel_dsim_orthologs"]]

col_names <- c("FBgn_ID", "GeneSymbol", "Arm/Scaffold", "Location",
               "Strand", "Ortholog_FBgn_ID", "Ortholog_GeneSymbol",
               "Ortholog_Arm/Scaffold", "Ortholog_Location",
               "Ortholog_Strand", "OrthoDB_Group_ID")

ortholog <-
    read_delim(ortholog_in, delim = "\t",
               col_names = col_names, comment = "#",
               col_types = cols(FBgn_ID = col_character()))

ortholog %>%
    filter(str_detect(Ortholog_GeneSymbol, "Dsim")) %>%
    select(Ortholog_FBgn_ID, FBgn_ID) %>%
    unique() %>%
    write_csv(ortholog_out)
