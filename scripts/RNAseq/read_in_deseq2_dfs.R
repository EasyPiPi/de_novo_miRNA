library(tidyverse)

# root_dir <- "/Users/yizhao/我的坚果云/Wu_lab"
root_dir <- '/home/yixin/NutstoreFiles/Nutstore/Wu_lab'
deseq2_dir <- file.path(root_dir, 'Snakemake_projects/miR983_975/output/deseq2')
output_dir <- file.path(root_dir, 'Snakemake_projects/miR983_975/output/misregulated_genes/figure')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# read in tables
files <- list.files(deseq2_dir, pattern = ".csv")
df_names <- str_split(files, "_DE", simplify = T)[,1]

deseq2_dfs <- map(
    map(file.path(deseq2_dir, files), read_csv),
    ~ .x %>% dplyr::rename(gene_ID = X1)
) %>% 
    set_names(df_names)

deseq2_dfs <- tibble(df_name = df_names, deseq2_df = deseq2_dfs)
deseq2_dfs <- deseq2_dfs %>% separate(df_name, into = c("species", "miRNA"), sep = "_")