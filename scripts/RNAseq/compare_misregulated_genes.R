suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(VennDiagram)))

# root_dir <- "/home/yixin/Desktop/github_repo/de_novo_miRNA"
# ortholog_in <- "/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/flybase/dmel_dsim_orthologs.csv"
# output_dir <- "."

root_dir <- snakemake@params[["root_dir"]]
output_dir <- snakemake@params[["figure_out_dir"]]
ortholog_in <- snakemake@input[["ortholog"]]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
deseq2_dir <- file.path(root_dir, 'outputs/expression/table')

# read ortholog table
ortholog <-
    read_csv(ortholog_in, col_types = cols(Ortholog_FBgn_ID = col_character()))

one2mul <- ortholog %>% count(Ortholog_FBgn_ID) %>% filter(n != 1) %>% pull(Ortholog_FBgn_ID)
mul2one <-
    ortholog %>%
    filter(!Ortholog_FBgn_ID %in% one2mul) %>%
    count(FBgn_ID) %>% filter(n != 1) %>% pull(FBgn_ID)

ortholog_1to1 <-
    ortholog %>%
    filter(!Ortholog_FBgn_ID %in% one2mul, !FBgn_ID %in% mul2one)

# read DESeq2 tables
files <- list.files(deseq2_dir, pattern = ".csv")
df_names <- str_split(files, "_DE", simplify = T)[,1]

deseq2_dfs <- map(
    map(file.path(deseq2_dir, files),
        read_csv, col_types = cols(log2FoldChange = col_double()), col_names = TRUE),
    ~ .x %>% dplyr::rename(gene_ID = X1)
) %>%
    set_names(df_names)

deseq2_dfs <- tibble(df_name = df_names, deseq2_df = deseq2_dfs)
deseq2_dfs <-
    deseq2_dfs %>% separate(df_name, into = c("species", "miRNA"), sep = "_")

# map D. simulans ids to D. melanogaster ids for GO analysis
deseq2_dfs_dsim <- deseq2_dfs %>%
    filter(species == "dsi") %>%
    dplyr::mutate(deseq2_df = map(deseq2_df,
                                  ~ .x %>%
                                      dplyr::left_join(ortholog_1to1, by = c("gene_ID" = "Ortholog_FBgn_ID")) %>%
                                      dplyr::select(-gene_ID) %>%
                                      dplyr::rename("gene_ID" = "FBgn_ID") %>%
                                      dplyr::filter_all(all_vars(!is.na(.)))))

# replace the dsim dfs
deseq2_dfs <- deseq2_dfs %>%
    filter(species == "dme") %>%
    bind_rows(deseq2_dfs_dsim) %>%
    mutate(deseq2_df = map(deseq2_df,
                           ~ .x %>% dplyr::filter_all(all_vars(!is.na(.)))))

rm(deseq2_dfs_dsim)

# define significantly mis-regulated genes
deseq2_dfs <-
    deseq2_dfs %>%
    mutate(deseq2_df = map(deseq2_df, ~ .x %>%
                               mutate(
                                   significant = if_else(padj < 0.05,
                                                         "misregulated",
                                                         "others")
                               )))

# Compare mis-regulated genes due to different miRNA KOs
misregulated_genes <- map(deseq2_dfs$deseq2_df,
    ~ .x %>% filter(significant == "misregulated") %>% pull(gene_ID))
names(misregulated_genes) <- paste(deseq2_dfs$species, deseq2_dfs$miRNA, sep = "_")

# suppress output the log file
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

venn_plot <- venn.diagram(
    x = misregulated_genes,
    filename = file.path(output_dir, "misregulated_gene_venn.tiff"),
    col = "transparent",
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white",
                  "white", "white", "white", "white", "darkblue", "white",
                  "white", "white", "white", "darkgreen", "white"),
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 270,
    margin = 0.2
)

# miR-983
venn.plot <- venn.diagram(
    x = with(misregulated_genes, list(dme_miR983, dsi_miR983)),
    filename = file.path(output_dir, "misregulated_gene_miR983_venn.tiff"),
    lwd = 4,
    fill = c("cornflowerblue", "darkorchid1"),
    alpha = 0.75,
    label.col = "white",
    cex = 3,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("cornflowerblue", "darkorchid1"),
    cat.cex = 1.5,
    cat.fontfamily = "serif",
    cat.fontface = "bold",
    cat.dist = c(0.03, 0.03),
    cat.pos = c(-20, 14),
    category.names = c(expression(italic('dme-mir-983')), expression(italic('dsi-mir-983')))
)

# miR-975
venn.plot <- venn.diagram(
    x = with(misregulated_genes, list(dme_miR975, dsi_miR975)),
    filename = file.path(output_dir, "misregulated_gene_miR975_venn.tiff"),
    lwd = 4,
    fill = c("cornflowerblue", "darkorchid1"),
    alpha = 0.75,
    label.col = "black",
    cex = 3,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("cornflowerblue", "darkorchid1"),
    cat.cex = 1.5,
    cat.fontfamily = "serif",
    cat.fontface = "bold",
    cat.dist = c(0.03, 0.03),
    cat.pos = c(-10, 0),
    category.names = c(expression(italic('dme-mir-975')), expression(italic('dsi-mir-975')))
)
