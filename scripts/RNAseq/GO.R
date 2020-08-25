suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(topGO)))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!"org.Dm.eg.db" %in% rownames(installed.packages())) BiocManager::install("org.Dm.eg.db")
suppressMessages(suppressWarnings(library(org.Dm.eg.db)))

# root_dir <- "/home/yixin/Desktop/github_repo/de_novo_miRNA"
# ortholog_in <- "/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/flybase/dmel_dsim_orthologs.csv"
# go_out <- file.path(output_df_dir, "GO.csv")

root_dir <- snakemake@params[["root_dir"]]
ortholog_in <- snakemake@input[["ortholog"]]
go_out <- snakemake@output[["go"]]

deseq2_dir <- file.path(root_dir, 'outputs/expression/table')

output_df_dir <- file.path(root_dir, "outputs/GO/table")
dir.create(output_df_dir, showWarnings = FALSE, recursive = TRUE)
# output_fig_dir <- file.path(root_dir, 'outputs/GO/figure')
# dir.create(output_fig_dir, showWarnings = FALSE, recursive = TRUE)

# read ortholog table
ortholog <-
    read_csv(ortholog_in, col_types = cols(Ortholog_FBgn_ID = col_character()))

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
                                      dplyr::left_join(ortholog, by = c("gene_ID" = "Ortholog_FBgn_ID")) %>%
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

# define p-value for filtering
padj_sel <- function (padj) {
    return(padj < 0.05)
}

get_GOdata <- function(df, geneSel = padj_sel) {
    all_genes <- df$padj
    names(all_genes) <- df$gene_ID
    all_genes <- na.omit(all_genes)

    geneID2GO <- mapIds(org.Dm.eg.db, keys = names(all_genes), column="GO",
                        keytype="FLYBASE", multiVals="list")

    GOdata <- new("topGOdata", ontology = "BP",
                  allGenes = all_genes,
                  geneSel = geneSel,
                  nodeSize = 10, annotationFun = annFUN.gene2GO, gene2GO = geneID2GO)
}

get_GO_res <- function(GOdata){

    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultFisher.w01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       weight01Fisher = resultFisher.w01, orderBy = "weight01Fisher",
                       ranksOf = "classicFisher", topNodes = 10)

    return(allRes)
}

deseq2_dfs <- deseq2_dfs %>%
    mutate(godata = map(deseq2_df, get_GOdata)) %>%
    mutate(result_fisher = map(godata, runTest, algorithm = "classic", statistic = "fisher"),
           result_fisherw01 = map(godata, runTest, algorithm = "weight01", statistic = "fisher"))

deseq2_dfs <- deseq2_dfs %>%
    mutate(gores_padj = map2(godata, result_fisherw01,
                             function(a, c) GenTable(object = a, weight01Fisher = c,
                                                        orderBy = "weight01Fisher",
                                                        topNodes = 20)))

deseq2_dfs %>% dplyr::select(species, miRNA, gores_padj) %>%
    unnest(cols = c(gores_padj)) %>%
    write_csv(go_out)
