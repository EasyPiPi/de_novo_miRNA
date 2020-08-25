library(tidyverse)
library(topGO)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!"org.Dm.eg.db" %in% rownames(installed.packages())) BiocManager::install("org.Dm.eg.db")
library(org.Dm.eg.db)

root_dir <- '/home/yixin/Desktop/github_repo/de_novo_miRNA'

deseq2_dir <- file.path(root_dir, 'outputs/expression/table')

output_df_dir <- file.path(root_dir, "outputs/GO/table")
dir.create(output_df_dir, showWarnings = FALSE, recursive = TRUE)
output_fig_dir <- file.path(root_dir, 'outputs/GO/figure')
dir.create(output_fig_dir, showWarnings = FALSE, recursive = TRUE)

# read in tables
files <- list.files(deseq2_dir, pattern = ".csv")
df_names <- str_split(files, "_DE", simplify = T)[,1]

deseq2_dfs <- map(
    map(file.path(deseq2_dir, files),
        read_csv, col_types = cols(log2FoldChange = col_double())),
    ~ .x %>% dplyr::rename(gene_ID = X1)
) %>%
    set_names(df_names)

deseq2_dfs <- tibble(df_name = df_names, deseq2_df = deseq2_dfs)
deseq2_dfs <-
    deseq2_dfs %>% separate(df_name, into = c("species", "miRNA"), sep = "_")

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

# deseq2_dfs <- deseq2_dfs %>%
#     mutate(gores_padj = pmap(list(godata, result_fisher, result_fisherw01),
#                              function(a, b, c) GenTable(object = a, classicFisher = b, weight01Fisher = c,
#                                       orderBy = "weight01Fisher",
#                                       ranksOf = "classicFisher", topNodes = 20)))

deseq2_dfs <- deseq2_dfs %>%
    mutate(gores_padj = map2(godata, result_fisherw01,
                             function(a, c) GenTable(object = a, weight01Fisher = c,
                                                        orderBy = "weight01Fisher",
                                                        topNodes = 20)))

deseq2_dfs %>% dplyr::select(species, miRNA, gores_padj) %>%
    unnest(cols = c(gores_padj)) %>%
    write_csv(file.path(output_df_dir, "GO.csv"))

# showSigOfNodes(deseq2_dfs$godata[[1]], score(deseq2_dfs$result_fisherw01[[1]]), firstSigNodes = 5, useInfo = 'all')

