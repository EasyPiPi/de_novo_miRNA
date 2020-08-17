library(tidyverse)
library(skimr)
library(VennDiagram)
library(cowplot)

# read in dfs and generate dirs
source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "read_in_deseq2_dfs.R"))
source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "misregulated_genes_function.R"))

resource_dir <- file.path(root_dir, 'Snakemake_projects/miR983_975/output/resources')
sebida_dir <- file.path(root_dir, "Snakemake_projects/miR983_975/data/sebida")
testes_exp_dir <- file.path(root_dir, "Snakemake_projects/miR983_975/data/flybase_testes_expressed")

deseq2_dfs <- deseq2_dfs %>% mutate(deseq2_df = map(deseq2_df, get_category))
deseq2_dfs <- deseq2_dfs %>% arrange(desc(miRNA), species)

# mis-regulated gene chromosome location
tx2gene <- read_csv(file.path(resource_dir, "tx2gene/dmelanogaster_tx2gene.csv")) 
chr_loc <- tx2gene %>%
    select(ensembl_gene_id, chromosome_name) %>%
    unique()

deseq2_dfs <- deseq2_dfs %>% 
    mutate(chr_loc_tab = map(deseq2_df, get_chr_loc_tab, chr_loc = chr_loc)) %>% 
    mutate(chr_loc_p = map_dbl(chr_loc_tab, get_fisher_pvalue))

# sex bias genes
sebida_dmel <- 
    read_delim(file.path(sebida_dir, "sebida_melanogaster_3.2.txt"), delim = "\t", col_types = cols(`MetaM/F` = col_character(), MetaP = col_character(), MetaFDR = col_character(), ENC = col_character(), FOP = col_character())) %>%
    rename("gene_ID" = `Current FBgn`) %>%
    select(gene_ID, MetaClass)

deseq2_dfs <- deseq2_dfs %>% 
    mutate(sexbias_tab = map(deseq2_df, get_sexbias_tab, sebida = sebida_dmel),
           sexbias_p = map_dbl(sexbias_tab, get_fisher_pvalue),
           sexbias_updown_tab = map(deseq2_df, get_sexbias_updown_tab, sebida = sebida_dmel),
           sexbias_updown_p = map_dbl(sexbias_tab, get_fisher_pvalue)
           )

# testes specific expressed genes
testes_exp_genes <- scan(file.path(testes_exp_dir, "FlyBase_IDs.txt"), what = character())

deseq2_dfs <- deseq2_dfs %>%
    mutate(testes_tab = map(deseq2_df, get_testes_tab),
           testes_p = map_dbl(testes_tab, get_fisher_pvalue))

# Visualization
deseq2_dfs <- deseq2_dfs %>% 
    mutate(fig_prefix = str_c(species, "-", str_replace(miRNA, "miR", "mir-"), " KO"))

# cowplot
chr_plots <- pmap(list(deseq2_dfs$chr_loc_tab, deseq2_dfs$fig_prefix, c(rep(TRUE, 4))), get_chr_loc_plt)
chr_plot_legend <- get_legend(get_chr_loc_plt(deseq2_dfs$chr_loc_tab[[1]], "", FALSE))

sex_plots <- map2(deseq2_dfs$sexbias_tab, c(rep(TRUE, 4)), get_sexbias_plt)
sex_plot_legend <- get_legend(get_sexbias_plt(deseq2_dfs$sexbias_tab[[1]], FALSE))

upper_row <- plot_grid(plotlist = chr_plots, nrow = 1)
bottom_row <- plot_grid(plotlist = sex_plots, nrow = 1)
upper_row <- plot_grid(upper_row, chr_plot_legend, rel_widths = c(10, 2))
bottom_row <- plot_grid(bottom_row, sex_plot_legend, rel_widths = c(10, 2))

plot_grid(upper_row, bottom_row, ncol = 1,
          labels = c("A", "B"), rel_heights = c(1, 1.2))

ggsave(file.path(output_dir, "misregulated_gene_chr_loc_sex_bias.png"), width = 10, height = 6)

walk2(deseq2_dfs$sexbias_updown_tab, str_c(deseq2_dfs$fig_prefix, "_sex_bias_updown.png"), get_sexbias_updown_plt)
walk2(deseq2_dfs$testes_tab, str_c(deseq2_dfs$fig_prefix, "_testes_specific.png"), get_testes_plt)

# Compare mis-regulated genes due to different miRNA KOs
misregulated_genes <- map(deseq2_dfs$deseq2_df, 
    ~ .x %>% filter(significant == "Misregulated") %>% pull(gene_ID))
names(misregulated_genes) <- deseq2_dfs$fig_prefix

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
    x = misregulated_genes[1:2],
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
    x = misregulated_genes[3:4],
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
