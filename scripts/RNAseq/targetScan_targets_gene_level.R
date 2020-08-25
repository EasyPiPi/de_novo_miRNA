suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(naniar)))
suppressMessages(suppressWarnings(library(grid)))
# library(skimr)

# root_dir <- '/home/yixin/Desktop/github_repo/de_novo_miRNA'
root_dir <- snakemake@params[["root_dir"]]

deseq2_dir <- file.path(root_dir, 'outputs/expression/table')
fig_output_dir <- file.path(root_dir, 'outputs/miRNA_targets/figure/gene_level')
dir.create(fig_output_dir, showWarnings = F, recursive = T)

# read in tables
files <- list.files(deseq2_dir, pattern = ".csv")
df_names <- str_split(files, "_DE", simplify = T)[,1]

deseq2_dfs <- map(
    map(file.path(deseq2_dir, files),
        read_csv, col_types = c(X1 = col_character())),
    ~ .x %>% rename(gene_ID = X1)
) %>%
    set_names(df_names)

deseq2_dfs <- tibble(df_name = df_names, deseq2_df = deseq2_dfs)
deseq2_dfs <- deseq2_dfs %>% separate(df_name, into = c("species", "miRNA"), sep = "_")

################################################################################
# miRNA targets
species <- c("dme", "dsi")
targetscan_files <-
    file.path(root_dir, "outputs/miRNA_targets/table", species, "targetscan_targets.tab")

targetscan <-
    tibble(
        species = species,
        df = map(targetscan_files, read_delim,
                 delim = "\t", col_types = cols(a_Gene_ID = col_character()))
    )

targetscan <-
    targetscan %>%
    mutate(df = map(df, ~ .x %>%
                        rename(miRNA_in_this_species = "miRNA in this species")))

tx2gene <-
    tibble(
        species = species,
        tx2gene = map(file.path(root_dir, "external_resources", species, "tx2gene.csv"),
                      read_csv, col_types = cols(transcript_id = col_character())))

targetscan <-
    targetscan %>%
    left_join(tx2gene, by = "species") %>%
    mutate(df = map2(df, tx2gene, ~ .x %>% left_join(.y, by = c("a_Gene_ID" = "transcript_id"))))

get_targetNum <- function(df) {
    df %>%
        select(miRNA_family_ID, gene_id, Site_type) %>%
        # filter target sites
        # filter_at("Site_type", ~.x %in% c("8mer-1a")) %>%
        group_by(gene_id, miRNA_family_ID) %>%
        count(miRNA_family_ID) %>%
        spread(miRNA_family_ID, n) %>%
        ungroup()
}

# get miRNA target number per genes
targetscan_target <- targetscan %>%
    mutate(targetNum = map(df, get_targetNum)) %>%
    select(species, targetNum)

# group target tyoe by target site number
get_target_site_gp <- function(df){
    df %>% mutate_at(vars(contains("3p"), contains("5p")), ~ case_when(
        is.na(.) ~ "non-target",
        TRUE ~ "target"
    ))
}

plot_dfs <- deseq2_dfs %>%
    left_join(targetscan_target, by = "species") %>%
    mutate(miRNA_num = str_sub(miRNA, 4)) %>%
    mutate(plot_df = map2(deseq2_df, targetNum, left_join, by = c("gene_ID" = "gene_id"))) %>%
    mutate(plot_df = map(plot_df, get_target_site_gp),
           plot_df = map(plot_df, function(df) {
               colnames(df) <- str_replace(colnames(df), "mir", "miR")
               return(df)
           } ))

# CDF plot
get_plots <- function(df, species, miRNA, miRNA_num) {
    subdf <- df %>%
        # filter(baseMean > 10) %>%
        # select mutant miRNA for visualization
        select(log2FoldChange, contains(miRNA_num))

    for (miR_name in colnames(subdf)[2:length(colnames(subdf))]) {
        print(paste(species, miRNA, "KO", miR_name, "target"))
        target_vals <- filter(subdf, subdf[[miR_name]] == "target")$log2FoldChange
        nontarget_vals <- filter(subdf, subdf[[miR_name]] == "non-target")$log2FoldChange
        print(paste("target number:", length(target_vals)))
        print(paste("non-target number:", length(nontarget_vals)))
        print(paste("p value:", ks.test(target_vals, nontarget_vals)$`p.value`))
        # CDF plot
        fig_1 <- subdf %>%
            ggplot(aes(x = log2FoldChange, color = as.factor(subdf[[miR_name]]))) +
            stat_ecdf() +
            xlim(-0.5, 0.5) +
            labs(
                x = expression(log[2]*"(foldchange)"),
                color = str_replace_all(miR_name, "_", "-"),
                y = "CDF"
            ) +
            theme_classic() +
            theme(legend.position = c(0.8, 0.2),
                  text = element_text(size = 12))
        ggsave(fig_1,
               filename = file.path(
                   fig_output_dir,
                   str_c(species, "_", miRNA, "KO_", miR_name, "_target_CDF.png")
               ), width = 9, height = 8, units = "cm"
        )

        fig_2 <- subdf %>%
            ggplot(aes(
                y = log2FoldChange,
                x = as.factor(subdf[[miR_name]]),
                color = as.factor(subdf[[miR_name]])
            )) +
            geom_violin() +
            geom_boxplot(width=0.1) +
            labs(y = expression(log[2]*"(foldchange)"),
                 color = str_replace_all(miR_name, "_", "-")) +
            ylim(-0.5, 0.5) +
            annotate("text", x = 2.7, y = -0.4, size = 3, hjust = 0,
                     label = paste("non-target n =", length(nontarget_vals))) +
            annotate("text", x = 2.7, y = -0.45, size = 3, hjust = 0,
                     label = paste("target n =", length(target_vals))) +
            annotate("text", x = 2.7, y = -0.5, size = 3, hjust = 0,
                     label = "KS test p-value < 1e-5") +
            coord_cartesian(xlim = c(1, 2),
                            clip = 'off') +
            theme_classic() +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )
        ggsave(fig_2,
               filename = file.path(
                   fig_output_dir,
                   str_c(species, "_", miRNA, "KO_", miR_name, "_target_boxplot.png")
               ), width = 9, height = 8, units = "cm")
    }
}

pwalk(list(plot_dfs$plot_df, plot_dfs$species, plot_dfs$miRNA, plot_dfs$miRNA_num), get_plots)
