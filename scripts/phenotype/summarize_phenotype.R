# suppress all warnings
options(warn=-1)

suppressMessages(library(tidyverse))
library(broom)
library(rcompanion)
suppressMessages(library(cowplot))
#install.packages("magick")
suppressMessages(library(magick))

# root_dir <- "/home/yixin/Desktop/github_repo/de_novo_miRNA"
# output_dir <- "."

root_dir <- snakemake@params[["root_dir"]]
output_dir <- snakemake@params[["figure_out_dir"]]

####
phenotype_dir <- file.path(root_dir, 'data/phenotype')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(phenotype_dir, full.names = T)
file_names <- str_remove(list.files(phenotype_dir), ".csv")

#### Mating success ####
mating_success_dfs <-
    map(str_subset(files, "mating_success"),
        read_csv, col_types = cols(X1 = col_character())) %>%
    set_names(str_subset(file_names, "mating_success"))

# Process df
process_mating_df <- function(df) {
    df <- df %>%
        rename(genotype = X1, unsuccess = "No success") %>%
        gather("mating_success", "number_of_pairs", -genotype) %>%
        mutate(genotype = case_when(
            genotype == "Dmel-975KO" ~ "miR-975 KO",
            genotype == "Dmel-983KO" ~ "miR-983 KO",
            genotype == "Dsim-975KO" ~ "miR-975 KO",
            genotype == "Dsim-983KO" ~ "miR-983 KO",
            TRUE ~ "wildtype"
        )) %>%
        mutate(genotype = fct_relevel(genotype,
                                      "wildtype", "miR-983 KO", "miR-975 KO"),
               mating_success = fct_relevel(mating_success, "unsuccess", "success"))
    return(df)
}

# Plots
mating_success_plot <- function(df, rm_legend = TRUE){
    p <- df %>% ggplot(aes(y = number_of_pairs, fill = mating_success, x = genotype)) +
        geom_col(position="fill") +
        geom_text(aes(label = number_of_pairs), position = position_fill(vjust=0.5)) +
        labs(y = "Percentage of flies", fill = "Male mating", x = "") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_brewer(palette = "OrRd") +
        theme_cowplot() # + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    if (rm_legend) p <- p + theme(legend.position="None")
    if (!rm_legend) p <- p + theme(legend.position="right")

    p

}

mating_success_miR983_plots <- mating_success_dfs %>%
    map(process_mating_df) %>%
    map(~ .x %>% filter(genotype %in% c("wildtype", "miR-983 KO"))) %>%
    map(mating_success_plot)

mating_success_miR975_plots <- mating_success_dfs %>%
    map(process_mating_df) %>%
    map(~ .x %>% filter(genotype %in% c("wildtype", "miR-975 KO"))) %>%
    map(mating_success_plot)

# Statistical test
ratio_stat_test <- function(df){
    df <- df %>%
        rename(genotype = X1, unsuccess = "No success") %>%
        column_to_rownames("genotype")
    df %>% as.matrix() %>%
        pairwiseNominalIndependence(method = "none", fisher = TRUE, gtest = FALSE, chisq = FALSE)
}

map(mating_success_dfs, ratio_stat_test)

#### Suppression_of_female_remating ####
suppression_remating_dfs <-
    map(str_subset(files, "suppression_of_female_remating"),
        read_csv, col_types = cols(X1 = col_character())) %>%
    set_names(str_subset(file_names, "suppression_of_female_remating"))

process_remating_df <- function(df) {
    df <- df %>%
        rename(genotype = X1,
               rejection = "No success",
               acception = "success"
        ) %>%
        gather("mating_success", "number_of_pairs", -genotype) %>%
        mutate(genotype = case_when(
            genotype == "Dmel-975KO" ~ "miR-975 KO",
            genotype == "Dmel-983KO" ~ "miR-983 KO",
            genotype == "Dsim-975KO" ~ "miR-975 KO",
            genotype == "Dsim-983KO" ~ "miR-983 KO",
            TRUE ~ "wildtype"
        )) %>%
        mutate(genotype = fct_relevel(genotype,
                                      "wildtype", "miR-983 KO", "miR-975 KO"),
               mating_success = fct_relevel(mating_success, "acception", "rejection"))
    return(df)
}

supression_remating_plot <- function(df, rm_legend = TRUE){
    p <- df %>% ggplot(aes(y = number_of_pairs, fill = mating_success, x = genotype)) +
        geom_col(position="fill") +

        geom_text(aes(label = number_of_pairs), position = position_fill(vjust=0.5)) +
        labs(y = "Percentage of flies", fill = "Remating", x = "") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_brewer(palette = "BuGn") +
        theme_cowplot() # + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    if (rm_legend) p <- p + theme(legend.position="None")
    if (!rm_legend) p <- p + theme(legend.position="right")

    p

}

female_remating_miR983_plots <- suppression_remating_dfs %>%
    map(process_remating_df) %>%
    map(~ .x %>% filter(genotype %in% c("wildtype", "miR-983 KO"))) %>%
    map(supression_remating_plot)

female_remating_miR975_plots <- suppression_remating_dfs %>%
    map(process_remating_df) %>%
    map(~ .x %>% filter(genotype %in% c("wildtype", "miR-975 KO"))) %>%
    map(supression_remating_plot)

map(suppression_remating_dfs, ratio_stat_test)

#### Male fertility ####
fertility_dfs <-
    map(str_subset(files, "fertility"),
        read_csv, col_types = cols()) %>%
    set_names(str_subset(file_names, "fertility"))

fertility_melt <- function(df) {
    df <- df %>%
        gather("key", "number_of_offspring") %>%
        separate(key, into = c("species", "genotype", "round_of_mating"), sep = "_") %>%
        mutate(genotype = case_when(
            genotype == "miR975KO" ~ "miR-975 KO",
            genotype == "miR983KO" ~ "miR-983 KO",
            TRUE ~ "wildtype"
        )) %>%
        mutate(genotype = fct_relevel(genotype,
                                      "wildtype", "miR-983 KO", "miR-975 KO"))
}

fertility_dfs_melt <- map(fertility_dfs, fertility_melt)

# Plots
fertility_plot <- function(df, rm_legend = TRUE){
    p <- df %>%
        ggplot(aes(x = round_of_mating, y = number_of_offspring, fill = genotype)) +
        geom_boxplot() +
        scale_fill_manual(values=c("wildtype" = "#e41a1c",
                                   "miR-983 KO" = "#377eb8",
                                   "miR-975 KO" = "#4daf4a")) +
        labs(x = "Round of mating", y = "Number of progeny") +
        theme_cowplot()
    if (rm_legend) {p <- p + theme(legend.position="none")
    } else { p <- p + theme(legend.position="right")
    }

    return(p)
}

fertility_plots <- map2(fertility_dfs_melt, rep(TRUE, 4), fertility_plot)

# Statistical test
mean_stat_test <- function(df){
    df <- df %>% group_by(round_of_mating)
    genotype_n <- length(unique(df$genotype))

    if (genotype_n == 2) {
        # Using t.test for 2 genotypes
        df %>% do(t.test(number_of_offspring ~ genotype, data = .) %>% tidy())
    } else if (genotype_n > 2) {
        # Using anova for 3 genotypes
        df %>% do(aov(number_of_offspring ~ genotype, data = .) %>% tidy())
    }
}

map(fertility_dfs_melt, mean_stat_test)

#### Sperm competition ####
sperm_competetion_dfs <-
    map(str_subset(files, "sperm_competetion"), read_csv, col_types = cols()) %>%
    set_names(str_subset(file_names, "sperm_competetion")) %>%
    map(~.x %>% gather("genotype", "offspring_ratio"))

# process df
process_sperm_df <- function(df) {
    df <- df %>%
        mutate(genotype = case_when(
            genotype == "Dmel-975KO" ~ "miR-975 KO",
            genotype == "Dmel-983KO" ~ "miR-983 KO",
            genotype == "Dsim-975KO" ~ "miR-975 KO",
            genotype == "Dsim-983KO" ~ "miR-983 KO",
            TRUE ~ "wildtype"
        )) %>%
        mutate(genotype = fct_relevel(genotype,
                                      "wildtype", "miR-983 KO", "miR-975 KO"))
    return(df)
}

# Plots
sperm_competetion_plot <- function(df, cat_string, rm_legend = TRUE){
    p <- df %>% ggplot(aes(x = genotype, y = offspring_ratio, fill = genotype)) +
        geom_boxplot() +
        scale_fill_brewer(palette = "Set1") +
        labs(y = paste0(cat_string, " score"), x = "") +
        guides(fill = guide_legend(title = "Genotype")) +
        theme_cowplot() # + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    if (rm_legend) p <- p + theme(legend.position="none")
    p
}

sperm_competetion_miR983_plots <- sperm_competetion_dfs %>%
    map(process_sperm_df) %>%
    map(~ .x %>% filter(genotype %in% c("wildtype", "miR-983 KO"))) %>%
    map2(rep(c("Defense", "Offense"), 2), sperm_competetion_plot)

sperm_competetion_miR975_plots <- sperm_competetion_dfs %>%
    map(process_sperm_df) %>%
    map(~ .x %>% filter(genotype %in% c("wildtype", "miR-975 KO"))) %>%
    map2(rep(c("Defense", "Offense"), 2), sperm_competetion_plot)

# Statistical test
mean_stat_test2 <- function(df){
    # aov_pvalue <- aov(offspring_ratio ~ genotype, data = df) %>% tidy() %>%
    #     filter(term == "genotype") %>% pull(p.value)
    #
    # if (aov_pvalue >= 0.05){
    #     return(aov_pvalue)
    # } else{
    pairwise.wilcox.test(df$offspring_ratio, df$genotype, p.adjust.method = "none") %>% tidy()
    # }
}

map(sperm_competetion_dfs, mean_stat_test2)

###############################################################################
# Cow plot
# Male fertility + female remating
miR_deletion <- ggdraw() + draw_image(image_read(file.path(root_dir, "data/miR_sequence/miRNAs_deletion.png")))

mating_legend <- get_legend(
    mating_success_dfs$dmel_male_mating_success_with_non_virgin_female %>%
        process_mating_df() %>% mating_success_plot(FALSE)
)

remating_legend <-
    get_legend(
        suppression_remating_dfs$dmel_suppression_of_female_remating %>%
            process_remating_df() %>% supression_remating_plot(FALSE)
    )

col1 <- plot_grid(mating_success_miR975_plots$dmel_male_mating_success_with_virgin_female,
                  mating_success_miR975_plots$dsim_male_mating_success_with_virgin_female,
                  mating_success_miR983_plots$dmel_male_mating_success_with_virgin_female,
                  mating_success_miR983_plots$dsim_male_mating_success_with_virgin_female,
                  ncol = 1)

col2 <- plot_grid(mating_success_miR975_plots$dmel_male_mating_success_with_non_virgin_female,
                  mating_success_miR975_plots$dsim_male_mating_success_with_non_virgin_female,
                  mating_success_miR983_plots$dmel_male_mating_success_with_non_virgin_female,
                  mating_success_miR983_plots$dsim_male_mating_success_with_non_virgin_female,
                  ncol = 1)

legend1 <- plot_grid(mating_legend, mating_legend, mating_legend, mating_legend, ncol = 1)

col3 <- plot_grid(female_remating_miR975_plots$dmel_suppression_of_female_remating,
                  female_remating_miR975_plots$dsim_suppression_of_female_remating,
                  female_remating_miR983_plots$dmel_suppression_of_female_remating,
                  female_remating_miR983_plots$dsim_suppression_of_female_remating,
                  ncol = 1)

legend2 <- plot_grid(remating_legend, remating_legend, remating_legend, remating_legend, ncol = 1)

col4 <- plot_grid(sperm_competetion_miR975_plots$dmel_sperm_competetion_offense,
                  sperm_competetion_miR975_plots$dsim_sperm_competetion_offense,
                  sperm_competetion_miR983_plots$dmel_sperm_competetion_offense,
                  sperm_competetion_miR983_plots$dsim_sperm_competetion_offense,
                  ncol = 1)

phenotype_fig <- plot_grid(col1, col2, legend1, col3, legend2, col4, nrow = 1,
                           rel_widths = c(1, 1, 0.4, 1, 0.5, 1))

phenotype_fig <- plot_grid(miR_deletion, phenotype_fig, ncol = 1,
                           rel_heights = c(0.4, 1))

ggsave(plot = phenotype_fig, file.path(output_dir, "phenotype.png"),
       width = 34, height = 40, units = "cm")

# fertility_legend <- get_legend(sperm_competetion_plot(sperm_competetion_dfs$dsim_sperm_competetion_defense, "", FALSE))
# mating_legend <- get_legend(
#     mating_success_plot(mating_success_dfs$dmel_male_mating_success_with_non_virgin_female, FALSE))
# remating_legend <- get_legend(supression_remating_plot(suppression_remating_dfs$dmel_suppression_of_female_remating, FALSE))
#
# upper <-
#     plot_grid(mating_success_plots$dmel_male_mating_success_with_virgin_female,
#               mating_success_plots$dsim_male_mating_success_with_virgin_female,
#               mating_legend,
#               mating_success_plots$dmel_male_mating_success_with_non_virgin_female,
#               mating_success_plots$dsim_male_mating_success_with_non_virgin_female,
#               mating_legend,
#               nrow = 1, labels = c("A","B", "", "C","D"), rel_widths = c(2,2,1,2,2,1))
#
# bottom <- plot_grid(female_remating_plots$dmel_suppression_of_female_remating,
#                  female_remating_plots$dsim_suppression_of_female_remating,
#                  remating_legend,
#                  sperm_competetion_plots$dmel_sperm_competetion_offense,
#                  sperm_competetion_plots$dsim_sperm_competetion_offense,
#                  fertility_legend,
#                  nrow = 1, labels = c("E","F","", "G", "H"), rel_widths = c(2,2,1,2,2,1))
#
# phenotype_fig <- plot_grid(upper, bottom, ncol = 1)
# ggsave(plot = phenotype_fig, file.path(output_dir, "phenotype.png"),
#        width = 32, height = 18, units = "cm")
#
# # supplementary
# sperm_defense <-
#     plot_grid(sperm_competetion_plots$dmel_sperm_competetion_defense,
#               sperm_competetion_plots$dsim_sperm_competetion_defense,
#               fertility_legend,
#               nrow = 1, labels = c("A", "B"), rel_widths = c(2, 2, 1))
#
# ggsave(plot = sperm_defense, file.path(output_dir, "sperm_competetion_defense.png"),
#        width = 15, height = 9, units = "cm")
#
#
# male_fertility <- plot_grid(fertility_plots$dmel_fertility_miR983,
#                             fertility_plots$dsim_fertility_miR983,
#                             fertility_plots$dmel_fertility_miR975,
#                             fertility_plots$dsim_fertility_miR975,
#                             nrow = 2, labels = c("A","B","C","D"))
#
# male_fertility <- plot_grid(male_fertility, fertility_legend, rel_widths = c(4, 1))
#
# ggsave(plot = male_fertility, file.path(output_dir, "male_fertility.png"),
#        width = 15, height = 15, units = "cm")

