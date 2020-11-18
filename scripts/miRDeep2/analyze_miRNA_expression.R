suppressMessages(suppressWarnings(library(tidyverse)))

#### read in files ####
# fig_out <- '/home/yixin/Desktop/github_repo/de_novo_miRNA/outputs/miRDeep2/figure'
# dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)
# exp_dir <- '/home/yixin/Desktop/github_repo/de_novo_miRNA/outputs/miRDeep2/quantifier'
# metadata_df <-
#     read_csv(
#         '/home/yixin/Desktop/github_repo/de_novo_miRNA/metadata/small_RNAseq/metadata.csv',
#         col_types = cols(run = col_character())
#     )
# orth_df <-
#     readxl::read_xlsx(
#         '/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/Mohammed_et_al/Supplemental_Table_S4.xlsx'
#     )
# class_df <-
#     read_csv(
#         '/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/Lyu_et_al/table_S8.csv'
#     )

fig_out <- snakemake@params[["figure_out_dir"]]
dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)
exp_dir <- snakemake@params[["miR_expression_dir"]]

metadata_df <-
    read_csv(snakemake@input[["metadata"]], col_types = cols(run = col_character()))

orth_df <- readxl::read_xlsx(snakemake@input[["ortholog"]])
class_df <- read_csv(snakemake@input[["classification"]],
                     col_types = cols(miRNA = col_character()))

#### get ortholog names ####
orth_df <- orth_df %>% select(familyname, OrthologName)

class_df <- class_df %>%
    left_join(orth_df, by = c("miRNA" = "familyname")) %>%
    mutate(species = str_sub(OrthologName, 1, 3)) %>%
    spread(species, OrthologName) %>%
    select("miRNA", "Age", "Evolutionary mode", dme, dsi, dse, der, dvi)

class_df$dme <- class_df$miRNA

# dplyr mutate/replace several columns on a subset of rows
# https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
    condition <- eval(substitute(condition), .data, envir)
    .data[condition, ] <- .data[condition, ] %>% mutate(...)
    .data
}

# set columns to NA where the name mapping might be problematic
class_df <- class_df %>%
    mutate_cond(Age == "0 - 4 Myrs", der = NA, dvi = NA) %>%
    mutate_cond(Age == "4 - 30 Myrs", dvi = NA)

class_num <- class_df %>%
    filter(!across(c(dme, dsi, dse, der, dvi), is.na)) %>%
    count(Age, `Evolutionary mode`)

#### get expression ####
run <- list.files(exp_dir)
exp_dfs <-
    map(file.path(exp_dir, run, "miRNAs_expressed_all_samples_now.csv"),
        read_delim, delim = "\t", col_types = cols(`#miRNA` = col_character()))
exp_df <- tibble(run = run, df = exp_dfs)
exp_df <- exp_df %>% left_join(metadata_df, by = "run")

# sum total read count for quality check
exp_df <- exp_df %>% mutate(total_rc = map_dbl(df, ~ sum(.x$read_count)))

# filter libraries based on total read count
# exp_df <- exp_df %>% filter(total_rc > 1e5)
sample_num <- exp_df %>% count(species, tissue) %>% spread(tissue, n)

exp_df <- exp_df %>%
    group_by(species, tissue) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(id = str_c(species, tissue, sample_num, sep = "-"))

# clean up dfs for joining
exp_df <- exp_df %>%
    # add col to indicate 5p or 3p mature
    mutate(df = map(df, ~ .x %>%
                        mutate(arm = case_when(
                            str_detect(`#miRNA`, "5p") ~ "5p",
                            str_detect(`#miRNA`, "3p") ~ "3p",
                            TRUE ~ "undefined"
                        )) %>%
                        rename(rpm = `seq(norm)`) %>%
                        select(precursor, arm, rpm))) %>%
    # sum 5p and 3p to get total expression
    mutate(precursor_df = map(df, ~ .x %>% group_by(precursor) %>% summarise(rpm = sum(rpm))))

rename_cols <- function(df, first_col, last_col) {
    colnames(df)[1] <- first_col
    colnames(df)[ncol(df)] <- last_col
    return(df)
}

exp_df <- exp_df %>%
    mutate(df = pmap(list(df, species, id), rename_cols),
           precursor_df = pmap(list(precursor_df, species, id), rename_cols))

#### precursor expression ####
precursor_exp <-
    reduce(c(list(class_df), exp_df$precursor_df), left_join) %>%
    rename(Evolutionary_mode = "Evolutionary mode")

# calculate miRNA expression ranking in D. mel and D. sim
precursor_percentile <- precursor_exp %>%
    select(miRNA, Age, Evolutionary_mode, dme, dsi,
           contains("dme"), contains("dsi"), -contains("male_body")) %>%
    drop_na()

precursor_percentile_num <-
    precursor_percentile %>% count(Evolutionary_mode)

precursor_percentile <- precursor_percentile %>%
    mutate_if(is.numeric, percent_rank) %>%
    mutate_if(is.numeric, ~ .x >= 0.5) %>%
    mutate(top50_num = rowSums(.[grep("testes", names(.))]))

precursor_percentile_top50_num <- precursor_percentile %>%
    filter(top50_num == 8) %>%
    count(Evolutionary_mode)

make_plot_df <- function(df) {
    plot_df <- df %>%
        select(-c(dme, dsi, dse, der, dvi)) %>%
        gather(key = "sample", value = "rpm", -c(miRNA, Age, Evolutionary_mode)) %>%
        separate(sample, into = c("species", "tissue", "sample_num"), sep = "-") %>%
        filter(rpm > 0) %>%
        mutate(log2_rpm = log2(rpm)) %>%
        mutate(species = fct_relevel(species, "dme", "dsi", "dse", "der", "dvi")) %>%
        mutate(Category = case_when(
            Age == "0 - 4 Myrs" & Evolutionary_mode == "adaptive" ~ "Adaptive - < 4 Myrs",
            Age == "4 - 30 Myrs" & Evolutionary_mode == "adaptive" ~ "Adaptive - 4 to 30 Myrs",
            Age %in% c("30 - 60 Myrs", "60 - 250 Myrs", "> 250 Myrs") &
                Evolutionary_mode == "adaptive" ~ "Adaptive - > 30 Myrs",
            Age %in% c("30 - 60 Myrs", "60 - 250 Myrs", "> 250 Myrs") &
                Evolutionary_mode == "conservative" ~ "Conservative - > 30 Myrs",
            TRUE ~ "NA"
        )) %>%
        mutate(Category =
                   fct_relevel(Category, "Adaptive - < 4 Myrs", "Adaptive - 4 to 30 Myrs",
                               "Adaptive - > 30 Myrs", "Conservative - > 30 Myrs")) %>%
        mutate(microRNA = case_when(
            str_detect(miRNA, "bantam") ~ "bantam",
            str_detect(miRNA, "mir-184") ~ "miR-184",
            str_detect(miRNA, "mir-983") ~ "miR-983",
            str_detect(miRNA, "mir-975") ~ "miR-975"
        ))
    return(plot_df)
}

precursor_plot <- make_plot_df(precursor_exp)

precursor_plot <-
    precursor_plot %>% mutate(tissue = str_replace(tissue, "_", " "))

#### global expression for all miRNAs surveyed in Lyu et al. ####
p <- precursor_plot %>%
    filter(Evolutionary_mode != "transitional") %>%
    ggplot() +
    geom_boxplot(mapping = aes(x = sample_num, y = log2_rpm, color = Category)) +
    # geom_jitter(mapping = aes(x = sample_num, y = log2_rpm, color = Evolutionary_mode)) +
    facet_grid(tissue ~ species, scales = "free", space = "free", switch = "y") + # use  to keep same width
    labs(x = "Samples in different species",
         y = expression(log[2]*"(Reads Per Million)")) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_blank())

ggsave(file.path(fig_out, "all_miR_expressions.png"), plot = p, width = 10, height = 6)

#### plot bantam, miR-184, miR-975, miR-983 ####
precursor_subplot <- precursor_plot %>% filter(!is.na(microRNA) & miRNA != "dme-mir-983-2")

p <- precursor_plot %>%
    filter(Evolutionary_mode != "transitional") %>%
    ggplot() +
    geom_boxplot(mapping = aes(x = sample_num, y = log2_rpm)) +
    geom_jitter(mapping = aes(x = sample_num, y = log2_rpm, color = microRNA), data = precursor_subplot, size = 3) +
    facet_grid(tissue ~ species, scales = "free", space = "free", switch = "y") + # use  to keep same width
    labs(x = "Samples in different species",
         y = expression(log[2]*"(Reads Per Million)")) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_blank())

ggsave(file.path(fig_out, "miR983_975_expressions.png"), plot = p, width = 10, height = 4)

#### mature expression ####
mature_exp <-
    reduce(c(list(class_df), exp_df$df), left_join) %>%
    rename(Evolutionary_mode = "Evolutionary mode") %>%
    mutate(miRNA = str_replace(str_sub(miRNA, start = 5), "mir", "miR")) %>%
    select(-c(dme, dsi, dse, der, dvi)) %>%
    # rename miR-983 for visualization
    mutate(miRNA = str_replace(miRNA, "miR-983-1", "miR-983"))

mature_prop <- mature_exp %>%
    group_by(miRNA) %>%
    mutate(across(where(is.numeric), ~ . / sum(.))) %>%
    ungroup()

# analysis for co-expression
major_prop <-
    mature_prop %>%
    group_by(miRNA) %>%
    summarise_if(is.numeric, max)

major_prop <-
    major_prop %>%
    # only analyze testes samples
    select(miRNA, contains("testes")) %>%
    mutate_if(is.numeric, ~ .x < 0.75) %>%
    mutate(lib_n = rowSums(.[, 2:ncol(.)], na.rm = TRUE))

coexp_miRNA <- major_prop %>%
    filter(lib_n >= 2) %>%
    select(miRNA) %>%
    mutate(coexpression = "Y")

# analysis for arm switching
arm_idx <- mature_prop %>%
    group_by(miRNA) %>%
    mutate_if(is.numeric, ~ .x == max(.x)) %>%
    ungroup() %>%
    filter(arm == "5p")

arm_idx <- arm_idx %>%
    select(miRNA, contains("testes")) %>%
    mutate(n_5p = rowSums(.[, 2:ncol(.)], na.rm = TRUE),
           n_3p = rowSums(!.[, 2:ncol(.)], na.rm = TRUE))

armswitch_miRNA <- arm_idx %>%
    filter(n_5p >= 2 & n_3p >= 2) %>%
    select(miRNA) %>%
    mutate(arm_switching = "Y")

# generate Table 1 for miRNAs with co-expression or arm switching
table_1 <-
    mature_exp %>%
        select(miRNA, Evolutionary_mode) %>%
        unique() %>%
        left_join(coexp_miRNA, by = "miRNA") %>%
        left_join(armswitch_miRNA, by = "miRNA") %>%
        filter(Evolutionary_mode %in% c("adaptive", "conservative")) %>%
        filter(!is.na(coexpression) | !is.na(arm_switching))

# adaptive_miR <- mature_prop %>%
#     filter(Evolutionary_mode == "adaptive", Age != "0 - 4 Myrs") %>%
#     select(-contains("dvi"), -contains("male_body"))

adaptive_miR_young <- mature_prop %>%
    filter(Evolutionary_mode == "adaptive", Age %in% c("0 - 4 Myrs", "4 - 30 Myrs")) %>%
    select(-contains("dvi"), -contains("male_body"))

adaptive_miR_old <- mature_prop %>%
    filter(Evolutionary_mode == "adaptive", Age %in% c("30 - 60 Myrs", "60 - 250 Myrs")) %>%
    select(-contains("male_body"))

conserved_miR <- mature_prop %>%
    filter(Evolutionary_mode == "conservative") %>%
    select(-contains("male_body"))

selected_miR <- mature_prop %>%
    filter(miRNA %in% c("bantam", "miR-184", "miR-983", "miR-975")) %>%
    select(-contains("male_body"))

make_plot_df_2 <- function(df) {
    df %>%
        gather(key = "sample", value = "rpm", -c(miRNA, Age, Evolutionary_mode, arm)) %>%
        separate(sample, into = c("species", "tissue", "sample_num"), sep = "-") %>%
        filter(!is.na(rpm)) %>%
        mutate(species = fct_relevel(species, "dme", "dsi", "dse", "der", "dvi"))
}

# adaptive_miR_plot <- make_plot_df_2(adaptive_miR)
adaptive_miR_young_plot <- make_plot_df_2(adaptive_miR_young)
adaptive_miR_old_plot <- make_plot_df_2(adaptive_miR_old)
conserved_miR_plot <- make_plot_df_2(conserved_miR)
selected_miR_plot <- make_plot_df_2(selected_miR)

# # adaptive miRNAs which are annotated in miRBase and detected in more than 3 species
# adaptive_miR_sel <- adaptive_miR_plot %>%
#     select(miRNA, species) %>% unique() %>% count(miRNA) %>%
#     filter(n >= 3) %>% pull(miRNA)
#
# p <- adaptive_miR_plot %>%
#     filter(miRNA %in% adaptive_miR_sel) %>% ggplot() +
#     geom_point(mapping = aes(x = species, y = rpm, color = arm)) +
#     geom_hline(yintercept = c(0.75, 0.25), linetype="dashed", color = "grey", size = 1.2) +
#     facet_wrap("miRNA", nrow = 2) +
#     labs(x = "Drosophila species",
#          y = "Proportion") +
#     scale_color_brewer(palette = "Set2") +
#     cowplot::theme_cowplot()
#
# ggsave(file.path(fig_out, "adaptive_miRNA_5p_3p_proportion.png"), plot = p, width = 10, height = 5)

p <- selected_miR_plot %>%
    ggplot() +
    geom_point(mapping = aes(x = species, y = rpm, color = arm)) +
    geom_hline(yintercept = c(0.75, 0.25), linetype="dashed", color = "grey", size = 1.2) +
    facet_wrap("miRNA", nrow = 1) +
    labs(x = "Drosophila species",
         y = "Proportion") +
    scale_color_brewer(palette = "Set2") +
    cowplot::theme_cowplot()

ggsave(file.path(fig_out, "selected_miRNA_5p_3p_proportion.png"), plot = p, width = 8, height = 2.5)

adaptive_miR_young_sel <- adaptive_miR_young_plot %>%
    select(miRNA, species) %>% unique() %>% count(miRNA) %>%
    filter(n >= 3) %>% pull(miRNA)

p <- adaptive_miR_young_plot %>%
    filter(miRNA %in% adaptive_miR_young_sel) %>% ggplot() +
    geom_point(mapping = aes(x = species, y = rpm, color = arm)) +
    geom_hline(yintercept = c(0.75, 0.25), linetype="dashed", color = "grey", size = 1.2) +
    facet_wrap("miRNA", nrow = 2) +
    labs(x = "Drosophila species",
         y = "Proportion") +
    scale_color_brewer(palette = "Set2") +
    cowplot::theme_cowplot()

ggsave(file.path(fig_out, "adaptive_miRNA_young_5p_3p_proportion.png"), plot = p, width = 6, height = 5)

adaptive_miR_old_sel <- adaptive_miR_old_plot %>%
    select(miRNA, species) %>% unique() %>% count(miRNA) %>%
    filter(n >= 3) %>% pull(miRNA)

p <- adaptive_miR_old_plot %>%
    filter(miRNA %in% adaptive_miR_old_sel) %>% ggplot() +
    geom_point(mapping = aes(x = species, y = rpm, color = arm)) +
    geom_hline(yintercept = c(0.75, 0.25), linetype="dashed", color = "grey", size = 1.2) +
    facet_wrap("miRNA", nrow = 2) +
    labs(x = "Drosophila species",
         y = "Proportion") +
    scale_color_brewer(palette = "Set2") +
    cowplot::theme_cowplot()

ggsave(file.path(fig_out, "adaptive_miRNA_old_5p_3p_proportion.png"), plot = p, width = 6, height = 5)

# conserved miRNAs which are annotated in miRBase and detected in all 5 species
conserved_miR_sel_1 <-
    conserved_miR_plot %>% select(miRNA, species) %>% unique() %>% count(miRNA) %>%
    filter(n == 5) %>% pull(miRNA)

p <- conserved_miR_plot %>%
    #filter(miRNA %in% c("bantam", "miR-184", "miR-281-1", "miR-310",
    #                    "miR-iab-4")) %>%
    filter(miRNA %in% conserved_miR_sel_1) %>%
    ggplot() +
    geom_point(mapping = aes(x = species, y = rpm, color = arm)) +
    geom_hline(yintercept = c(0.75, 0.25), linetype="dashed", color = "grey", size = 1.2) +
    facet_wrap("miRNA", nrow = 1) +
    labs(x = "Drosophila species",
         y = "Proportion") +
    scale_color_brewer(palette = "Set2") +
    cowplot::theme_cowplot()

ggsave(file.path(fig_out, "conserved_miRNA_5p_3p_proportion_5_species.png"),
       plot = p, width = 10, height = 3)

conserved_miR_sel_2 <-
    conserved_miR_plot %>% select(miRNA, species) %>%
    unique() %>% filter(species %in% c("dme", "dvi")) %>% count(miRNA) %>%
    filter(!(miRNA %in% conserved_miR_sel_1), n == 2) %>% pull(miRNA)

# conserved miRNAs which are annotated in miRBase and detected in dme and dvi
p <- conserved_miR_plot %>%
    filter(miRNA %in% conserved_miR_sel_2, species %in% c("dme", "dvi")) %>%
    ggplot() +
    geom_point(mapping = aes(x = species, y = rpm, color = arm)) +
    geom_hline(yintercept = c(0.75, 0.25), linetype="dashed", color = "grey", size = 1.2) +
    facet_wrap("miRNA", nrow = 6) +
    labs(x = "Drosophila species",
         y = "Proportion") +
    scale_color_brewer(palette = "Set2") +
    cowplot::theme_cowplot()

ggsave(file.path(fig_out, "conserved_miRNA_5p_3p_proportion_2_species.png"),
       plot = p, width = 20, height = 12)
