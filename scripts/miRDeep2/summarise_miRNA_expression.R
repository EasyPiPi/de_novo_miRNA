library(tidyverse)

exp_dir <- '/home/yixin/Desktop/github_repo/de_novo_miRNA/outputs/miRDeep2/quantifier'
metadata_df <- read_csv('/home/yixin/Desktop/github_repo/de_novo_miRNA/metadata/small_RNAseq/metadata.csv',
                        col_types = cols(run = col_character()))
orth_df <- readxl::read_xlsx('/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/miRNA_ortholog/Supplemental_Table_S4.xlsx')
class_df <- read_csv('/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/miRNA_classification/miRNA_category.csv')

#### get ortholog names ####
orth_df <- orth_df %>% select(familyname, OrthologName)

class_df <- class_df %>%
    left_join(orth_df, by = c("miRNA" = "familyname")) %>%
    mutate(species = str_sub(OrthologName,1,3)) %>%
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
exp_df <- exp_df %>% filter(total_rc > 1e5)
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

precursor_exp <-
    reduce(c(list(class_df), exp_df$precursor_df), left_join) %>%
    rename(Evolutionary_mode = "Evolutionary mode")

precursor_plot <- precursor_exp %>%
    select(-c(dme, dsi, dse, der, dvi)) %>%
    gather(key = "sample", value = "rpm", -c(miRNA, Age, Evolutionary_mode)) %>%
    separate(sample, into = c("species", "tissue", "sample_num"), sep = "-") %>%
    filter(rpm > 0) %>%
    mutate(log2_rpm = log2(rpm)) %>%
    mutate(species = fct_relevel(species, "dme", "dsi", "dse", "der", "dvi")) %>%
    mutate(Category = case_when(
        Age == "0 - 4 Myrs" & Evolutionary_mode == "adaptive" ~ "Adaptive - < 4 Myrs",
        Age == "4 - 30 Myrs" & Evolutionary_mode == "adaptive" ~ "Adaptive - 4 to 30 Myrs",
        Age %in% c("30 - 60 Myrs", "60 - 250 Myrs", "> 250 Myrs") & Evolutionary_mode == "adaptive" ~ "Adaptive - > 30 Myrs",
        Age %in% c("30 - 60 Myrs", "60 - 250 Myrs", "> 250 Myrs") & Evolutionary_mode == "conservative" ~ "Conservative - > 30 Myrs",
        TRUE ~ "NA"
    )) %>%
    mutate(Category = fct_relevel(Category, "Adaptive - < 4 Myrs", "Adaptive - 4 to 30 Myrs",
                                  "Adaptive - > 30 Myrs", "Conservative - > 30 Myrs")) %>%
    mutate(microRNA = case_when(
        str_detect(miRNA, "bantam") ~ "bantam",
        str_detect(miRNA, "mir-184") ~ "miR-184",
        str_detect(miRNA, "mir-983") ~ "miR-983",
        str_detect(miRNA, "mir-975") ~ "miR-975"
    ))

#### global expression for all miRNAs surveyed in Lyu et al. ####
precursor_plot %>%
    filter(Evolutionary_mode != "transitional") %>%
    ggplot() +
    geom_boxplot(mapping = aes(x = sample_num, y = log2_rpm, color = Category)) +
    # geom_jitter(mapping = aes(x = sample_num, y = log2_rpm, color = Evolutionary_mode)) +
    facet_grid(tissue ~ species, scales = "free", space = "free", switch = "y") + # use  to keep same width
    labs(x = "Samples in different species",
         y = expression(log[2]*"(Reads Per Million)")) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_blank())

#### plot bantam, miR-184, miR-975, miR-983 ####
precursor_subplot <- precursor_plot %>% filter(!is.na(microRNA) & miRNA != "dme-mir-983-2")

precursor_plot %>%
    filter(Evolutionary_mode != "transitional") %>%
    ggplot() +
    geom_boxplot(mapping = aes(x = sample_num, y = log2_rpm)) +
    geom_jitter(mapping = aes(x = sample_num, y = log2_rpm, color = microRNA), data = precursor_subplot) +
    facet_grid(tissue ~ species, scales = "free", space = "free", switch = "y") + # use  to keep same width
    labs(x = "Samples in different species",
         y = expression(log[2]*"(Reads Per Million)")) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_blank())
