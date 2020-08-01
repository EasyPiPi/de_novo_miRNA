library(tidyverse)

exp_dir <- '/home/yixin/Desktop/github_repo/de_novo_miRNA/outputs/miRDeep2/quantifier'
metadata_df <- read_csv('/home/yixin/Desktop/github_repo/de_novo_miRNA/metadata/small_RNAseq/metadata.csv',
                        col_types = cols(run = col_character()))

run <- list.files(exp_dir)

exp_dfs <-
    map(file.path(exp_dir, run, "miRNAs_expressed_all_samples_now.csv"),
        read_delim, delim = "\t", col_types = cols(`#miRNA` = col_character()))
exp_df <- tibble(run = run, df = exp_dfs)
exp_df <- exp_df %>% left_join(metadata_df, by = "run")

# sum total read count for quality check
exp_df <- exp_df %>% mutate(total_rc = map_dbl(df, ~ sum(.x$read_count)))

# filter libraries based on total read count
exp_df <- exp_df %>% filter(total_rc > 0)
sample_num <- exp_df %>% count(species, tissue) %>% spread(tissue, n)



exp_df %>% group_by(species, tissue)

exp_df$df[[1]]
