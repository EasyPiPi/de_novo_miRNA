# build significant and up/down-regulated categories for fisher exact test
get_category <- function(df){
    df %>% 
        mutate(log2FoldChange_cat = ifelse(log2FoldChange > 0, "Up-regulaed", "Down-regulated"),
               significant = ifelse(padj > 0.05 | is.na(padj), "Other", "Misregulated"))
}

# get fisher exact test p values
get_fisher_pvalue <- function(df) {
    df %>% 
        select(c(2:3)) %>% 
        fisher.test() %>% broom::tidy() %>%
        select(p.value) %>% pull()
}

# get contingency table for statistical test
# chromosome location
get_chr_loc_tab <- function(df, chr_loc) {
    df %>%
        left_join(chr_loc, by = c("gene_ID" = "ensembl_gene_id")) %>%
        mutate(chromosome = ifelse(chromosome_name == "X", "X", "Autosome")) %>%
        count(significant, chromosome) %>%
        spread(significant, n) 
}

# sex bias genes
get_sexbias_tab <- function(df, sebida){
    df %>% 
        left_join(sebida, by = "gene_ID") %>%
        count(MetaClass, significant) %>%
        filter(MetaClass %in% c("Female", "Male")) %>%
        spread(MetaClass, n)
} 

get_sexbias_updown_tab <- function(df, sebida){
    df %>% 
        left_join(sebida, by = "gene_ID") %>%
        filter(significant == "Misregulated") %>%
        count(MetaClass, log2FoldChange_cat) %>%
        filter(MetaClass %in% c("Female", "Male")) %>%
        spread(MetaClass, n)
} 

# testes specific expressed genes
get_testes_tab <- function(df) {
    df %>% 
        mutate(testes_exp = ifelse(gene_ID %in% testes_exp_genes, "Testes-specific", "Other")) %>% 
        count(significant, testes_exp) %>% 
        spread(testes_exp, n)
}

####################################################################################
# Visualization
get_chr_loc_plt <- function(df, title, rm_legend = TRUE) {
    p <- df %>%
        gather(gene, n, -chromosome) %>%
        ggplot(aes(y = n, fill = chromosome, x = gene)) +
        geom_col(position="fill") +
        geom_text(aes(label = n), position = position_fill(vjust=0.5)) + 
        labs(y = "", fill = "Chromosome", x = "",
             title = title) +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_brewer(palette = "Pastel1") + 
        theme_cowplot() +
        theme(axis.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 9))
    if (rm_legend) p <- p + theme(legend.position="None")
    # ggsave(file.path(output_dir, file_name), width = 5, height = 5)
    if (!rm_legend) p <- p + theme(legend.position="right")
    p
}

get_sexbias_plt <- function(df, rm_legend = TRUE) {
    p <- df %>% 
        gather(Sex_biased, n, -significant) %>%
        ggplot(aes(y = n, x = significant, fill = Sex_biased)) +
        geom_col(position="fill") +
        geom_text(aes(label = n), position = position_fill(vjust=0.5)) + 
        labs(y = "", x = "", fill = "Sex biased") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_brewer(palette = "Pastel2", direction = -1) + 
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (rm_legend) p <- p + theme(legend.position="None")
    if (!rm_legend) p <- p + theme(legend.position="right")
    p
    # ggsave(file.path(output_dir, file_name), width = 5, height = 5)
}

get_sexbias_updown_plt <- function(df, file_name) {
    df %>% 
        gather(Sex_biased, n, -log2FoldChange_cat) %>%
        ggplot(aes(y = n, x = Sex_biased, fill = log2FoldChange_cat)) +
        geom_col(position="fill") +
        geom_text(aes(label = n), position = position_fill(vjust=0.5)) + 
        labs(y = "Proportion", x = "Sex biased", fill = "Misregulated gene") +
        scale_y_continuous(labels = scales::percent)
    ggsave(file.path(output_dir, file_name), width = 5, height = 5)
}

get_testes_plt <- function(df, file_name) {
    df %>% 
        gather(testes_specific, n, -significant) %>%
        ggplot(aes(y = n, fill = testes_specific, x = significant)) +
        geom_col(position="fill") +
        geom_text(aes(label = n), position = position_fill(vjust=0.5)) + 
        labs(y = "Proportion", fill = "Testes-sepcific", x = "Gene") +
        scale_y_continuous(labels = scales::percent)
    ggsave(file.path(output_dir, file_name), width = 5, height = 5)
}