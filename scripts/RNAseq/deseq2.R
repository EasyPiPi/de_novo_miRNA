suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(tximport)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(apeglm)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GGally)))

# root_dir <- '/home/yixin/Desktop/github_repo/de_novo_miRNA'
# setDTthreads(threads = 4)
setDTthreads(threads = snakemake@threads[[1]])
root_dir <- snakemake@params[["root_dir"]]
external_dir <- file.path(root_dir, "external_resources")
salmon_dir <- file.path(root_dir, 'outputs/salmon')
output_dir <- file.path(root_dir, 'outputs/expression')
table_out <- file.path(output_dir, "table")
figure_out <- file.path(output_dir, "figure")

dir.create(table_out, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_out, showWarnings = FALSE, recursive = TRUE)

# generate tx2gene from cdna files for later use in deseq2
make_tx2gene <- function(species) {
    tx2gene_raw <-
        fread(cmd = paste(
            "grep",
            "'>'",
            file.path(external_dir, species, "cdna.fasta")
        ), sep = ";")

    tx2gene <- tx2gene_raw[, c("V3", "V8")]
    colnames(tx2gene) <- c("transcript_id", "gene_id")
    tx2gene[, transcript_id := str_remove(transcript_id, "ID=")]
    tx2gene[, gene_id := str_remove(gene_id, "parent=")]
    return(tx2gene)
}

tx2gene_dme <- make_tx2gene("dme")
tx2gene_dsi <- make_tx2gene("dsi")

# read in salmon files and make ddsTxi (gene-level analysis)
make_ddsTxi <- function(species, tx2gene, sample_rm){
    # prepare file diectory for reading expression files
    salmon_subdirs <- list.files(file.path(salmon_dir, species))
    salmon_subdirs <- salmon_subdirs[!salmon_subdirs %in% sample_rm]
    # RNA-seq expression
    salmon_files <- file.path(salmon_dir, species, salmon_subdirs, "quant.sf")
    # read in rnaseq files
    names(salmon_files) <- salmon_subdirs
    # tx_exp <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene, txOut = T)
    gn_exp <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene)
    # generate sample data table
    samples <- data.table(id = salmon_subdirs)
    samples[, species := fifelse(str_detect(id, "DM"), "dme", "dsi")]
    samples[, genotype := as.factor(fcase(
        str_detect(id, "975"), "miR975KO",
        str_detect(id, "983"), "miR983KO",
        str_detect(id, "CT"), "wildtype"
    ))]
    samples[, replicate := str_extract(id, "[1-9]$")]
    samples <- as.data.frame(samples)
    rownames(samples) <- samples$id

    # read in files
    ddsTxi <- DESeqDataSetFromTximport(gn_exp,
                                       colData = samples,
                                       design = ~ genotype)
    return(ddsTxi)
}

# samples with abnormal GC content are removed
# remove dm975 KO 2 from further analysis
ddsTxi_dme <-
    make_ddsTxi(species = "dme", tx2gene = tx2gene_dme,
                sample_rm = c("DM975KO2"))
# remove ds975 KO 5 and ds983 KO 4 from further analysis
ddsTxi_dsi <-
    make_ddsTxi(species = "dsi", tx2gene = tx2gene_dsi,
                sample_rm = sample_rm <- c("DS975KO5", "DS983KO4"))

#### data visualization for replicates ####
# Pairwise scatter plot
# Revise the lower plot, change its transparency
my_point <- function(data, mapping, ...) {
    ggplot(data = data, mapping=mapping) +
        geom_point(..., alpha=0.01)}

make_scatter <- function(ddsTxi){
    vsd <- vst(ddsTxi, blind=F)
    df <- as.data.table(assay(vsd))
    ggpairs(df, mapping = aes(), lower = list(continuous = my_point),
                upper = list(continuous = wrap(ggally_cor, size=3,
                                               hjust=0.15, alignPercent=0.9)))
}

# Data transformations and visualization
p <- make_scatter(ddsTxi = ddsTxi_dme)
ggsave(file.path(figure_out, "dme_replicate_scatter.png"), plot = p, width = 10, height = 8)

p <- make_scatter(ddsTxi = ddsTxi_dsi)
ggsave(file.path(figure_out, "dsi_replicate_scatter.png"), plot = p, width = 10, height = 8)

#### Run DEseq2 ####
make_DEseq2result <- function(ddsTxi){
    # filter low total count genes
    keep <- rowSums(counts(ddsTxi)) >= 1
    ddsTxi <- ddsTxi[keep, ]

    # relevel genotype, use control as ref level
    ddsTxi$genotype <- relevel(ddsTxi$genotype, ref = "wildtype")

    # Run DEseq2
    ddsTxi <- DESeq(ddsTxi)

    resLFC_983 <- lfcShrink(ddsTxi, coef="genotype_miR983KO_vs_wildtype", type = "apeglm")
    resLFC_975 <- lfcShrink(ddsTxi, coef="genotype_miR975KO_vs_wildtype", type = "apeglm")

    list(resLFC_983, resLFC_975, ddsTxi)
}

dme_results <- make_DEseq2result(ddsTxi = ddsTxi_dme)
dsi_results <- make_DEseq2result(ddsTxi = ddsTxi_dsi)

write.csv(dme_results[[1]], file.path(table_out, "dme_miR983_DEseq2.csv"))
write.csv(dme_results[[2]], file.path(table_out, "dme_miR975_DEseq2.csv"))
write.csv(dsi_results[[1]], file.path(table_out, "dsi_miR983_DEseq2.csv"))
write.csv(dsi_results[[2]], file.path(table_out, "dsi_miR975_DEseq2.csv"))

#### data visualization for mis-regulated genes ####
save_fig <- function(file_name, plt_cmd){
    png(filename = file.path(figure_out, file_name),
        units = "in", width = 4, height = 4, res = 300)
    p <- plt_cmd
    invisible(dev.off())
}

# par(c(5, 4, 4, 2) + 0.1)
plotMA <- function(df, main) {
    DESeq2::plotMA(df, ylab = expression(log[2]*"(fold change)"),
                   main = main,
                   ylim = c(-2, 2), alpha = 0.05)
}

save_fig("dme-miR-983KO_mis-regulated_genes.png", plotMA(dme_results[[1]], main = "dme-miR-983 KO"))
save_fig("dme-miR-975KO_mis-regulated_genes.png", plotMA(dme_results[[2]], main = "dme-miR-975 KO"))
save_fig("dsi-miR-983KO_mis-regulated_genes.png", plotMA(dsi_results[[1]], main = "dsi-miR-983 KO"))
save_fig("dsi-miR-975KO_mis-regulated_genes.png", plotMA(dsi_results[[2]], main = "dsi-miR-975 KO"))

