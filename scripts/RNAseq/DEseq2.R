#library(data.table)

library(tximport)
library(DESeq2)

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

root_dir <- "/home/yixin/NutstoreFiles/Nutstore/Wu_lab"
# root_dir <- "/Users/yizhao/我的坚果云/Wu_lab"
# root_dir <- "D:/NutstoreFiles/Wu_lab"

output_dir <- file.path(root_dir, 'Snakemake_projects/miR983_975/output/deseq2')
exp_dir <- file.path(root_dir, 'Snakemake_projects/miR983_975/output/expression')

make_DEseq2result <- function(rds_name){
    ddsTxi <- readRDS(file.path(exp_dir, rds_name))
    
    # filter low total count genes
    keep <- rowSums(counts(ddsTxi)) >= 1
    ddsTxi <- ddsTxi[keep, ]
    
    # relevel genotype, use control as ref level
    ddsTxi$genotype <- relevel(ddsTxi$genotype, ref = "control")
    
    # Run DEseq2
    ddsTxi <- DESeq(ddsTxi)
    
    resLFC_983 <- lfcShrink(ddsTxi, coef="genotype_miR983KO_vs_control")
    resLFC_975 <- lfcShrink(ddsTxi, coef="genotype_miR975KO_vs_control")
    
    list(resLFC_983, resLFC_975, ddsTxi)
    }

dme_results <- make_DEseq2result("dmel_Txi.rds")
dsi_results <- make_DEseq2result("dsim_Txi.rds")

write.csv(dme_results[[1]], file.path(output_dir, "dme_miR983_DEseq2.csv"))
write.csv(dme_results[[2]], file.path(output_dir, "dme_miR975_DEseq2.csv"))

write.csv(dsi_results[[1]], file.path(output_dir, "dsi_miR983_DEseq2.csv"))
write.csv(dsi_results[[2]], file.path(output_dir, "dsi_miR975_DEseq2.csv"))

save_fig <- function(file_name, plt_cmd){
    png(filename = file.path(output_dir, file_name))
    plt_cmd
    dev.off()
    }

par(mar=c(1, 1, 1, 1))

save_fig("dme-miR-983KO_mis-regulated_genes.png",
         plotMA(dme_results[[1]], main = "", xlab = "", ylab = "", ylim = c(-2, 2), alpha = 0.05, cex.main=1.8, cex.lab=1.8, cex.axis = 1.8))
save_fig("dme-miR-975KO_mis-regulated_genes.png",
         plotMA(dme_results[[2]], main = "", xlab = "", ylab = "", ylim = c(-2, 2), alpha = 0.05, cex.main=1.8, cex.lab=1.8, cex.axis = 1.8))
save_fig("dsi-miR-983KO_mis-regulated_genes.png",
         plotMA(dsi_results[[1]], main = "", xlab = "", ylab = "", ylim = c(-2, 2), alpha = 0.05, cex.main=1.8, cex.lab=1.8, cex.axis = 1.8))
save_fig("dsi-miR-975KO_mis-regulated_genes.png",
         plotMA(dsi_results[[2]], main = "", xlab = "", ylab = "", ylim = c(-2, 2), alpha = 0.05, cex.main=1.8, cex.lab=1.8, cex.axis = 1.8))


# # get DEseq results
# res_983 <- results(ddsTxi, contrast=c("genotype","miR983KO","control"), alpha = 0.05)
# res_975 <- results(ddsTxi, contrast=c("genotype","miR975KO","control"), alpha = 0.05)
# 
# mcols(res_983)$description
# 
# summary(res_983)
# summary(res_975)
# 
# # Log fold change shrinkage for visualization and ranking
# # get coef names
# resultsNames(ddsTxi)
# 
# summary(resLFC_983)
# summary(resLFC_975)
# 
# # p-values and adjusted p-values
# sum(res_983$padj < 0.05, na.rm=TRUE)
# sum(res_975$padj < 0.05, na.rm=TRUE)
# 
# # MA plot 
# # plotMA(res_983)
# 
# plotMA(resLFC_983, ylim=c(-2,2))
# plotMA(resLFC_975, ylim=c(-2,2))
# 
# # Plot counts
# plotCounts(ddsTxi, gene=which.min(res_983$padj), intgroup="genotype")
# plotCounts(ddsTxi, gene=which.min(res_975$padj), intgroup="genotype")

