library(data.table)
library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(GGally)

root_dir <- '/home/yixin/NutstoreFiles/Nutstore/Wu_lab'
ortholog_dir <- file.path(root_dir, "Snakemake_projects/miR983_975/output/resources/ortholog")
output_dir <- file.path(root_dir, 'Snakemake_projects/miR983_975/output/expression')

# read in transcript info
# dmel_tx2gene <- fread("dmel_tx2gene.csv", drop = "V1", header = T)
# dsim_tx2gene <- fread("dsim_tx2gene.csv")
# dsim_tx2gene <- dsim_tx2gene[,c(2,1)]

tx2gene <- fread(file.path(output_dir, "fbgn_fbtr_fbpp_fb_2019_04_hearder_removed.tsv"))
tx2gene <- tx2gene[, c(2,1)]

# read in orthology table
dsim2dmel_1to1_ortholog <- fread(file.path(ortholog_dir, "dsim2dmel_1to1_ortholog.csv"))
# only analyze 1 to 1 orthologs
tx2gene <- tx2gene[tx2gene$V1 %in% c(dsim2dmel_1to1_ortholog$ensembl_gene_id, dsim2dmel_1to1_ortholog$dmelanogaster_eg_homolog_ensembl_gene)]

# generate tx2gene table for dsim expression quantification, using dmel gene names
dsim_tx2gene <- merge.data.frame(tx2gene, dsim2dmel_1to1_ortholog[,c("ensembl_gene_id", "dmelanogaster_eg_homolog_ensembl_gene")], by.x = "V1", by.y = "ensembl_gene_id")
dsim_tx2gene <- dsim_tx2gene[, c("V2", "dmelanogaster_eg_homolog_ensembl_gene")]


library(tximport)
library(DESeq2)

################################################################################
# D. melanogaster
# prepare file diectory for reading expression files
dir <- '/home/yixin/Desktop/github_repo/de_novo_miRNA/outputs/salmon/dme'
dmfiles <- grep("DM", list.files(dir), value = T)
# remove dm975 KO 2 from further analysis
# dmfiles <- dmfiles[dmfiles != "DM975KO2"]
# RNA-seq expression
files <- file.path(dir, dmfiles, "quant.sf")
# read in rnaseq files
names(files) <- dmfiles
tx_exp <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = T)
write.csv(tx_exp$abundance, file.path(output_dir, "dmel_tx_TPM.csv"))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# generate sample data table
samples <- data.frame(row.names = dmfiles)
samples$genotype <- as.factor(c(rep("miR975KO", 2), rep("miR983KO", 3), rep("control", 3)))
samples$species <- rep("dmel", 8)
samples$replicate <- as.factor(c(1, 3, rep(1:3,2)))
# read in files
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ genotype)

# Data transformations and visualization
vsd <- vst(ddsTxi, blind=F)

# save files
saveRDS(ddsTxi, file.path(output_dir, "dmel_Txi.rds"))
write.csv(assay(vsd), file.path(output_dir, "dmel_vst.csv"))

################################################################################
# D. simulans
# prepare file diectory for reading expression files
dir <- '/home/yixin/Desktop/project_data/miR983_975/all/salmon/dsim_index'
# dsfiles_bac <- grep("^DS+.*[1|2|3]$",  list.files(dir), value = T)
dsfiles_nonbac <- grep("^DS+.*[4|5|6]$",  list.files(dir), value = T)

# filter two samples
dsfiles_nonbac <- dsfiles_nonbac[! dsfiles_nonbac %in% c("DS975KO5", "DS983KO4")]

files <- file.path(dir, dsfiles_nonbac, "quant.sf")
# read in rnaseq files
names(files) <- dsfiles_nonbac
tx_exp <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = T)
write.csv(tx_exp$abundance, file.path(output_dir, "dsim_tx_TPM.csv"))

txi <- tximport(files, type = "salmon", tx2gene = dsim_tx2gene)

samples <- data.frame(row.names = dsfiles_nonbac)

samples$genotype <- as.factor(c(rep("miR975KO", 2), rep("miR983KO", 2), rep("control", 3)))
samples$species <- rep("dsim", 7)
samples$replicate <- as.factor(c(4, 6, 5, 6, 4:6))

# read in files
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ genotype)

# Data transformations and visualization
vsd <- vst(ddsTxi, blind=F)

# save files
saveRDS(ddsTxi, file.path(output_dir, "dsim_Txi.rds"))
write.csv(assay(vsd), file.path(output_dir, "dsim_vst.csv"))
