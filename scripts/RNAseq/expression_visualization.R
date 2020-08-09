library(tidyverse)
library(skimr)
library(naniar)

library(RColorBrewer)
library(pheatmap)
library(GGally)

root_dir <- '/home/yixin/NutstoreFiles/Nutstore/Wu_lab'
exp_dir <- file.path(root_dir, 'Snakemake_projects/miR983_975/output/expression')
fig_output_dir <- file.path(exp_dir, "figure")
if (!dir.exists(fig_output_dir)) dir.create(fig_output_dir)

dme_vsd <- read_csv(file.path(exp_dir, "dmel_vst.csv")) %>% rename(gene_ID = X1)
dsi_vsd <- read_csv(file.path(exp_dir, "dsim_vst.csv")) %>% rename(gene_ID = X1)

################################################################################
# expression distance
make_heatmap <- function(df){
    colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
    df %>%
        select_if(is.numeric) %>% 
        cor() %>% 
        pheatmap(col=colors)
}

# save files
png_partial <- partial(png, width = 600, height = 500)
save_fig <- function(file_path, fig_fun){
    png_partial(file_path)
    fig_fun
    dev.off()}

save_fig(file.path(fig_output_dir, "dme_vsd_heatmap.png"), make_heatmap(dme_vsd))
save_fig(file.path(fig_output_dir, "dsi_vsd_heatmap.png"), make_heatmap(dsi_vsd))

# Pair plot
# Revise the lower plot, change its transparency
my_point <- function(data, mapping, ...) {
    ggplot(data = data, mapping=mapping) +
        geom_point(..., alpha=0.02)}

make_scatter <- function(df){
    df %>%
        select_if(is.numeric) %>% 
        ggpairs(mapping = aes(), lower = list(continuous = my_point), 
                upper = list(continuous = wrap(ggally_cor, size=3, 
                                               hjust=0.15, alignPercent=0.9)))
}

make_scatter(dme_vsd)
make_scatter(dsi_vsd)



