library(seqinr)
library(stringr)
library(tibble)
library(purrr)

# hairpin_file <- "/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/miRNA/hairpin.fa"
# mature_file <- "/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/miRNA/mature.fa"
# species <- c("dme", "dsi", "dse", "der", "dvi")
# output_dir <- "/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources"

hairpin_file <- snakemake@params[["hairpin"]]
mature_file <- snakemake@params[["mature"]]
species <- snakemake@params[["species"]]
output_dir <- snakemake@params[["output_dir"]]

hairpin <- read.fasta(hairpin_file, forceDNAtolower = FALSE, as.string = TRUE)
mature <- read.fasta(mature_file, forceDNAtolower = FALSE, as.string = TRUE)

hairpin_df <- tibble(species = species, miRNA = map(species, ~ hairpin[str_detect(names(hairpin), .x)]))
mature_df <- tibble(species = species, miRNA = map(species, ~ mature[str_detect(names(mature), .x)]))

pwalk(list(hairpin_df$miRNA,
          map(hairpin_df$miRNA, names),
          file.path(output_dir, hairpin_df$species, "hairpin.fasta")),
     write.fasta
     )

pwalk(list(mature_df$miRNA,
           map(mature_df$miRNA, names),
           file.path(output_dir, mature_df$species, "mature.fasta")),
      write.fasta
)
