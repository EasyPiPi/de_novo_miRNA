suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(data.table)))

# utr3_in <- '/home/yixin/Desktop/github_repo/de_novo_miRNA/external_resources/dme/utr3.fasta'
# taxonomy_id <- 7227
# utr3_out <- "utr.tab"

utr3_in <- snakemake@input[["utr3_fasta"]]
taxonomy_id <- snakemake@params[["taxonomy_id"]]
utr3_out <- snakemake@output[["utr3_tab"]]

utr3_fa <- readDNAStringSet(utr3_in)

seq_df <- as.data.table(utr3_fa)
colnames(seq_df) <- "sequence"

seq_df[, "transcript_id" := substr(names(utr3_fa), 1, 11)]
seq_df[, "taxonomy_id" := taxonomy_id]

seq_df <- seq_df[, c("transcript_id", "taxonomy_id", "sequence")]

fwrite(seq_df, file = utr3_out, sep = "\t", col.names = FALSE)
