suppressWarnings(suppressMessages(library(SRAdb)))

# sra_accessions <-
#     c('SRR069836', 'SRR016854', 'SRR060643', 'SRR060645', 'SRR060650',
#       'SRR060652', 'SRR060657', 'SRR060659', 'SRR060664', 'SRR060666',
#       'SRR060670', 'SRR060683', 'SRR016854', 'SRR453265', 'SRR453266',
#       'SRR902009', 'SRR1205790', 'SRR1205791', 'SRR1205793', 'SRR5410380',
#       'SRR5461108', 'SRR5461110', 'SRR5461113', 'SRR6667443', 'SRR6667444')
# raw_data_dir <- "/home/yixin/Desktop/github_repo/de_novo_miRNA/raw_data"

sra_accessions <- snakemake@params[["sra_accessions"]]
raw_data_dir <- snakemake@params[["raw_data_dir"]]

fastq_dir <- file.path(raw_data_dir, "small_RNAseq")
sql_dir <- file.path(raw_data_dir, "sql")

dir.create(sql_dir, showWarnings = FALSE, recursive = TRUE)

# check local file before downloading sql file
if (file.exists(file.path(sql_dir, "SRAmetadb.sqlite"))) {
    sqlfile <- file.path(destdir = sql_dir, "SRAmetadb.sqlite")
} else {
    sqlfile <- getSRAdbFile(destdir = sql_dir)
}

sra_con <- dbConnect(SQLite(),sqlfile)
# download files one by one
for (accession in sra_accessions) {
    if (!file.exists(paste0(file.path(fastq_dir, accession), ".fastq.gz"))) {
        print(paste('Now is downloading', accession))
        try(getSRAfile(accession, destDir = fastq_dir,
                   sra_con, fileType = 'fastq', makeDirectory = TRUE))

    }
}

dbDisconnect(sra_con)
