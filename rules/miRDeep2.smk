rule trim_galore_se:
    input:
        "raw_data/small_RNAseq/{run}.fastq.gz"
    output:
        "raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq.gz",
        "raw_data/small_RNAseq/trim_galore/{run}.fastq.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20 --length 18"
    log:
        "logs/trim_galore/{run}.log"
    wrapper:
        "0.63.0/bio/trim_galore/se"
