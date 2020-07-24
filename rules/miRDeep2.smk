rule trim_galore_se:
    input:
        "raw_data/small_RNAseq/{run}.fastq.gz"
    output:
        "raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq",
        "raw_data/small_RNAseq/trim_galore/{run}.fastq.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20 --length 18 --dont_gzip"
    log:
        "logs/trim_galore/{run}.log"
    wrapper:
        "0.63.0/bio/trim_galore/se"

rule bowtie_build:
    input:
        genome = rules.download_genome_miRNA_fasta.output.genome
    output:
        touch("indicator/miRDeep2/bowtie_build/{species}.done")
    threads:4
    shell:
        "bowtie-build --threads {threads} {input.genome} genome"
