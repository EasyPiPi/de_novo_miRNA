# Trim adapters
rule cutadapt:
    input:
        ["/local1/home/yizhao/download/miRNA/{sample}_combined_R1.fastq.gz", "/local1/home/yizhao/download/miRNA/{sample}_combined_R2.fastq.gz"]
    output:
        fastq1="/local1/home/yizhao/download/miRNA/trimmed/{sample}_combined_R1.fastq.gz",
        fastq2="/local1/home/yizhao/download/miRNA/trimmed/{sample}_combined_R2.fastq.gz",
        qc="trimmed/{sample}.qc.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters_r1 = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        adapters_r2 = "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        others = "--minimum-length 1 -q 20 --cores=4"
    log:
        "logs/cutadapt/{sample}.log"
    wrapper:
        "0.36.0/bio/cutadapt/pe"

# Salmon quantification
# Build index
rule salmon_index:
    input:
        "/sonas-hs/siepel/nlsas/data/home/yizhao/projects/Snakemake_projects/miRNA/miR983_975/cdna/{cdna}.cdna.all.filter.fa"
    output:
        directory("salmon/{cdna}_index")
    log:
        "logs/salmon/{cdna}_index.log"
    threads: 2
    params:
        # optional parameters
        extra=""
    wrapper:
        "0.36.0/bio/salmon/index"

# Quantification
rule salmon_quant_reads:
    input:
        r1 = "/local1/home/yizhao/download/miRNA/trimmed/{sample}_combined_R1.fastq.gz",
        r2 = "/local1/home/yizhao/download/miRNA/trimmed/{sample}_combined_R2.fastq.gz",
        index = "/sonas-hs/siepel/nlsas/data/home/yizhao/projects/Snakemake_projects/miRNA/miR983_975/salmon/{species_index}"
    output:
        quant = '/sonas-hs/siepel/nlsas/data/home/yizhao/projects/Snakemake_projects/miRNA/miR983_975/salmon/{species_index}/{sample}/quant.sf',
        lib = '/sonas-hs/siepel/nlsas/data/home/yizhao/projects/Snakemake_projects/miRNA/miR983_975/salmon/{species_index}/{sample}/lib_format_counts.json'
    log:
        'logs/salmon/{species_index}/{sample}.log'
    params:
        # optional parameters
        libtype ="A",
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra="--seqBias --gcBias"
    threads: 2
    wrapper:
        "0.36.0/bio/salmon/quant"
