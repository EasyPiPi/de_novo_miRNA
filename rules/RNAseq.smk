# Trim adapters
rule cutadapt:
    input:
        ["raw_data/RNAseq/{sample}_combined_R1.fastq.gz", "raw_data/RNAseq/{sample}_combined_R2.fastq.gz"]
    output:
        fastq1="raw_data/RNAseq/trimmed/{sample}_combined_R1.fastq.gz",
        fastq2="raw_data/RNAseq/trimmed/{sample}_combined_R2.fastq.gz",
        qc="outputs/qc/cutadapt/{sample}.qc.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        others = "--minimum-length 1 -q 20"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 4 # set desired number of threads here
    wrapper:
        "0.64.0/bio/cutadapt/pe"

# Salmon quantification
# Build index
rule salmon_index:
    input:
        "external_resources/{species}/cdna.fasta"
    output:
        directory("external_resources/{species}/salmon_index")
    log:
        "logs/salmon/{species}_index.log"
    threads: 4
    params:
        # optional parameters
        extra=""
    wrapper:
        "0.64.0/bio/salmon/index"

# Quantification
rule salmon_quant_reads:
    input:
        r1 = rules.cutadapt.output.fastq1,
        r2 = rules.cutadapt.output.fastq2,
        index = "external_resources/{species}/salmon_index"
    output:
        quant = 'outputs/salmon/{species}/{sample}/quant.sf',
        lib = 'outputs/salmon/{species}/{sample}/lib_format_counts.json'
    log:
        'logs/salmon/{species}/{sample}.log'
    params:
        # optional parameters
        libtype ="A",
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra = "" # --seqBias --gcBias
    threads: 4
    wrapper:
        "0.64.0/bio/salmon/quant"
