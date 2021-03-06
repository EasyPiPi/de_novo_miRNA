#### Trim adapters ####
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
        others = "--minimum-length 10 -q 20"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 4 # set desired number of threads here
    wrapper:
        "0.64.0/bio/cutadapt/pe"

#### Salmon quantification ####
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

# run DESeq2
rule deseq2:
    input:
        expand("outputs/salmon/{df.species}/{df.sample}/quant.sf", df = metadata_RNAseq.itertuples()),
        rules.download_annotation_complete.output
    params:
        root_dir = config["root_dir"]
    output:
        touch("indicator/DESeq2/all.done")
    threads: 4
    script:
        "../scripts/RNAseq/run_deseq2.R"

# get dmel ortholog for GO and mis-regulated gene analysis
rule get_ortholog:
    input:
        dmel_orthologs = rules.download_flybase_info.output.ortholog
    output:
        dmel_dsim_orthologs = "external_resources/flybase/dmel_dsim_orthologs.csv"
    script:
        "../scripts/RNAseq/get_ortholog.R"

# run GO analysis
rule run_GO_analysis:
    input:
        ortholog = rules.get_ortholog.output,
        complete = rules.deseq2.output
    params:
        root_dir = config["root_dir"]
    output:
        go = "outputs/GO/table/GO.csv"
    script:
        "../scripts/RNAseq/run_GO.R"

# compare mis-regulated genes between D. mel and D.sim
rule compare_misregulated_genes:
    input:
        ortholog = rules.get_ortholog.output,
        complete = rules.deseq2.output
    params:
        root_dir = config["root_dir"],
        figure_out_dir = "outputs/misregulated_genes/figure"
    output:
        touch("indicator/DESeq2/compare_misregulated_genes.done")
    script:
        "../scripts/RNAseq/compare_misregulated_genes.R"

#### hisat2 ####
rule hisat2_index:
    input:
        fasta = "external_resources/{species}/genome.fasta"
    output:
        directory("external_resources/{species}/hisat2_index")
    params:
        prefix = "external_resources/{species}/hisat2_index/genome"
    log:
        "logs/hisat2/index/{species}.log"
    threads: 4
    wrapper:
        "0.65.0/bio/hisat2/index"

def get_hisat2_genome_index(wildcards):
    species = metadata_RNAseq.loc[wildcards.sample, "species"]
    return os.path.join("external_resources", species + "/hisat2_index/genome")

rule hisat2_align:
    input:
        expand("external_resources/{species}/hisat2_index", species = ["dme", "dsi"]),
        fastq1 = rules.cutadapt.output.fastq1,
        fastq2 = rules.cutadapt.output.fastq2
    output:
        temp("raw_data/tmp/sam/{sample}.sam")
    threads: 4
    params:
        genome = get_hisat2_genome_index,
        others = ""
    log:
        "logs/hisat2/align/{sample}.log"
    shell:
        "hisat2 -p {threads} {params.others} -x {params.genome} -1 {input.fastq1} -2 {input.fastq2} -S {output} --summary-file {log}"

#### sort sam files using samtools ####
rule samtools_sort:
    input:
        rules.hisat2_align.output
    output:
        "outputs/hisat2/bam/{sample}.bam"
    threads:4
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule samtools_build_index:
    input:
        rules.samtools_sort.output
    output:
        "outputs/hisat2/bam/{sample}.bam.bai"
    threads:1
    shell:
        "samtools index {input} {output}"

rule generate_bam_complete:
    input:
        expand("outputs/hisat2/bam/{sample}.bam.bai", sample = metadata_RNAseq.index)
    output:
        touch("indicator/hisat2/bam/all.done")
