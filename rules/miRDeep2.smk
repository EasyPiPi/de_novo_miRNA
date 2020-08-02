rule fastqc:
    input:
        "raw_data/small_RNAseq/{run}.fastq.gz"
    output:
        html="outputs/qc/fastqc/{run}.html",
        zip="outputs/qc/fastqc/{run}_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/{run}.log"
    threads: 1
    wrapper:
        "0.63.0/bio/fastqc"

# rule trim_galore_se:
#     input:
#         "raw_data/small_RNAseq/{run}.fastq.gz"
#     output:
#         "raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq",
#         "raw_data/small_RNAseq/trim_galore/{run}.fastq.gz_trimming_report.txt"
#     params:
#         extra="--illumina -q 20 --length 18 --dont_gzip"
#     log:
#         "logs/trim_galore/{run}.log"
#     wrapper:
#         "0.63.0/bio/trim_galore/se"

rule trimmomatic:
    input:
        "raw_data/small_RNAseq/{run}.fastq.gz"  # input and output can be uncompressed or compressed
    output:
        "raw_data/small_RNAseq/trimmomatic/{run}.fastq"
    log:
        "logs/trimmomatic/{run}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:Adapter.fa:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15",
        "MINLEN:10"],
        # ,
        # optional parameters
        extra="",
        # optional compression levels from -0 to -9 and -11
        compression_level="-9"
    threads:4
    wrapper:
        "0.63.0/bio/trimmomatic/se"

# rule bowtie_build:
#     input:
#         genome = rules.download_genome_miRNA_fasta.output.genome
#     params:
#         "external_resources/{species}/genome"
#     output:
#         touch("indicator/miRDeep2/bowtie_build/{species}.done")
#     threads:4
#     shell:
#         "bowtie-build --threads {threads} {input.genome} {params}"

rule extract_miRNA_sequences:
    params:
        hairpin = os.path.join(config["root_dir"], "external_resources/miRNA/hairpin.fa"),
        mature = os.path.join(config["root_dir"], "external_resources/miRNA/mature.fa"),
        species = metadata_annotation.species,
        output_dir = os.path.join(config["root_dir"], "external_resources")
    output:
        touch("indicator/miRDeep2/extract_miRNA_sequences/all.done")
    script:
        "../scripts/miRDeep2/extract_miRNA_sequences.R"

rule miRDeep2_mapper:
    input:
        # fastq = "raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq",
        fastq = rules.trimmomatic.output
    params:
        # genome = lambda wildcards: os.path.join("external_resources", metadata_small_RNAseq.loc[wildcards.run, "species"], "genome"),
        others = "-e -h -i -j -l 18 -m -v -n"
    output:
        fasta = "raw_data/small_RNAseq/miRDeep2/{run}/reads.fa"
        # arf = "raw_data/small_RNAseq/miRDeep2/{run}/reads_vs_genome.arf"
    shell:
        # "mapper.pl {input.fastq} {params.others} -p {params.genome} -s {output.fasta}"
        "mapper.pl {input.fastq} {params.others} -s {output.fasta}"

rule miRDeep2_preparation:
    input:
        # expand("indicator/miRDeep2/bowtie_build/{species}.done", species = metadata_annotation.index),
        # expand("raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq", run = metadata_small_RNAseq.run),
        expand("raw_data/small_RNAseq/miRDeep2/{run}/reads.fa", run = metadata_small_RNAseq.run),
        complete = rules.extract_miRNA_sequences.output
    output:
        touch("indicator/miRDeep2/preparation_for_miRDeep2/all.done")

rule miRDeep2_quantifier:
    input:
        rules.miRDeep2_preparation.output
    params:
        reads = os.path.join(config["root_dir"], "raw_data/small_RNAseq/miRDeep2/{run}/reads.fa"),
        output_dir = os.path.join(config["root_dir"], "outputs/miRDeep2/quantifier/{run}"),
        hairpin = lambda wildcards: os.path.join(config["root_dir"], "external_resources", metadata_small_RNAseq.loc[wildcards.run, "species"], "hairpin.fasta"),
        mature = lambda wildcards: os.path.join(config["root_dir"], "external_resources", metadata_small_RNAseq.loc[wildcards.run, "species"], "mature.fasta"),
        species = lambda wildcards: metadata_small_RNAseq.loc[wildcards.run, "species"],
        others = "-y now -g 3"
    output:
        touch("indicator/miRDeep2/quantifier/{run}.done")
    shell:
        """
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        quantifier.pl -p {params.hairpin} -m {params.mature} -r {params.reads} -t {params.species} {params.others}
        """
