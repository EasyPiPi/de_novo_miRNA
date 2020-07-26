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
        hairpin = "external_resources/miRNA/hairpin.fa",
        mature = "external_resources/miRNA/mature.fa",
        species = metadata_annotation.species,
        output_dir = "de_novo_miRNA/external_resources"
    output:
        touch("indicator/miRDeep2/extract_miRNA_sequences/all.done")
    shell:
        "../scripts/miRDeep2/extract_miRNA_sequences.R"

rule miRDeep2_preparation:
    input:
        # expand("indicator/miRDeep2/bowtie_build/{species}.done", species = metadata_annotation.index),
        # expand("raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq", run = metadata_small_RNAseq.run),
        "indicator/miRDeep2/extract_miRNA_sequences/all.done"
    output:
        touch("indicator/miRDeep2/preparation_for_miRDeep2/all.done")

rule miRDeep2_mapper:
    input:
        fastq = "raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq"
    params:
        # genome = lambda wildcards: os.path.join("external_resources", metadata_small_RNAseq.loc[wildcards.run, "species"], "genome"),
        others = "-e -h -i -j -l 18 -m -v -n"
    output:
        fasta = "raw_data/small_RNAseq/miRDeep2/{run}/reads.fa"
        # arf = "raw_data/small_RNAseq/miRDeep2/{run}/reads_vs_genome.arf"
    shell:
        # "mapper.pl {input.fastq} {params.others} -p {params.genome} -s {output.fasta}"
        "mapper.pl {input.fastq} {params.others} -s {output.fasta}"

rule miRDeep2_quantifier:
    input:
        complete = rules.miRDeep2_preparation.output
    params:
        output_dir = os.path.join(config["root_dir"], "outputs/miRDeep2/quantifier/{run}"),
        hairpin = lambda wildcards: os.path.join(config["root_dir"], "external_resources", metadata_small_RNAseq.loc[wildcards.run, "species"], "hairpin.fasta"),
        mature = lambda wildcards: os.path.join(config["root_dir"], "external_resources", metadata_small_RNAseq.loc[wildcards.run, "species"], "mature.fasta"),
        read = os.path.join(config["root_dir"], "raw_data/small_RNAseq/miRDeep2/{run}/reads.fa"),
        species = lambda wildcards: metadata_small_RNAseq.loc[wildcards.run, "species"],
        others = "-y now -j -g 4 -n"
    output:
        touch("indicator/miRDeep2/quantifier/{run}.done")
    shell:
        """
        mkdir -p {params.output_dir} && cd {params.output_dir}
        quantifier.pl {params.others} -p {params.hairpin} -m {params.mature} -r {params.read} -t {params.species}
        """
