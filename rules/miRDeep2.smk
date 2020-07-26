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
    params:
        "external_resources/{species}/genome"
    output:
        touch("indicator/miRDeep2/bowtie_build/{species}.done")
    threads:4
    shell:
        "bowtie-build --threads {threads} {input.genome} {params}"

rule miRDeep2_preparation:
    input:
        expand("indicator/miRDeep2/bowtie_build/{species}.done", species = metadata_annotation.index),
        expand("raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq", run = metadata_small_RNAseq.run)
    output:
        touch("indicator/miRDeep2/preparation_for_miRDeep2/all.done")

rule miRDeep2_mapper:
    input:
        complete = "indicator/miRDeep2/preparation_for_miRDeep2/all.done",
        fastq = "raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq"
    params:
        genome = lambda wildcards: os.path.join("external_resources", metadata_small_RNAseq.loc[wildcards.run, "species"], "genome"),
        others = "-e -h -i -j -l 18 -m -v -n"
    output:
        fasta = "raw_data/small_RNAseq/miRDeep2/{run}/reads.fa"
        # arf = "raw_data/small_RNAseq/miRDeep2/{run}/reads_vs_genome.arf"
    shell:
        "mapper.pl {input.fastq} {params.others} -p {params.genome} -s {output.fasta}"
