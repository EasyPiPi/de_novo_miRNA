# # download small RNA-seq libraries
# rule download_small_RNAseq:
#     params:
#         raw_data_dir = "/home/yixin/Desktop/github_repo/de_novo_miRNA/raw_data",
#         sra_accessions = metadata_small_RNAseq.run.values
#     output:
#         "indicator/download/small_RNAseq.complete"
#     script:
#         "../scripts/download/small_RNAseq.R"

rule download_small_RNAseq:
    output:
        touch("indicator/download/small_RNAseq/{run}.done")
    params:
        dir = os.path.join(config["root_dir"], "raw_data/small_RNAseq"),
        run = "{run}"
    shell:
        """
        mkdir -p {params.dir}
        fastq-dump -O {params.dir} --gzip -A {params.run}
        """

rule download_small_RNAseq_complete:
    input:
        expand("indicator/download/small_RNAseq/{run}.done", run = metadata_small_RNAseq.run)
    output:
        "indicator/download/small_RNAseq/all.done"
    shell:
        "touch {output}"

rule download_miRNA_sequence:
    output:
        touch("indicator/download/miRNA_sequence/miRNA_sequence.done")
    params:
        mature = config["miRBase_mature"],
        hairpin = config["miRBase_hairpin"]
    shell:
        """
        mkdir -p external_resources/miRNA
        wget {params.mature} -O external_resources/miRNA/mature.fa.gz
        gzip -d external_resources/miRNA/mature.fa.gz
        wget {params.hairpin} -O external_resources/miRNA/hairpin.fa.gz
        gzip -d external_resources/miRNA/hairpin.fa.gz
        """

rule download_miRNA_gff:
    output:
        "external_resources/{species}/miRNA.gff3"
    params:
        url = lambda wildcards: metadata_annotation.loc[wildcards.species, "miRNA_gff"]
    threads:1
    shell:
        "wget -c {params.url} -O {output}"

rule download_genome_miRNA_fasta:
    output:
        genome = "external_resources/{species}/genome.fasta",
        miRNA = "external_resources/{species}/miRNA.fasta"
    params:
        genome = lambda wildcards: metadata_annotation.loc[wildcards.species, "genome"],
        miRNA = lambda wildcards: metadata_annotation.loc[wildcards.species, "miRNA"]
    threads:1
    shell:
        """
        wget -O - {params.genome} | gunzip -c > {output.genome}
        wget -O - {params.miRNA} | gunzip -c > {output.miRNA}
        """

rule download_utr_gff_gtf:
    output:
        utr3 = "external_resources/{species}/utr3.fasta",
        gff = "external_resources/{species}/all.gff",
        gtf = "external_resources/{species}/all.gtf"
    params:
        utr3 = lambda wildcards: metadata_annotation.loc[wildcards.species, "utr3"],
        gff = lambda wildcards: metadata_annotation.loc[wildcards.species, "gff"],
        gtf = lambda wildcards: metadata_annotation.loc[wildcards.species, "gtf"]
    threads:1
    shell:
        """
        wget -O - {params.utr3} | gunzip -c > {output.utr3}
        wget -O - {params.gff} | gunzip -c > {output.gff}
        wget -O - {params.gtf} | gunzip -c > {output.gtf}
        """

rule download_annotation_complete:
    input:
        expand("external_resources/{species}/miRNA.gff3", species = metadata_annotation.index),
        expand("external_resources/{species}/genome.fasta", species = metadata_annotation.index),
        expand("external_resources/{species}/all.gff", species = ["dme", "dsi"])
    output:
        touch("indicator/download/annotation/all.done")
