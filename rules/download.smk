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
        dir = config["small_RNAseq"],
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
