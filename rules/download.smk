# download small RNA-seq libraries
rule download_small_RNAseq:
    params:
        raw_data_dir = "/home/yixin/Desktop/github_repo/de_novo_miRNA/raw_data",
        sra_accessions = metadata_small_RNAseq.run.values
    output:
        "indicator/download/small_RNAseq.complete"
    script:
        "../scripts/download/small_RNAseq.R"
