rule parse_UTR_sequence:
    input:
        utr3_fasta = rules.download_utr_gff_gtf.output.utr3
    params:
        taxonomy_id = lambda wildcards: metadata_annotation.loc[wildcards.species, "taxonomy_id"]
    output:
        utr3_tab = "external_resources/{species}/utr3.tab"
    script:
        "../scripts/targetscan_70/parse_UTR_sequence.R"
