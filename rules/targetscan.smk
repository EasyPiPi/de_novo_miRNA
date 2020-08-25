rule parse_UTR_sequence:
    input:
        utr3_fasta = rules.download_utr_gff_gtf.output.utr3
    params:
        taxonomy_id = lambda wildcards: metadata_annotation.loc[wildcards.species, "taxonomy_id"]
    output:
        utr3_tab = "external_resources/{species}/utr3.tab"
    script:
        "../scripts/targetscan_70/parse_UTR_sequence.R"

rule run_targetscan:
    input:
        utr3 = rules.parse_UTR_sequence.output.utr3_tab,
        miR_seed = config["miR_seed"]
    output:
        "outputs/miRNA_targets/table/{species}/targetscan_targets.tab"
    shell:
        "scripts/targetscan_70/targetscan_70.pl {input.miR_seed} {input.utr3} {output}"

rule analyze_miRNA_target_expression:
    input:
        expand("outputs/miRNA_targets/table/{species}/targetscan_targets.tab", species = ["dme", "dsi"])
    params:
        root_dir = config["root_dir"]
    output:
        touch("indicator/targetScan/all.done")
    script:
        "../scripts/RNAseq/targetScan_targets_gene_level.R"
