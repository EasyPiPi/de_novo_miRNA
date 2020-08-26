rule summarize_phenotype:
    params:
        root_dir = config["root_dir"],
        figure_out_dir = "outputs/phenotype/figure"
    output:
        touch("indicator/phenotype/all.done")
    script:
        "../scripts/phenotype/summarize_phenotype.R"
