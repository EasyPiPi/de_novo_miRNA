configfile: "config.yaml"

import os
import pandas as pd
import numpy as np

##### metadata #####
metadata_annotation = pd.read_csv("metadata/annotation/metadata.csv", dtype=str)
metadata_annotation.set_index("species", drop = False, inplace = True)
metadata_annotation.sort_index(inplace = True)

metadata_small_RNAseq = pd.read_csv("metadata/small_RNAseq/metadata.csv", dtype=str)
metadata_small_RNAseq.set_index("run", drop = False, inplace = True)
metadata_small_RNAseq.sort_index(inplace = True)

metadata_RNAseq = pd.read_csv("metadata/RNAseq/metadata.csv", dtype=str)
metadata_RNAseq.set_index("sample", drop = False, inplace = True)
metadata_RNAseq.sort_index(inplace = True)

#### main rule ####
rule all:
    input:
        #### download ####
        "indicator/download/small_RNAseq/all.done",
        "indicator/download/miRNA_sequence/miRNA_sequence.done",
        "indicator/download/annotation/all.done",
        #### QC ####
        expand("outputs/qc/fastqc/{run}.html", run = metadata_small_RNAseq.run),
        #### small RNA-seq analysis ####
        "indicator/miRDeep2/analyze_miRNA_expression/all.done",
        #### RNA-seq analysis ####
        "indicator/DESeq2/all.done",
        "external_resources/flybase/dmel_dsim_orthologs.csv",
        #### targetScan ####
        "indicator/targetScan/all.done",

##### load rules #####
include: "rules/download.smk"
include: "rules/miRDeep2.smk"
include: "rules/RNAseq.smk"
include: "rules/targetscan.smk"
