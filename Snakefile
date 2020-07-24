configfile: "config.yaml"

import os
import pandas as pd
import numpy as np

##### metadata #####
metadata_small_RNAseq = pd.read_csv("metadata/small_RNAseq/metadata.csv", dtype=str)
metadata_annotation = pd.read_csv("metadata/annotation/metadata.csv", dtype=str)
metadata_annotation.set_index("species", drop = False, inplace = True)
metadata_annotation.sort_index(inplace = True)

#### main rule ####
rule all:
    input:
        #### download ####
        "indicator/download/small_RNAseq/all.done",
        "indicator/download/miRNA_sequence/miRNA_sequence.done",
        "indicator/download/annotation/all.done",
        #### miRDeep2 ####
        expand("raw_data/small_RNAseq/trim_galore/{run}_trimmed.fq", run = metadata_small_RNAseq.run),
        expand("indicator/miRDeep2/bowtie_build/{species}.done", species = metadata_annotation.index)

##### load rules #####
include: "rules/download.smk"
include: "rules/miRDeep2.smk"
