configfile: "config.yaml"

import os
import pandas as pd
import numpy as np

##### metadata #####
metadata_small_RNAseq = pd.read_csv("metadata/small_RNAseq/metadata.csv", dtype=str)

rule all:
    input:
        #### download ####
        "indicator/download/small_RNAseq/all.done"

##### load rules #####
include: "rules/download.smk"
