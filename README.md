# Run or die in the evolution of new microRNAs
This is a repository of the data analysis for [Zhao et al. (2020)](https://www.biorxiv.org/content/10.1101/345769v2).

## Abstract
The Red Queen hypothesis depicts evolution as the continual struggle to adapt. According to this hypothesis, new genes, especially those originating from non-genic sequences (i.e., de novo genes), are eliminated unless they evolve continually in adaptation to a changing environment. Here, we analyze two Drosophila de novo miRNAs that are expressed in a testis-specific manner with very high rates of evolution in their DNA sequence. We knocked out these miRNAs in two sibling species and investigated their contributions to different fitness components. We observed that the fitness contributions of miR-975 in D. simulans seem positive, in contrast to its neutral contributions in D. melanogaster, while miR-983 appears to have negative contributions in both species, as the fitness of the knockout mutant increases. As predicted by the Red Queen hypothesis, the fitness difference of these de novo miRNAs indicates their different fates.

## Re-run the analysis
The analysis pipeline is managed using [Snakemake](https://snakemake.readthedocs.io/en/stable/).

To run the analysis, first revise config.yml to set working directory.

Then create a conda environment by
`conda env create --file environment.yml`

Activate conda environment by
`conda activate de_novo_miRNA`

Run the whole analysis using 8 cores by
`snakemake --cores 8 --use-conda`  
