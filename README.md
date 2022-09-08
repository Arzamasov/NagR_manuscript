# Overview of the repository
This repository contains code and analysis conducted for the manuscript "Human Milk Oligosaccharide Utilization in Intestinal Bifidobacteria is Governed by a Global Transcriptional Regulator NagR".

## Abstract
*Bifidobacterium longum* subsp. *infantis* is a prevalent beneficial bacterium that colonizes the human neonatal gut and is uniquely adapted to efficiently use human milk oligosaccharides (HMOs) as a carbon and energy source. Multiple studies have focused on characterizing the elements of HMO utilization machinery in *B. longum* subsp. *infantis*; however, the regulatory mechanisms governing the expression of these catabolic pathways remain poorly understood. A bioinformatic regulon reconstruction approach used in this study implicated NagR, a transcription factor from the ROK family, as a negative global regulator of gene clusters encoding lacto-*N*-biose/galacto-*N*-biose (LNB/GNB), lacto-*N*-tetraose (LNT), and lacto-*N*-neotetraose (LNnT) utilization pathways in *B. longum* subsp. *infantis*. This conjecture was corroborated by transcriptome profiling upon *nagR* genetic inactivation and experimental assessment of binding of recombinant NagR to predicted DNA operators. The latter approach also implicated *N*-acetylglucosamine (GlcNAc), a universal intermediate of LNT and LNnT catabolism, and its phosphorylated derivatives as plausible NagR transcriptional effectors. Reconstruction of NagR regulons in various *Bifidobacterium* lineages revealed multiple potential regulon expansion events, suggesting evolution from a local regulator of GlcNAc catabolism in ancestral bifidobacteria to a global regulator controlling the utilization of mixtures of GlcNAc-containing host glycans in *B. longum* subsp. *infantis* and *Bifidobacterium bifidum*.

## Directory structure
 - `code/` - scripts with custom R functions used in the analysis
 - `data/growth` - raw growth data
 - `data/metabolomics` - raw HMO consumption and organic acid production data
 - `data/rnaseq/` - RNA-seq data: (i) Kallisto mapping output and log files, (ii) study design file, (iii) mcSEED genome annotation files
 - `data/emsa/` - raw EMSA (gel quantification) data
 - `data/regulons/` - NagR sequences and regulon reconstruction data
 - `results/figures/` - figure drafts produced using the Rmarkdown file
 - `results/tables/` - table drafts produced using the Rmarkdown file
 - `results/evolution/` - aligments and phylogenetic trees
 - `renv/` - snapshot of the R environment captured via the renv package
 - `NagR_manuscript.Rproj` - RProject file
 - `Arzamasov_SuppCode_File.Rmd` - Rmarkdown file that combines all code and outputs and was used to generate the supplementary code file included in the manuscript

## Reproducibility and accessibility
1. All code used in the analysis is available in this repository. Once the repository has been downloaded, navigate to `NagR_manuscript/` to find the Rmarkdown document as well as the RProject file. This should be your working directory for executing code. To fully reproduce the processing of raw fastq files (RNA-seq data analysis), you will need to download them from Gene Expression Omnibus, under accession [GSE196064](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196064). Downloaded fastq files should be placed to `NagR_manuscript/data/rnaseq/fastq/`. Otherwise, `data/rnaseq/kallisto/` already includes Kallisto mapping outputs
2. To reproduce the R environment used in the analysis, use the [renv](https://rstudio.github.io/renv/articles/renv.html) package and `renv::restore()` to restore the environment from `renv.lock`

