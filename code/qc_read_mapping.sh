#!/bin/bash
echo $BASH_VERSION
set -ex

# required software: fastqc (v0.11.9), cutadapt (v3.4), bowtie2 (v2.4.4), kallisto (v0.46.2), multiqc (1.11), and parallel (v20210222)
# sample names should be in data/runids.txt

# activate conda enviroment with required software
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" # initializes conda in sub-shell
conda activate transcriptomics
conda info|egrep "conda version|active environment"

# create directories
mkdir data/qc1 # qc results for raw reads
mkdir data/qc2 # qc results for filtered reads
mkdir data/fq_trim # trimmed reads
mkdir data/fq_filt # filtered reads
mkdir data/sam # sam files produced during bowtie2 alignment; will be deleted
mkdir data/kallisto # kallisto mapping results

# run fastqc on raw reads
cat data/runids.txt | parallel "fastqc data/fastq/{}.fastq.gz --outdir data/qc1"

# trim adapters using cutadapt
cat data/runids.txt | parallel "cutadapt -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTC \
-o data/fq_trim/{}.fastq.gz data/fastq/{}.fastq.gz \
&> data/fq_trim/{}.fastq.qz.log"

# filter reads mapping to rRNA and tRNA genes
# build bowtie2 index
bowtie2-build data/refs/Binfantis_ATCC15697_rRNA_tRNA.fasta data/refs/Binfantis_ATCC15697_rRNA_tRNA
# align reads via bowtie2; save ones that did not align to a separate file
cat data/runids.txt | parallel "bowtie2 -x data/refs/Binfantis_ATCC15697_rRNA_tRNA \
-U data/fq_trim/{}.fastq.gz \
-S data/sam/{}.sam \
--un data/fq_filt/{}.fastq \
&> data/fq_filt/{}.log"
rm -rf data/sam
cat data/runids.txt | parallel "gzip data/fq_filt/{}.fastq"

# run fastqc on filtered reads
cat data/runids.txt | parallel "fastqc data/fq_filt/{}.fastq.gz --outdir data/qc2"

# pseudolalign reads to transcriptome
# build kallisto index
kallisto index -i data/refs/Binfantis_ATCC15697_transcriptome.index data/refs/Binfantis_ATCC15697_transcriptome.fasta
# map reads to indexed reference via kallisto
cat data/runids.txt | parallel "kallisto quant -i data/refs/Binfantis_ATCC15697_transcriptome.index \
-o data/kallisto/{} \
--single \
-l 200 \
-s 20 \
data/fq_filt/{}.fastq.gz \
&> data/kallisto/{}_2.log"

# run multiqc
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
multiqc -d . -o data
