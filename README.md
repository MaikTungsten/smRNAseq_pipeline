# smRNA-seq Processing Pipeline

## Background: smRNA-seq

Single-microbe RNA sequencing (smRNA-seq) has emerged as a technology to monitor up to thousands of single microbial transcriptomes. For data integration and further utilization

## Principle and Possibilites

This pipeline employs STARsolo for processing (mapping and quantifying) data from:
- microSPLiT
- PETRI-seq

For processing MATQ-seq data, raw reads are trimmed with bbduk, followed by mapping with either STAR or HISAT and quantification with featureCounts.

## Installation & Usage

### Installation

Install using conda or mamba, e.g., 

```bash
mamba env create -f smRNAseq.yml
```

Please make sure to provide the correct path to the conda/mamba environments in the Nextflow config file and that the environment is actually named smRNAseq.

### Usage

Typical usage at the example of MATQ-seq:

```bash
nextflow run process_smRNAseq --reads '*.fastq.gz' --matqseq \
--genomeFasta "genome.fasta" --genomeGFF "genome.gff" \
--use_hisat2 --outputDir results_MATQseq
```

Typical usage at the example of microSPLiT using two reference genomes
from two different species:

```bash
nextflow run process_smRNAseq --reads '*_{1,2}.fastq.gz' --microsplit \
--genomeFasta "genome1.fasta" "genome2.fasta" \
--genomeGFF "genome1.gff" "genome2.gff" \
--outputDir results_microSPLiT
```

### Arguments

All arguments of the pipeline are listed below. Importantly, there are general and dataset-specific parameters to specify:

#### Settings
```bash
pubDirMode                  = "copy"
```

#### General input
```bash
genomeFASTA                 = []
genomeGFF                   = []
reads                       = "*.fastq.gz"
outputDir                   = "output"
```

#### Flag for protocol-specific processing
Select depending on data supplied
Then, specific parameters and commands will be executed for processing (see below)

```bash
microsplit                  = false
matqseq                     = false
petriseq                    = false
```

#### microSPLiT-specific input

```bash
barcodes_microSPLiT_R1      = "/bin/barcodes_microSPLiT_R1.txt"
barcodes_microSPLiT_R2      = "/bin/barcodes_microSPLiT_R2.txt"
barcodes_microSPLiT_R3      = "/bin/barcodes_microSPLiT_R3.txt"
```

#### MATQ-seq-specific input

```bash
use_hisat2                  = false
use_star                    = false
bbduk_opts                  = "-Xmx16g t=8 minlen=18 qtrim=rl trimq=20 k=17 mink=11 hdist=1"
bbduk_ref_left              = "/bin/matseq_primers.fa"
bbduk_ref_right             = "/bin/nextera_and_primers.fa"
bbduk_polyA                 = 30
countFeature                = "gene"
featureIdentifier           = "ID"
```

#### PETRI-seq-specific input
```bash
barcodes_PETRIseq_R1        = "/bin/barcodes_PETRIseq_R1.txt"
barcodes_PETRIseq_R2        = "/bin/barcodes_PETRIseq_R2.txt"
barcodes_PETRIseq_R3        = "/bin/barcodes_PETRIseq_R3.txt"
```

#### Conda environments
The environment containing all required packages must be named "smRNAseq"; otherwise rename in processes accordingly
```bash
conda_path                  = "USERPATH/mamba/envs"
```

Arguments can also be set in the config file (./conf/params.config)

## Acknowledgements

We want to thank the authors of the MATQ-seq, microSPLiT and PETRI-seq workflows and corresponding datasets for their careful documentation and provided metadata and code, which allowed use to process their original data:

- microSPLiT: Kuchina et al., Science, 2020; Gaisser et al., Nature Protocols, 2024
- MATQ-seq: Homberger et al., mBio, 2023; Homberger et al., Nature Protocols, 2025
- PETRI-seq: Blattman et al., Science, 2020
