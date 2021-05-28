# WoodchuckGenomeAnnotation

The code in the HealthyWoodchuck repository outlines the workflow for the woodchuck genome annotation project. The R and bash scripts provided are not intended to be run, but instead are guidelines of the process that contain functional code. This is to insure reproducibility of the annotation and may also be helpful for others undergoing similar projects.

## Workflow

The logical flow of the scripts is as follows:
1. AnnotateGenome.sh
2. FindOrthologLabels.sh
3. getOrthologs.R
4. makeHomologeneFile.R
5. AssignLabelsToGenome.sh

## Explanations of each script

1. AnnotateGenome.sh

This begins with the raw woodchuck genome sequence and uses annotation and transcript filtering tools to obtain the final woodchuck genome annotation. It uses the files with the .yaml extensions in the Github repo and produces a GFF file.

2. FindOrthologLabels.sh

This takes the output from AnnotateGenome.sh, the final woodchuck GFF file, and translates the gene sequences. These gene sequences are used as input with the peptide sequences of human, mouse, yellow-bellied marmot, Alpine marmot, and 13-lined ground squirrel for OrthoFinder2.

3. getOrthologs.R

This script reads in the pairwise orthologs from OrthoFinder2, which are then collected into a table to label the woodchuck genes. This table is also provided in the Github repo as collectedOrthofinderPairings.tsv.

4. makeHomologeneFile.R

This uses the pairwise human-woodchuck ortholog output from OrthoFinder2 to create a homolgene-formatted file which is used to convert human gene pathway files into woodchuck gene pathway files.

5. AssignLabelsToGenome.sh

This uses the table collectedOrthofinderPairings.tsv to label the woodchuck genes with a unique orthologous gene symbol. The mitochondrial genome annotations are also added here.

## Dependencies

The following tools are used:
- LiftOff v1.5.1
- Minimap2 v2.17-r941
- Python v3.6.10, v3.6.11, v3.8.2
- hisat2 v2.2.1
- SRA-toolkit v2.10.8
- Hisat2 v2.2.1
- Samtools v1.12
- Stringtie v2.1.3
- Portcullis v1.2.0
- Mikado v2.0rc2
- Transdecoder v5.5.0
- Blast+ v2.11.0
- Gffread v0.12.1, v0.12.3
- OrthoFinder v2.5.2
- Diamond v2.0.4
- Mcl v14.137
- Fastme v2.1.6.2
- R v4.0.4
- Dplyr v1.0.6
- Tidyverse v1.3.1

## Data

This workflow takes advantage of publicly available data. This includes the genome fasta, GFF/GTF, and translated CDS RefSeq files for the following organisms:
- Human, _Homo sapiens_, GRCh38.p13/GCF_000001405.39
- Mouse, _Mus musculus_, GRCm39/GCF_000001635.27
- Yellow-bellied marmot, _Marmota flaviventris_, GSC_YBM_2.0/GCF_003676075.2
- Alpine marmot, _Marmota marmota marmota_, marMar2.1/GCF_001458135.1
- 13-lined ground squirrel, _Ictidomys tridecemlineatus_, SpeTri2.0/GCF_000236235.1

In addition, the following RNA-seq data were used: woodchuck 250bp paired-end RNA-seq from GEO (GSE137911). One run was selected per available tissue type:
- SRR10172922 (Liver)
- SRR10172923 (Kidney)
- SRR10172924 (Spleen)
- SRR10172925 (Lung)
- SRR10172926 (Heart)
- SRR10172930 (Pancreas)
- SRR10172931 (Thymus)

