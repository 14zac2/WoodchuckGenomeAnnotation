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

## Data used
