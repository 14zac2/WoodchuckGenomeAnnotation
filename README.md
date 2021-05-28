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
- LiftOff vXX
