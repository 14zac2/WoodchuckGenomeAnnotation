#!/bin/bash

# This is the main workflow for the woodchuck genome annotation.
# Calls upon programs and files that are either listed in the README file or located in the repository.

# The woodchuck genome is WCK01_AAH20201022_F8-SCF.fasta

# First, perform annotation liftover using LiftOff from closely related species

# Download RefSeq Alpine marmot genome and output to a .fa file (recognized by other tools whereas .fna extensions sometimes are not)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/458/135/GCF_001458135.1_marMar2.1/GCF_001458135.1_marMar2.1_genomic.fna.gz \
 -O GCF_001458135.1_marMar2.1_genomic.fa.gz
 gunzip GCF_001458135.1_marMar2.1_genomic.fa.gz
 # Download the corresponding annotation file
 wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/458/135/GCF_001458135.1_marMar2.1/GCF_001458135.1_marMar2.1_genomic.gff.gz
 gunzip GCF_001458135.1_marMar2.1_genomic.gff.gz

# Perform liftover of the Alpine marmot genome annotation to the woodchuck genome sequence
# Uses LiftOff v1.5.1, Minimap2 v2.17-r941 and Python v3.6.11
# Arguments:
# -g is the Alpine marmot annotation file to liftover
# -o is the name of the output annotation
# -p is number of threads
# -f are the feature types to liftover; includes gene, mRNA, exon, CDS, and lnc_RNA
# -flank: amount of flanking sequence to align as a fraction [0.0-1.0] of gene length. This can improve gene alignment where gene structure differs between target and reference
# -copies: look for extra gene copies in the target genome
liftoff \
 WCK01_AAH20201022_F8-SCF.fasta \
 GCF_001458135.1_marMar2.1_genomic.fa \
 -g GCF_001458135.1_marMar2.1_genomic.gff \
 -o from_marMar_copies.gff -p 15 \
 -f liftoffFeatures \
 -m /path/to/minimap2 -flank 0.5 -copies
