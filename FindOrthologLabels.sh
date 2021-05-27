#!/bin/bash

# This script labels the annotated woodchuck genome with orthologs from
# human, mouse, yellow-bellied marmot, Alpine marmot, and 13-lined ground squirrel

# Use gffread to translate the woodchuck protein-coding genes into their predicted peptide sequences
# Not including mitochondrial genes for ortholog detection
# Used gffread v0.12.3
# -y is the output in a peptide fasta file
# -g is the reference genome
# -S indicates to use * as the stop codon annotation
# The last argument is the input gff file
gffread -y mikado_final_sc2_stringent_noMito.protein.faa -g WCK01_AAH20201022_F8-SC2.fasta -S mikado_final_sc2_stringent_noMito.gff

# Acquire the predicted peptide sequences of the other organisms
# These will be used as the input with the translated woodchuck sequences for OrthoFinder2
# Human (Assembly GRCh38.p13):
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_protein.faa.gz
gunzip GRCh38_latest_protein.faa.gz
# Mouse:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_protein.faa.gz
gunzip GCF_000001635.27_GRCm39_protein.faa.gz
# Yellow-bellied marmot:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/676/075/GCF_003676075.2_GSC_YBM_2.0/GCF_003676075.2_GSC_YBM_2.0_protein.faa.gz
gunzip GCF_003676075.2_GSC_YBM_2.0_protein.faa.gz
# Alpine marmot:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/458/135/GCF_001458135.1_marMar2.1/GCF_001458135.1_marMar2.1_protein.faa.gz
gunzip GCF_001458135.1_marMar2.1_protein.faa.gz
# 13-lined ground squirrel:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/236/235/GCF_000236235.1_SpeTri2.0/GCF_000236235.1_SpeTri2.0_protein.faa.gz 
gunzip GCF_000236235.1_SpeTri2.0_protein.faa.gz

# OrthoFinder2 performs better when a phylogeny of input species is provided
# Used https://phylot.biobyte.de/ to output a phylogenetic tree from the RefSeq identifiers of the different species
# The tree needed to be rooted, so manually modified this tree to create preOrthoFinderTreeProtNames3.txt which is in the Github repo

# Used OrthoFinder v2.5.2
# Relied on many dependencies; used the following:
# Diamond v2.0.4, mcl v14.137, fastme v2.1.6.2, python v3.8.2
# -t and -a are parallelization specifications
# -s points towards the species tree
# -f points towards the directory in which the peptide sequence fasta files were contained
python ~/OrthoFinder_source/orthofinder.py -t 48 -a 12 \
 -s ~/orthofinder_sc2/preOrthoFinderTreeProtNames3.txt -f protein_seqs
