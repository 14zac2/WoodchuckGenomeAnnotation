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

# Before the results can be processed by OrthoFinder2, additional files need to be created associated protein IDs
# with their associated gene symbols. This is because OrthoFinder2 gives results in terms of peptide-peptide relationships
# and single-cell analysis relies on gene symbols.
# Databases like Uniprot and BioMart were not used because they do not contain complete gene symbol and protein ID information
# for all of the organisms involved in this analysis.
# Instead, conversion files created by manipulating the headers of translated CDS fasta files from Refseq that contained
# the required gene + protein information
# First, download the translated CDS fasta file for a certain species
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_translated_cds.faa.gz
gunzip GCF_000001405.39_GRCh38.p13_translated_cds.faa.gz
# Grab header of this file
grep '>' GCF_000001405.39_GRCh38.p13_translated_cds.faa > GRCh38.p13_cdsHeader.tsv
# Get rid of everything up to gene name
sed -i -e 's/>\(.*\)gene=//' GRCh38.p13_cdsHeader.tsv
# Get rid of everything between gene and protein id
sed -i -e 's/]\(.*\)protein_id=/\t/' GRCh38.p13_cdsHeader.tsv
# Get rid of everything after protein id
sed -i -e 's/].*//' GRCh38.p13_cdsHeader.tsv
# Make sure all have gene and protein
awk -F'\t' 'NF==2' GRCh38.p13_cdsHeader.tsv > temp.tsv && mv temp.tsv GRCh38.p13_cdsHeader.tsv
# This has created a file that has the gene id in the first column and the protein id in the second column

# Repeat this for mouse
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_translated_cds.faa.gz
gunzip GCF_000001635.27_GRCm39_translated_cds.faa.gz
grep '>' GCF_000001635.27_GRCm39_translated_cds.faa > GRCm39_cdsHeader.tsv
sed -i -e 's/>\(.*\)gene=//' GRCm39_cdsHeader.tsv
sed -i -e 's/]\(.*\)protein_id=/\t/' GRCm39_cdsHeader.tsv
sed -i -e 's/].*//' GRCm39_cdsHeader.tsv
awk -F'\t' 'NF==2' GRCm39_cdsHeader.tsv > temp.tsv && mv temp.tsv GRCm39_cdsHeader.tsv

# Repeat this for yellow-bellied marmot
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/676/075/GCF_003676075.2_GSC_YBM_2.0/GCF_003676075.2_GSC_YBM_2.0_translated_cds.faa.gz
gunzip GCF_003676075.2_GSC_YBM_2.0_translated_cds.faa.gz
grep '>' GCF_003676075.2_GSC_YBM_2.0_translated_cds.faa > GSC_YBM_2.0_cdsHeader.tsv
sed -i -e 's/>\(.*\)gene=//' GSC_YBM_2.0_cdsHeader.tsv
sed -i -e 's/]\(.*\)protein_id=/\t/' GSC_YBM_2.0_cdsHeader.tsv
sed -i -e 's/].*//' GSC_YBM_2.0_cdsHeader.tsv
awk -F'\t' 'NF==2' GSC_YBM_2.0_cdsHeader.tsv > temp.tsv && mv temp.tsv GSC_YBM_2.0_cdsHeader.tsv

# Repeat this for Alpine marmot
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/458/135/GCF_001458135.1_marMar2.1/GCF_001458135.1_marMar2.1_translated_cds.faa.gz
gunzip GCF_001458135.1_marMar2.1_translated_cds.faa.gz
grep '>' GCF_001458135.1_marMar2.1_translated_cds.faa > marMar2.1_cdsHeader.tsv
sed -i -e 's/>\(.*\)gene=//' marMar2.1_cdsHeader.tsv
sed -i -e 's/]\(.*\)protein_id=/\t/' marMar2.1_cdsHeader.tsv
sed -i -e 's/].*//' marMar2.1_cdsHeader.tsv
awk -F'\t' 'NF==2' marMar2.1_cdsHeader.tsv > temp.tsv && mv temp.tsv marMar2.1_cdsHeader.tsv

# Repeat this for 13-lined ground squirrel
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/236/235/GCF_000236235.1_SpeTri2.0/GCF_000236235.1_SpeTri2.0_translated_cds.faa.gz
gunzip GCF_000236235.1_SpeTri2.0_translated_cds.faa.gz
grep '>' GCF_000236235.1_SpeTri2.0_translated_cds.faa > SpeTri2.0_cdsHeader.tsv
sed -i -e 's/>\(.*\)gene=//' SpeTri2.0_cdsHeader.tsv
sed -i -e 's/]\(.*\)protein_id=/\t/' SpeTri2.0_cdsHeader.tsv
sed -i -e 's/].*//' SpeTri2.0_cdsHeader.tsv
awk -F'\t' 'NF==2' SpeTri2.0_cdsHeader.tsv > temp.tsv && mv temp.tsv SpeTri2.0_cdsHeader.tsv
# All of these files will now be read into an R script to convert the OrthoFinder peptide relationships to gene relationships

# Another file that needs to be made is a "Homologene" style file for human-woodchuck
# This will be used to convert human pathway analysis files (GMT files) to woodchuck pathway analysis files
# This file will also be made in R, but a precursor file will be made in bash
# The precursor is similar to the gene-protein cdsHeader files, but also includes the gene ID (a unique number representing the gene)
grep '>' ../../../refseq_blast/human_blast/GCF_000001405.39_GRCh38.p13_translated_cds.faa > GRCh38.p13_cdsHeader_ID.tsv
# Get rid of everything up to gene name
sed -i -e 's/>\(.*\)gene=//' GRCh38.p13_cdsHeader_ID.tsv
# Get rid of everything between gene name and gene id
sed -i -e 's/]\(.*\)GeneID:/\t/' GRCh38.p13_cdsHeader_ID.tsv
# Get rid of everything between gene id and protein id
sed -i -e 's/]\(.*\)protein_id=/\t/' GRCh38.p13_cdsHeader_ID.tsv
# Get rid of everything after protein id
sed -i -e 's/].*//' GRCh38.p13_cdsHeader_ID.tsv
# Make sure all have gene name, gene id, and protein by ensuring there are 3 columns
awk -F'\t' 'NF==3' GRCh38.p13_cdsHeader_ID.tsv > temp.tsv && mv temp.tsv GRCh38.p13_cdsHeader_ID.tsv
