#!/bin/bash

# This is the main workflow for the woodchuck genome annotation.
# Calls upon programs and files that are either listed in the README file or located in the repository.

# The woodchuck genome is WCK01_AAH20201022_F8-SCF.fasta

# First, perform annotation liftover using LiftOff from closely related species

# Download RefSeq Alpine marmot genome (marMar2.1)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/458/135/GCF_001458135.1_marMar2.1/GCF_001458135.1_marMar2.1_genomic.fna.gz
 gunzip GCF_001458135.1_marMar2.1_genomic.fna.gz
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
 GCF_001458135.1_marMar2.1_genomic.fna \
 -g GCF_001458135.1_marMar2.1_genomic.gff \
 -o from_marMar_copies_scf.gff -p 15 \
 -f liftoffFeatures.txt \
 -m /path/to/minimap2 -flank 0.5 -copies

# Repeat the annotation liftover with the yellow-bellied marmot RefSeq genome annotation (GSC_YBM_2.0)

# Download genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/676/075/GCF_003676075.2_GSC_YBM_2.0/GCF_003676075.2_GSC_YBM_2.0_genomic.fna.gz
gunzip GCF_003676075.2_GSC_YBM_2.0_genomic.fna.gz
# Download the annotation file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/676/075/GCF_003676075.2_GSC_YBM_2.0/GCF_003676075.2_GSC_YBM_2.0_genomic.gff.gz
gunzip GCF_003676075.2_GSC_YBM_2.0_genomic.gff.gz

# Run LiftOff using the same parameters
liftoff \
 WCK01_AAH20201022_F8-SCF.fasta \
 GCF_003676075.2_GSC_YBM_2.0_genomic.fna \
 -g GCF_003676075.2_GSC_YBM_2.0_genomic.gff \
 -o from_gsc_ybm_scf.gff -p 15 \
 -f liftoffFeatures.txt \
 -m /path/to/minimap2 -flank 0.5 -copies
 
 # Repeat this again with the 13-lined ground squirrel Refseq genome annotation (SpeTri2.0)
 
 # Download genome
 wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/236/235/GCF_000236235.1_SpeTri2.0/GCF_000236235.1_SpeTri2.0_genomic.fna.gz
 gunzip GCF_000236235.1_SpeTri2.0_genomic.fna.gz
 # Download annotation
 wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/236/235/GCF_000236235.1_SpeTri2.0/GCF_000236235.1_SpeTri2.0_genomic.gff.gz
 gunzip GCF_000236235.1_SpeTri2.0_genomic.gff.gz
 
 # Run LiftOff using the same parameters
liftoff \
 WCK01_AAH20201022_F8-SCF.fasta \
 GCF_000236235.1_SpeTri2.0_genomic.fna \
 -g GCF_000236235.1_SpeTri2.0_genomic.gff \
 -o from_SpeTri_scf.gff -p 15 \
 -f liftoffFeatures.txt \
 -m /path/to/minimap2 -flank 0.5 -copies
