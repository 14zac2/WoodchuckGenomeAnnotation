#!/bin/bash

# This is the main workflow for the woodchuck genome annotation.
# Calls upon programs and files that are either listed in the README file or located in the repository.

# The woodchuck genome is WCK01_AAH20201022_F8-SC2.fasta

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
# -flank: amount of flanking sequence to align as a fraction [0.0-1.0] of gene length. This can improve gene alignment where gene structure differs between target and reference
# -copies: look for extra gene copies in the target genome
liftoff \
 WCK01_AAH20201022_F8-SC2.fasta \
 GCF_001458135.1_marMar2.1_genomic.fna \
 -g GCF_001458135.1_marMar2.1_genomic.gff \
 -o from_marMar_copies_sc2.gff -p 15 \
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
 WCK01_AAH20201022_F8-SC2.fasta \
 GCF_003676075.2_GSC_YBM_2.0_genomic.fna \
 -g GCF_003676075.2_GSC_YBM_2.0_genomic.gff \
 -o from_gsc_ybm_sc2.gff -p 15 \
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
 WCK01_AAH20201022_F8-SC2.fasta \
 GCF_000236235.1_SpeTri2.0_genomic.fna \
 -g GCF_000236235.1_SpeTri2.0_genomic.gff \
 -o from_SpeTri_sc2.gff -p 15 \
 -m /path/to/minimap2 -flank 0.5 -copies

# Second, use transcriptional evidence to locate transcribed genes

# The first step of this process is downloading publically available woodchuck RNA-seq data
# These reads are paired-end and 200bp which improves the accuracy of gene calling compared to short read sequences that
# seem to predict more fragmented gene structure

# Using SRA-toolkit v2.10.8
# Retrieve and the liver reads
prefetch SRR10172922
# Unpack - this provides the two paired fastq files within the reads folder SRR10172922
fasterq-dump --split-files SRR10172922.sra
# Retrieve and unpack kidney
prefetch SRR10172923
fasterq-dump --split-files SRR10172923.sra
# Retrieve and unpack spleen
prefetch SRR10172924
fasterq-dump --split-files SRR10172924.sra
# Retrieve and unpack lung
prefetch SRR10172925
fasterq-dump --split-files SRR10172925.sra
# Retrieve and unpack heart
prefetch SRR10172926
fasterq-dump --split-files SRR10172926.sra
# Retrieve and unpack pancreas
prefetch SRR10172930
fasterq-dump --split-files SRR10172930.sra
# Retrieve and unpack thymus
prefetch SRR10172931
fasterq-dump --split-files SRR10172931.sra

# Now align these reads to the woodchuck genome sequence using hisat2 v2.2.1
# Index the woodchuck genome; -p is the number of threads
# hisat2-build-s specifies the building of a small index library and creates .ht2 extensions
# This is used over the creation of a large index because the genome is less than 4 billion base pairs in length
hisat2-build-s -p 48 WCK01_AAH20201022_F8-SC2.fasta WCK01_AAH20201022_F8-SC2
# Align the reads
# hisat2-align-s specifies that a small reference index is used
# -p is the number of threads
# --dta reports alignments tailored for transcript assemblers (e.g. stringtie)
# -x is the base name of the reference genome index that was specified in hisat2-build-s (precedes the .ht2 extensions)
# -1 and -2 specify the first and second mates of paired-end reads
# -S is the output SAM alignment file
hisat2-align-s -p 48 --dta -x WCK01_AAH20201022_F8-SC2 \
 -1 SRR10172922/SRR10172922.sra_1.fastq \
 -2 SRR10172922/SRR10172922.sra_2.fastq \
 -S SRR10172922_sc2.sam
hisat2-align-s -p 48 --dta -x WCK01_AAH20201023_F8-SC2 \
 -1 SRR10172923/SRR10172923.sra_1.fastq \
 -2 SRR10172923/SRR10172923.sra_2.fastq \
 -S SRR10172923_sc2.sam
hisat2-align-s -p 48 --dta -x WCK01_AAH20201022_F8-SC2 \
 -1 SRR10172924/SRR10172924.sra_1.fastq \
 -2 SRR10172924/SRR10172924.sra_2.fastq \
 -S SRR10172924_sc2.sam
hisat2-align-s -p 48 --dta -x WCK01_AAH20201022_F8-SC2 \
 -1 SRR10172925/SRR10172925.sra_1.fastq \
 -2 SRR10172925/SRR10172925.sra_2.fastq \
 -S SRR10172925_sc2.sam
hisat2-align-s -p 48 --dta -x WCK01_AAH20201022_F8-SC2 \
 -1 SRR10172926/SRR10172926.sra_1.fastq \
 -2 SRR10172926/SRR10172926.sra_2.fastq \
 -S SRR10172926_sc2.sam
hisat2-align-s -p 48 --dta -x WCK01_AAH20201022_F8-SC2 \
 -1 SRR10172930/SRR10172930.sra_1.fastq \
 -2 SRR10172930/SRR10172930.sra_2.fastq \
 -S SRR10172930_sc2.sam
hisat2-align-s -p 48 --dta -x WCK01_AAH20201022_F8-SC2 \
 -1 SRR10172931/SRR10172931.sra_1.fastq \
 -2 SRR10172931/SRR10172931.sra_2.fastq \
 -S SRR10172931_sc2.sam
 
# These alignments now needed to be converted into sorted BAM files to be used as input for stringtie
# Using samtools v1.12
# -@ is number of threads
# -S specifies SAM input
# -h include header in SAM output; might be irrelevant in this scenario
# -u indicates uncompressed BAM output
# -o is name of output file
samtools view -@ 48 -Shu -o SRR10172922_sc2.bam SRR10172922_sc2.sam
samtools sort -@ 48 -o SRR10172922_sc2.sorted.bam SRR10172922_sc2.bam
samtools view -@ 48 -Shu -o SRR10172923_sc2.bam SRR10172923_sc2.sam
samtools sort -@ 48 -o SRR10172923_sc2.sorted.bam SRR10172923_sc2.bam
samtools view -@ 48 -Shu -o SRR10172924_sc2.bam SRR10172924_sc2.sam
samtools sort -@ 48 -o SRR10172924_sc2.sorted.bam SRR10172924_sc2.bam
samtools view -@ 48 -Shu -o SRR10172925_sc2.bam SRR10172925_sc2.sam
samtools sort -@ 48 -o SRR10172925_sc2.sorted.bam SRR10172925_sc2.bam
samtools view -@ 48 -Shu -o SRR10172926_sc2.bam SRR10172926_sc2.sam
samtools sort -@ 48 -o SRR10172926_sc2.sorted.bam SRR10172926_sc2.bam
samtools view -@ 48 -Shu -o SRR10172930_sc2.bam SRR10172930_sc2.sam
samtools sort -@ 48 -o SRR10172930_sc2.sorted.bam SRR10172930_sc2.bam
samtools view -@ 48 -Shu -o SRR10172931_sc2.bam SRR10172931_sc2.sam
samtools sort -@ 48 -o SRR10172931_sc2.sorted.bam SRR10172931_sc2.bam

# Use the sorted bam files as input for stringtie
# The output is a genome annotation file
# Using stringtie v2.1.3
# -o is the name of the output gtf file
# -p is the number of threads
# -l is the name prefix for output transcripts
stringtie SRR10172922_sc2.sorted.bam \
 -o stringtie_SRR10172922_sc2.gtf -p 48 -l STRG22
stringtie SRR10172923_sc2.sorted.bam \
 -o stringtie_SRR10172923_sc2.gtf -p 48 -l STRG23
stringtie SRR10172924_sc2.sorted.bam \
 -o stringtie_SRR10172924_sc2.gtf -p 48 -l STRG24
stringtie SRR10172925_sc2.sorted.bam \
 -o stringtie_SRR10172925_sc2.gtf -p 48 -l STRG25
stringtie SRR10172926_sc2.sorted.bam \
 -o stringtie_SRR10172926_sc2.gtf -p 48 -l STRG26
stringtie SRR10172930_sc2.sorted.bam \
 -o stringtie_SRR10172930_sc2.gtf -p 48 -l STRG30
stringtie SRR10172931_sc2.sorted.bam \
 -o stringtie_SRR10172931_sc2.gtf -p 48 -l STRG31
 
# To prepare for filtering the stringtie annotation sets with Mikado, validate splice junctions with Portcullis
# Combine separate sorted bam files for input to Portcullis
# -@ is the number of threads
samtools merge -@ 48 alioto_sc2_merged.sorted.bam \
 SRR10172922_sc2.sorted.bam \
 SRR10172923_sc2.sorted.bam \
 SRR10172924_sc2.sorted.bam \
 SRR10172925_sc2.sorted.bam \
 SRR10172926_sc2.sorted.bam \
 SRR10172930_sc2.sorted.bam \
 SRR10172931_sc2.sorted.bam

# Run Portcullis on the combined BAM to validate splice junctions
# Using portcullis v1.2.0
# The full pipeline includes prep, junc, filt and outputs validated junctions
# -t is number of threads
portcullis full -t 1 WCK01_AAH20201022_F8-SC2.fasta alioto_sc2_merged.sorted.bam
# Use the filt output portcullis.pass.junctions.bed for Mikado

# Filter the stringtie annotation with Mikado
# Using Mikado v2.0rc2 and Python v3.6.10
# conf.yaml is a separate configuration file used in the healthy woodchuck Github repo
# Mikado prepare collects the assemblies, removes redundant transcripts, and extracts their sequences
# -p is the number of threads
# --start-method indicates how to begin the multiprocessing
# --json-conf points to the configuration file
mikado prepare -p 48 --start-method spawn --json-conf mikado_round1_configurationFile.yaml

# Now BLAST+ and Transdecoder must be run on the output of mikado prepare
# The output was a fasta file and GTF file of the non-redundant transcripts

# Using Transdecoder v5.5.0; identifies likely coding sequences based on open reading frames (ORFs)
# -t is the transcript fasta file
# -m is the minimum protein length of 30 amino acids
TransDecoder.LongOrfs -t mikado_prepared.fasta -m 30
# --retain_long_orfs_length retains ORFs equal to or longer than 30 nucleotides
TransDecoder.Predict -t mikado_prepared.fasta --retain_long_orfs_length 30
# Outputs valid ORFs in a bed file, in this case: mikado_prepared.fasta.transdecoder.bed

# Use blast+ v2.11.0 to look for sequence homology in SwissProt database
# First download SwissProt database; this was accessed on Oct 4th, 2020
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
# Index the protein database so that it can be used for blast
# -dbtype indicates that this is a protein database
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot
# Blast the transcript sequences against the SwissProt database
# blastx specifies blasting nucleotides against amino acids
# -max_target_seqs: keep a maximum of 5 hits
# -query is the transcript file from Mikado that will be blasted against the protein database
# -outfmt specifies the format required by the next step of Mikado
# -db is the SwissProt database
# -evalue is a minimum measure of significance to consider a protein sequence in the SwissProt database a hit against the query
blastx -max_target_seqs 5 -num_threads 48 -query mikado_prepared.fasta \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" \
 -db uniprot_sprot -evalue 0.000001 -out blast_results.tsv

# The outputs of blast+, Transdecoder, and Portcullis are used for the next step of Mikado
# -p is the number of threads and --start-method indicates how to begin the multi-processing
# --orfs points to the output from Transdecoder
# --transcripts points to the output from Mikado prepare
# --tsv points to the results from blast+
# --json-conf points to the Mikado configuration file that was used for Mikado prepare and is in the Github repo
# --genome_fai points to the index file of the reference genome
# --log specifies the name of the output log file
# --blast-targets points towards the SwissProt fasta file that was used for the blast+ homology search
# --max-target-seqs 5 indicates the maximum target sequences allowed for the homology search
# --junctions points to the Portcullis output
mikado serialise -p 48 --start-method spawn \
 --orfs mikado_prepared.fasta.transdecoder.bed \
 --transcripts mikado_prepared.fasta --tsv blast_results.tsv \
 --json-conf mikado_round1_configurationFile.yaml --genome_fai WCK01_AAH20201022_F8-SC2.fasta.fai \
 --log mikado_serialise.log --blast-targets uniprot_sprot.fasta --max-target-seqs 5 \
 --junctions portcullis.pass.junctions.bed
# The important output is the .db file. It needs to be deleted if mikado serialise fails and needs to be run again.

# The final step of Mikado uses the information from the Mikado serialise .db file to "pick" the final annotation set
# --json-conf points towards the Mikado configuration file
# -db points towards the output from Mikado serialise
# --start-method and -p specify the multi-processing requirements
# --loci-out specifies the name of the GFF output file
# --log the output log file
# --scoring points towards the scoring file used to prioritise the transcripts in the final annotation;
# different species typically have different scoring files, as optimal exon lengths, numbers, etc. vary across species
# --mode indicates how transcripts should be split in the case of multiple blast hits or ORFs in a single transcript;
# permissive presumes that transcripts with multiple ORFs are a chimeras, and splits them, unless two ORFs share a hit in the protein database
# This is recommended for noisy RNA-seq data (https://mikado.readthedocs.io/en/stable/Tutorial/Adapting/)
mikado pick --json-conf mikado_round1_configurationFile.yaml -db mikado.db --start-method spawn -p 48 \
 --loci-out mikado_stringtie_sc2_permissive.gff --log mikado_pick_stringtie_sc2_permissive.log mikado_prepared.gtf \
 --scoring mikado_round1_scoringFile.yaml --mode permissive
 # This created mikado_stringtie_sc2_permissive.gff which can be used as input for the next round of Mikado
 
# Now use the LiftOff annotations, the filtered stringtie annotations, and the Ovaltine annotation as input for a second round of Mikado
# mikado_round2_configurationFile.yaml is the new Mikado configuration file for this round of Mikado, located in the Github repo
mikado prepare -p 48 --start-method spawn --json-conf mikado_round2_configurationFile.yaml
# Creates the transcript fasta output from the various input assemblies listed in mikado_round2_configurationFile.yaml
# Now run both blastx and Transdecoder, but not Portcullis because none of the inputs are raw transcript assemblies
# This creates the same intermediate files as before (e.g. log files, mikado_prepared.fasta, etc.)
TransDecoder.LongOrfs -t mikado_prepared.fasta -m 30
TransDecoder.Predict -t mikado_prepared.fasta --retain_long_orfs_length 30
blastx -max_target_seqs 5 -num_threads 48 -query mikado_prepared.fasta \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" \
 -db uniprot_sprot -evalue 0.000001 -out blast_results.tsv
# Now use the output of these tools and the different annotation set scores (in mikado_round2_configurationFile.yaml) for Mikado serialise
mikado serialise -p 48 --start-method spawn \
 --orfs mikado_prepared.fasta.transdecoder.bed \
 --transcripts mikado_prepared.fasta --tsv blast_results.tsv \
 --json-conf mikado_round2_configurationFile.yaml --genome_fai WCK01_AAH20201022_F8-SC2.fasta.fai \
 --log mikado_serialise.log --blast-targets uniprot_sprot.fasta --max-target-seqs 5
 # Pick the final annotation set
 # --mode is changed to stringent which means that transcripts are only split if two consecutive ORFs have both blast hits
 # and none of those hits is against the same target.
 # mikado_round2_scoringFile.yaml is the new scoring file to use to pick the transcripts
 mikado pick --json-conf conf.yaml -db mikado.db --start-method spawn -p 48 --mode stringent \
 --loci-out mikado_final_sc2_stringent.gff --log mikado_final_sc2_stringent.log mikado_prepared.gtf \
 --scoring mikado_round2_scoringFile.yaml
 # This generates the final annotation
 
 # Remove any mitochondrial annotations from the final gff file, as these will be annotated separately
 # The mitochondrial genome is the contig WCK01_MT20190704
 grep -v "WCK01_MT20190704" mikado_final_sc2_stringent.gff > mikado_final_sc2_stringent_noMito.gff
 # mikado_final_sc2_stringent_noMito.gff is the final annotation that will be used for ortholog detection
