# This script reads in the output of OrthoFinder2 and creates a table of orthologs to the
# genes that have been annotated in the woodchuck genome

library(dplyr)
library(tidyverse)

# Navigate to the directory of the OrthoFinder results that contain the species-species orthologs (not orthogroups)
setwd("~/OrthoFinder/Results_May24/Orthologues/Orthologues_mikado_final_sc2_stringent_noMito_protein/")

# Read in woodchuck and human pairwise orthologs from OrthoFinder
humanPairs <- read.table("mikado_final_sc2_stringent_noMito_protein__v__GRCh38_latest_protein.tsv",
                         header = TRUE,
                         sep = "\t")
# Expand by both human and woodchuck
humanPairs <- separate_rows(humanPairs, mikado_final_sc2_stringent_noMito_protein, sep = ", ")
humanPairs <- separate_rows(humanPairs, GRCh38_latest_protein, sep = ", ")
# Remove decimals from woodchuck only, because the decimals represent splice variants
humanPairs$mikado_final_sc2_stringent_noMito_protein <- gsub("\\.[0-9]+$","",
                                             humanPairs$mikado_final_sc2_stringent_noMito_protein)
# Remove orthogroup column as this is not used
humanPairs <- dplyr::select(humanPairs,
                            mikado_final_sc2_stringent_noMito_protein,
                            GRCh38_latest_protein)
# Remove duplicates
humanPairs <- dplyr::distinct(humanPairs)
# Read in conversion file that was created in the bash script
humanConversion <- read.table("GRCh38.p13_cdsHeader.tsv",
                              header = FALSE,
                              sep = "\t")
colnames(humanConversion) <- c("humanGene", "GRCh38_latest_protein")
# Combine with pairs
humanPairs <- dplyr::left_join(humanPairs, humanConversion, by = "GRCh38_latest_protein")
# Remove col of protein IDs
humanPairs <- dplyr::select(humanPairs, mikado_final_sc2_stringent_noMito_protein, humanGene)
# Remove duplicate rows
humanPairs <- dplyr::distinct(humanPairs)
# Remove NAs or blanks
humanPairs <- humanPairs[!(is.na(humanPairs$humanGene) | humanPairs$humanGene==""), ]
# Collapse by woodchuck
humanPairs <- plyr::ddply(humanPairs, "mikado_final_sc2_stringent_noMito_protein", summarize,
                          humanGene = paste(humanGene, collapse = ";"))
# 18,578 woodchuck genes have human orthologs
write.table(humanPairs, file = "humanPairedOrthos.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# Read in woodchuck and mouse pairwise orthologs from OrthoFinder
mousePairs <- read.table("mikado_final_sc2_stringent_noMito_protein__v__GCF_000001635.27_GRCm39_protein.tsv",
                         header = TRUE,
                         sep = "\t")
# Expand by both mouse and woodchuck
mousePairs <- separate_rows(mousePairs, mikado_final_sc2_stringent_noMito_protein, sep = ", ")
mousePairs <- separate_rows(mousePairs, GCF_000001635.27_GRCm39_protein, sep = ", ")
# Remove decimals from woodchuck only
mousePairs$mikado_final_sc2_stringent_noMito_protein <- gsub("\\.[0-9]+$","",
                                             mousePairs$mikado_final_sc2_stringent_noMito_protein)
# Remove orthogroup column
mousePairs <- dplyr::select(mousePairs,
                            mikado_final_sc2_stringent_noMito_protein,
                            GCF_000001635.27_GRCm39_protein)
# Remove duplicates
mousePairs <- dplyr::distinct(mousePairs)
# Read in conversion file
mouseConversion <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/GRCm39_cdsHeader.tsv",
                              header = FALSE,
                              sep = "\t")
colnames(mouseConversion) <- c("mouseGene", "GCF_000001635.27_GRCm39_protein")
# Combine with pairs
mousePairs <- dplyr::left_join(mousePairs, mouseConversion, by = "GCF_000001635.27_GRCm39_protein")
# Remove col of protein IDs
mousePairs <- dplyr::select(mousePairs, mikado_final_sc2_stringent_noMito_protein, mouseGene)
# Remove duplicate rows
mousePairs <- dplyr::distinct(mousePairs)
# Remove NAs or blanks
mousePairs <- mousePairs[!(is.na(mousePairs$mouseGene) | mousePairs$mouseGene==""), ]
# Collapse by woodchuck
mousePairs <- plyr::ddply(mousePairs, "mikado_final_sc2_stringent_noMito_protein", summarize,
                          mouseGene = paste(mouseGene, collapse = ";"))
# 18,352 woodchuck genes have mouse orthologs
write.table(mousePairs, file = "~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/mousePairedOrthos.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
            
# Read in woodchuck and ybmar pairwise orthologs from OrthoFinder
ybmarPairs <- read.table("mikado_final_sc2_stringent_noMito_protein__v__GCF_003676075.2_GSC_YBM_2.0_protein.tsv",
                         header = TRUE,
                         sep = "\t")
# Expand by both ybmar and woodchuck
ybmarPairs <- separate_rows(ybmarPairs, mikado_final_sc2_stringent_noMito_protein, sep = ", ")
ybmarPairs <- separate_rows(ybmarPairs, GCF_003676075.2_GSC_YBM_2.0_protein, sep = ", ")
# Remove decimals from woodchuck only
ybmarPairs$mikado_final_sc2_stringent_noMito_protein <- gsub("\\.[0-9]+$","",
                                             ybmarPairs$mikado_final_sc2_stringent_noMito_protein)
# Remove orthogroup column
ybmarPairs <- dplyr::select(ybmarPairs,
                            mikado_final_sc2_stringent_noMito_protein,
                            GCF_003676075.2_GSC_YBM_2.0_protein)
# Remove duplicates
ybmarPairs <- dplyr::distinct(ybmarPairs)
# Read in conversion file
ybmarConversion <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/GSC_YBM_2.0_cdsHeader.tsv",
                              header = FALSE,
                              sep = "\t")
colnames(ybmarConversion) <- c("ybmarGene", "GCF_003676075.2_GSC_YBM_2.0_protein")
# Combine with pairs
ybmarPairs <- dplyr::left_join(ybmarPairs, ybmarConversion, by = "GCF_003676075.2_GSC_YBM_2.0_protein")
# Remove col of protein IDs
ybmarPairs <- dplyr::select(ybmarPairs, mikado_final_sc2_stringent_noMito_protein, ybmarGene)
# Remove duplicate rows
ybmarPairs <- dplyr::distinct(ybmarPairs)
# Remove NAs or blanks
ybmarPairs <- ybmarPairs[!(is.na(ybmarPairs$ybmarGene) | ybmarPairs$ybmarGene==""), ]
# Collapse by woodchuck
ybmarPairs <- plyr::ddply(ybmarPairs, "mikado_final_sc2_stringent_noMito_protein", summarize,
                          ybmarGene = paste(ybmarGene, collapse = ";"))
# 19,567 woodchuck genes have ybmar orthologs
write.table(ybmarPairs, file = "~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/ybmarPairedOrthos.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# Read in woodchuck and almar pairwise orthologs from OrthoFinder
almarPairs <- read.table("mikado_final_sc2_stringent_noMito_protein__v__GCF_001458135.1_marMar2.1_protein.tsv",
                         header = TRUE,
                         sep = "\t")
# Expand by both almar and woodchuck
almarPairs <- separate_rows(almarPairs, mikado_final_sc2_stringent_noMito_protein, sep = ", ")
almarPairs <- separate_rows(almarPairs, GCF_001458135.1_marMar2.1_protein, sep = ", ")
# Remove decimals from woodchuck only
almarPairs$mikado_final_sc2_stringent_noMito_protein <- gsub("\\.[0-9]+$","",
                                             almarPairs$mikado_final_sc2_stringent_noMito_protein)
# Remove orthogroup column
almarPairs <- dplyr::select(almarPairs,
                            mikado_final_sc2_stringent_noMito_protein,
                            GCF_001458135.1_marMar2.1_protein)
# Remove duplicates
almarPairs <- dplyr::distinct(almarPairs)
# Read in conversion file
almarConversion <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/marMar2.1_cdsHeader.tsv",
                              header = FALSE,
                              sep = "\t")
colnames(almarConversion) <- c("almarGene", "GCF_001458135.1_marMar2.1_protein")
# Combine with pairs
almarPairs <- dplyr::left_join(almarPairs, almarConversion, by = "GCF_001458135.1_marMar2.1_protein")
# Remove col of protein IDs
almarPairs <- dplyr::select(almarPairs, mikado_final_sc2_stringent_noMito_protein, almarGene)
# Remove duplicate rows
almarPairs <- dplyr::distinct(almarPairs)
# Remove NAs or blanks
almarPairs <- almarPairs[!(is.na(almarPairs$almarGene) | almarPairs$almarGene==""), ]
# Collapse by woodchuck
almarPairs <- plyr::ddply(almarPairs, "mikado_final_sc2_stringent_noMito_protein", summarize,
                          almarGene = paste(almarGene, collapse = ";"))
# 19,067 woodchuck genes have almar orthologs
write.table(almarPairs, file = "~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/almarPairedOrthos.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# Read in woodchuck and squir pairwise orthologs from OrthoFinder
squirPairs <- read.table("mikado_final_sc2_stringent_noMito_protein__v__GCF_000236235.1_SpeTri2.0_protein.tsv",
                         header = TRUE,
                         sep = "\t")
# Expand by both squir and woodchuck
squirPairs <- separate_rows(squirPairs, mikado_final_sc2_stringent_noMito_protein, sep = ", ")
squirPairs <- separate_rows(squirPairs, GCF_000236235.1_SpeTri2.0_protein, sep = ", ")
# Remove decimals from woodchuck only
squirPairs$mikado_final_sc2_stringent_noMito_protein <- gsub("\\.[0-9]+$","",
                                             squirPairs$mikado_final_sc2_stringent_noMito_protein)
# Remove orthogroup column
squirPairs <- dplyr::select(squirPairs,
                            mikado_final_sc2_stringent_noMito_protein,
                            GCF_000236235.1_SpeTri2.0_protein)
# Remove duplicates
squirPairs <- dplyr::distinct(squirPairs)
# Read in conversion file
squirConversion <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/SpeTri2.0_cdsHeader.tsv",
                              header = FALSE,
                              sep = "\t")
colnames(squirConversion) <- c("squirGene", "GCF_000236235.1_SpeTri2.0_protein")
# Combine with pairs
squirPairs <- dplyr::left_join(squirPairs, squirConversion, by = "GCF_000236235.1_SpeTri2.0_protein")
# Remove col of protein IDs
squirPairs <- dplyr::select(squirPairs, mikado_final_sc2_stringent_noMito_protein, squirGene)
# Remove duplicate rows
squirPairs <- dplyr::distinct(squirPairs)
# Remove NAs or blanks
squirPairs <- squirPairs[!(is.na(squirPairs$squirGene) | squirPairs$squirGene==""), ]
# Collapse by woodchuck
squirPairs <- plyr::ddply(squirPairs, "mikado_final_sc2_stringent_noMito_protein", summarize,
                          squirGene = paste(squirGene, collapse = ";"))
# 18,973 woodchuck genes have squir orthologs
write.table(squirPairs, file = "~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/spetriPairedOrthos.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# Now to combine all of this, bring in gene IDs
setwd("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/")
mikado <- read.table("mikado_geneList.txt",
                     header = FALSE,
                     sep = "\t")
# Rename column to match
colnames(mikado) <- "mikado_final_sc2_stringent_noMito_protein"
# Now join everything together
fullOrthos <- dplyr::left_join(mikado, humanPairs, by = "mikado_final_sc2_stringent_noMito_protein")
fullOrthos <- dplyr::left_join(fullOrthos, mousePairs, by = "mikado_final_sc2_stringent_noMito_protein")
fullOrthos <- dplyr::left_join(fullOrthos, ybmarPairs, by = "mikado_final_sc2_stringent_noMito_protein")
fullOrthos <- dplyr::left_join(fullOrthos, almarPairs, by = "mikado_final_sc2_stringent_noMito_protein")
fullOrthos <- dplyr::left_join(fullOrthos, squirPairs, by = "mikado_final_sc2_stringent_noMito_protein")

# Loop through to create hierarchical genes
# Human -> mouse -> yellow -> alpine -> squirrel
# Doing yellow before alpine because more recent annotation project (2019 vs 2015) and
# has more orthologs with woodchuck, better liftoff scores
fullOrthos$hierGene <- rep(NA, nrow(fullOrthos))
for (j in 1:nrow(fullOrthos)) {
  # First: if there is a human gene in table, use it
  if (!is.na(fullOrthos$humanGene[j])) {
    fullOrthos$hierGene[j] <- as.character(fullOrthos$humanGene[j])
  }
  # If no human gene, use mouse gene
  else if (!is.na(fullOrthos$mouseGene[j])) {
    fullOrthos$hierGene[j] <- as.character(fullOrthos$mouseGene[j])
  }
  # If no mouse gene, use yellow-bellied marmot
  else if (!is.na(fullOrthos$ybmarGene[j])) {
    fullOrthos$hierGene[j] <- as.character(fullOrthos$ybmarGene[j])
  }
  # If no yellow-bellied marmot, use Alpine marmot
  else if (!is.na(fullOrthos$almarGene[j])) {
    fullOrthos$hierGene[j] <- as.character(fullOrthos$almarGene[j])
  }
  # If no Alpine marmot, use 13-lined ground squirrel
  else if (!is.na(fullOrthos$squirGene[j])) {
    fullOrthos$hierGene[j] <- as.character(fullOrthos$squirGene[j])
  }
}

# Make hierarchical genes unique
fullOrthos$uniqueHier <- make.unique(fullOrthos$hierGene, sep = "-")
fullOrthos$uniqueHier <- gsub('^NA-[0-9]+$', 'NA', fullOrthos$uniqueHier)

# Get one-to-one orthologs between human and mouse for future conversions

# Isolate mouse and human columns
humanOnes <- dplyr::select(fullOrthos, mikado_final_sc2_stringent_noMito_protein, humanGene)
mouseOnes <- dplyr::select(fullOrthos, mikado_final_sc2_stringent_noMito_protein, mouseGene)
# Count how many have genes
sum(!is.na(humanOnes$humanGene)) # 18,578
sum(!is.na(mouseOnes$mouseGene)) # 18,352
# Sort these to make sure distinct is working properly
humanOnes <- humanOnes[order(humanOnes$humanGene),]
mouseOnes <- mouseOnes[order(mouseOnes$mouseGene),]
# Only keep genes that have no duplicates
humanOnes <- humanOnes[!(duplicated(humanOnes$humanGene)|duplicated(humanOnes$humanGene,
                                                                    fromLast = TRUE)),, drop = FALSE]
dim(humanOnes) # Kept 16,173
mouseOnes <- mouseOnes[!(duplicated(mouseOnes$mouseGene)|duplicated(mouseOnes$mouseGene,
                                                                    fromLast = TRUE)),, drop = FALSE]
dim(mouseOnes) # Kept 16,099
# Rename so distinct from colnames in orthoTable
colnames(humanOnes) <- c("mikado_final_sc2_stringent_noMito_protein", "humanOneToOne")
colnames(mouseOnes) <- c("mikado_final_sc2_stringent_noMito_protein", "mouseOneToOne")
# Add these to table
fullOrthos <- dplyr::left_join(fullOrthos, humanOnes, by = "mikado_final_sc2_stringent_noMito_protein")
fullOrthos <- dplyr::left_join(fullOrthos, mouseOnes, by = "mikado_final_sc2_stringent_noMito_protein")

# Now save this; also available in the Github repo
write.table(fullOrthos, file = "collectedOrthofinderPairings.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
