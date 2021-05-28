# Original script made Jan 15, 2021
# Modified to suit SC2 May 28, 2021

library(dplyr)
library(tidyverse)

setwd("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/protein_seqs/OrthoFinder/Results_May24/Orthologues/Orthologues_mikado_final_sc2_stringent_noMito_protein/")
# Read in woodchuck and human pairwise orthologs from OrthoFinder
orthoFinderPairs <- read.table("mikado_final_sc2_stringent_noMito_protein__v__GRCh38_latest_protein.tsv",
                               header = TRUE,
                               sep = "\t")
# Use rows as homologroups
humanGroups <- dplyr::select(orthoFinderPairs, mikado_final_sc2_stringent_noMito_protein, GRCh38_latest_protein)
humanGroups$homoloGroup <- row.names(humanGroups)
# Expand by both human and woodchuck
humanGroups <- separate_rows(humanGroups, GRCh38_latest_protein, sep = ", ")
humanGroups <- separate_rows(humanGroups, mikado_final_sc2_stringent_noMito_protein, sep = ", ")
# Remove decimals from woodchuck only
humanGroups$mikado_final_sc2_stringent_noMito_protein <- gsub("\\.[0-9]+$","",
                                                  humanGroups$mikado_final_sc2_stringent_noMito_protein)
# Remove duplicate rows
humanGroups <- dplyr::distinct(humanGroups)
# Get entrez id for human genes
# Read in conversion file
humanConversion <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/GRCh38.p13_cdsHeader_ID.tsv",
                              header = FALSE,
                              sep = "\t",
                              na.strings=c("", "NA"))
humanConversion[rowSums(is.na(humanConversion)) > 0, ] # Check for NAs
colnames(humanConversion) <- c("humanGene", "entrezID", "GRCh38_latest_protein")
# Separate human and woodchuck info
woodchuckInfo <- dplyr::select(humanGroups, homoloGroup, mikado_final_sc2_stringent_noMito_protein)
humanInfo <- dplyr::select(humanGroups, homoloGroup, GRCh38_latest_protein)
# Add entrez ID and such to human info
humanInfo <- dplyr::left_join(humanInfo, humanConversion, by = "GRCh38_latest_protein")
# Only keep desired columns
humanInfo <- dplyr::select(humanInfo, homoloGroup, humanGene, entrezID)
# Remove duplicates
humanInfo <- dplyr::distinct(humanInfo)
# Add species ID
humanInfo$taxID <- rep(9606, nrow(humanInfo))
# Now get rid of woodchuck duplicates
woodchuckInfo <- dplyr::distinct(woodchuckInfo)
# Add speciesID
woodchuckInfo$taxID <- rep(9995, nrow(woodchuckInfo))
# Do a full join by row
allGroups <- dplyr::inner_join(humanInfo, woodchuckInfo, by = "homoloGroup")
# Remove any row that has a single NA
allGroups[rowSums(is.na(allGroups)) > 0,] # There are a few hundred NAs
allGroups <- allGroups[complete.cases(allGroups), ]
# Remove all rows that have a repeat woodchuck gene ID or human entrezID
allDistinct <- dplyr::distinct(allGroups, entrezID, mikado_final_sc2_stringent_noMito_protein, .keep_all = TRUE)
# Now separate back into human and woodchuck
humanDistinct <- dplyr::select(allDistinct, homoloGroup, taxID.x,
                               entrezID, humanGene)
woodchuckDistinct <- dplyr::select(allDistinct, homoloGroup,
                                   taxID.y, mikado_final_sc2_stringent_noMito_protein)
# Remove duplicates
humanDistinct <- dplyr::distinct(humanDistinct)
woodchuckDistinct <- dplyr::distinct(woodchuckDistinct)
# Read in woodchuck file with unique ortholog names
woodchuckNames <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/collectedOrthofinderPairings.tsv",
                             sep = "\t",
                             header = TRUE)
woodchuckNames <- dplyr::select(woodchuckNames, mikado_final_sc2_stringent_noMito_protein, uniqueHier)
# Add names to woodchuck table
woodchuckDistinct <- dplyr::left_join(woodchuckDistinct, woodchuckNames,
                                      by = "mikado_final_sc2_stringent_noMito_protein")
# Check for NAs just in case
woodchuckDistinct[rowSums(is.na(woodchuckDistinct)) > 0,] # Nope!
# Rename human and woodchuck so they will stack
colnames(humanDistinct) <- c("homoloGroup", "taxID", "geneID", "geneSymbol")
colnames(woodchuckDistinct) <- c("homoloGroup", "taxID", "geneID", "geneSymbol")
# Join by column names
allFinal <- dplyr::bind_rows(humanDistinct, woodchuckDistinct)
# Write table
write.table(allFinal, file= "~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/homologroupHumanWoodchuck.tsv",
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")
