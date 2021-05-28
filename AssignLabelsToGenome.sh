#!/bin/bash
# This script takes the output of getOrthologLabels.R (collectedOrthoFinderPairings.tsv) and uses this 
# to label the woodchuck annotation GFF file with gene symbols

# First, the GFF annotation file will be converted to a simpler, standardized GTF file with Gffread v0.12.1
gffread -T mikado_final_sc2_stringent.gff -o mikado_final_sc2_stringent.gffread.gtf

# Then, collectedOrthoFinderPairings.tsv needs to be modified into a file that can be used to add
# gene symbols to the appropriate gene labels in the GTF file. This will label the genes with
# hierarchically assigned, unique gene symbols (e.g. a gene with no human or mouse gene will be
# labeled with a yellow-bellied marmot gene if it's available, and these genes have unique extensions
# to ensure the same symbol does not label two different genes).
# Get rid of header
tail -n +2 collectedOrthofinderPairings.tsv > robustOrthoNames.tsv
# Add backslashes to things that also are symbols with meaning for something else
# (dashes and periods)
sed -i 's/\-/\\-/g' robustOrthoNames.tsv
sed -i 's/\./\\./g' robustOrthoNames.tsv
# Only keep hier names and gene names
cut -f 1,8 robustOrthoNames.tsv > robustHierOrthos.tsv
# Remove NAs
grep -v -P '\tNA$' robustHierOrthos.tsv > test && mv test robustHierOrthos.tsv

# Make a copy of the GTF file to be modified with the addition of the orthologous gene symbols
cp mikado_final_sc2_stringent.gffread.gtf mikado_sc2_orthoNames.gtf

# Loop through robustHierOrthos.tsv to add a gene_name attribute to genes that have the appropriate
# gene ID. This will likely take more than a day to run.
while IFS=$'\t' read -r mikado ortho
do
     echo "Replacing" $mikado "with" $ortho
     # The following line is the one that adds the gene_name attribute after the appropriate gene_id
     sed -i -e "s/gene_id \"${mikado}\";/& gene_name \"${ortho}\";/" mikado_sc2_orthoNames.gtf
done < robustHierOrthos.tsv

# At this point, the GTF file still doesn't have the mitochondrial genes in it. This exists in
# a separate file and will be concatenated along with the full genome annotation.
