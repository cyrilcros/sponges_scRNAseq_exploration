# Profiling cellular diversity in sponges informs animal cell type and nervous system evolution

## *Spongilla* proteome/transcriptome

This directory contains proteome and transcriptome fasta files, and a gene lookup table, used in the study Musser et al. 2019,  "Profiling cellular diversity in sponges informs animal cell type and nervous system evolution."

The following files contain:

1) `spongilla_transcriptome_annotated_200501.fasta` - *Spongilla* transcriptome. Fasta header is "trinity transcriptID|automated gene name|manually curated gene name"
2) `spongilla_proteome_100AA_annotated_210501.fasta`- *Spongilla* proteome. Fasta header is "transdecoder proteinID|automated gene name|manually curated gene name"
3) `spongilla_merged_genes_and_names.tsv` - lookup table with trinity gene ID, automated gene name, and manually-curated gene name. Trinity gene IDs were merged into a final gene ID using the *Spongilla* phylome, which identified Trinity genes that were in fact likely 5' and 3' ends of the same gene. The automated gene name thus also represents a merged gene identifier, which is found in expression matrix and Seurat object. Note, merged gene IDs were only used for building the expression count matrix, and we did not attempt to merge spongilla transcript or protein sequences into a single representative sequence. Manually-curated names for each gene (of which there are about 600), were used in the manuscript.
4) `spongilla_trees.zip` - zip archive of directory containing all gene trees from the spongilla phylome.
5) `spongilla_orthology_table.2.1_nobrackets.tsv` - orthology table generated via parsing phylome gene trees.