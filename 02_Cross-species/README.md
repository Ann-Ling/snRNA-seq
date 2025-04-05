# Introduction
This folder contains scripts for evaluating transcriptional similarities between cell clusters across species using MetaNeighbor. The analysis computes AUROC (Area Under the Receiver Operating Characteristic) scores to identify transcriptionally similar clusters between species. The key steps include:
1. Homology comparison
2. Data preprocessing
3. Pseudo-cell generation
4. Cross-species correlation analysis
5. Visualization

# Pipeline overview
## Step 1: Homology comparison (01Diamond.sh)
To support our analysis, FASTA-formatted protein sequences for *Apis mellifera* and *Drosophila melanogaster* were obtained from the NCBI database.

## Step 2: Data preprocessing (02Data preprocessing.md)
To support our analysis, we performed the following steps:
1. Downloaded the genome annotation file (.gff) from NCBI and extracted the corresponding relationships between protein IDs and gene IDs based on the annotation information.
2. Extracted cell ordering, metadata, and gene expression count matrices from the Seurat .RDS file to enable downstream integration and analysis.
