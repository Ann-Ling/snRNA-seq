# Introduction
This folder contains scripts for evaluating transcriptional similarities between cell clusters across species using MetaNeighbor. The analysis computes AUROC (Area Under the Receiver Operating Characteristic) scores to identify transcriptionally similar clusters between species. The key steps include:
1. Homology comparison
2. Data preprocessing
3. Pseudo-cell generation
4. Cross-species correlation analysis and Visualization

# Pipeline overview
## Step 1: Homology comparison (01Diamond.sh)
To support our analysis, FASTA-formatted protein sequences for *Apis mellifera* and *Drosophila melanogaster* were obtained from the NCBI database.
We will get homologous genes in *Drosophila melanogaster* and *Apis mellifera*

## Step 2: Data preprocessing (02Data preprocessing.md)
To support our analysis, we performed the following steps:
1. Downloaded the genome annotation file (.gff) from NCBI and extracted the corresponding relationships between protein IDs and gene IDs based on the annotation information.
2. Extracted cell ordering, metadata, and gene expression count matrices from the Seurat .RDS file to enable downstream integration and analysis.

We will get:
1. gene expression count matrices
2. metadata
3. cell ordering

## Step 3: Pseudo-cell generation (03Pseudo-cell generation.sh)
This step requires one additional scripts: `random.combine.pl`

We will get:
1. Generate pseudo-cell UMI matrices by summing the UMI counts of ten randomly selected cells per cluster.
2. Construct metadata for pseudo-cells.
3. Identify orthologous genes between species to create a common pseudo-cell UMI matrix.

## Step 4: Cross-species correlation analysis and Visualization (04Auroc.analyais.R)
To support our analysis, we performed the following steps:
1. Perform z-score normalization.
2. Compute AUROC scores to quantify transcriptional similarity.
3. Generate a heatmap for visualizing.

We used homologous genes in *Drosophila melanogaster* and *Apis mellifera* for `marker.list`

This step requires two additional scripts for AUROC calculations:
1. `2017-08-28-runMN-US.R`
2. `2017-08-28-runMN-US.pearson.R`
