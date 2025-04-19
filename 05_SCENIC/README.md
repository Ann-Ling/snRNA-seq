# Introduction
To elucidate how Relish influences transcriptional programs in ISC/EBs, we constructed a SCENIC database for *A. mellifera* and predicted GRN based on coexpression and motif enrichment.
The analysis involved three key steps:
1. Creating cisTarget databases
2. Formatting input
3. SCENIC
4. Integrate SCENIC analysis results into the `.RDS` file before SCENIC analysis

# Pipeline overview
## Step 1: Creating cisTarget databases
Since single-cell regulatory network inference and clustering (SCENIC) does not provide pre-built motif databases for *A. mellifera*, we first constructed the motif databases for cisTarget and SCENIC analyses specific to this species.
To support this process, the following files were prepared:
1. The genome annotation file `.gtf` of *A. mellifera*
2. Motifs in Cluster-Buster format
3. A file containing the corresponding motif IDs

## Step2: Formatting input
Extracts data from a `.RDS` object and exports it as a `.loom` file for compatibility with downstream tools.

## Step3: SCENIC
To complete this step：
1. Ensure that pySCENIC is installed in advance. For installation and usage, refer to the official pySCENIC tutorial [pySCENIC](https://pyscenic.readthedocs.io/en/latest/index.html).
2. Depending on your data and workflow preference, you may either follow the standard pySCENIC pipeline [pySCENIC](https://pyscenic.readthedocs.io/en/latest/index.html) or execute the analysis using `03SCENIC.sh`.

## Step4: Integrate results into the `.RDS` file
In this step, the inferred activity of each transcription factor is incorporated into the cell-level metadata to facilitate further visualization and interpretation.
