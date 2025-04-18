# Introduction
To elucidate how Relish influences transcriptional programs in ISC/EBs, we constructed a SCENIC database for *A. mellifera* and predicted GRN based on coexpression and motif enrichment. The key steps include:
1. Creating cisTarget databases
2. Formatting input
3. SCENIC

# Pipeline overview
## Step 1: Creating cisTarget databases
Since single-cell regulatory network inference and clustering (SCENIC) does not provide pre-built motif databases for *A. mellifera*, we first constructed the motif databases for cisTarget and SCENIC analyses specific to this species.
To support this analysis, we need to prepare the following documents:
1. The genome annotation file (.gtf) of *A. mellifera*
2. Motifs in Cluster-Buster format
3. File with motif IDs

## Step2: Formatting input
This step is used to extract data from the .RDS file and convert it into .loom file format

## Step3: SCENIC
To successfully complete this step：
1. You need to install pySCENIC in advance. Refer to https://pyscenic.readthedocs.io/en/latest/index.html
2.  According to your own data, follow the tutorial ([pySCENIC](https://pyscenic.readthedocs.io/en/latest/index.html)) or run `03SCENIC.R`

