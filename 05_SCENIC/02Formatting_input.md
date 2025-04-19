# Data extraction
## SCENIC_annotation.txt
```txt
ISC/EB
EE
EC
VM
```
## Data extraction
```R
rm(list=ls())
library(Matrix)
#library(Matrix.utils)
library(dplyr)
library(data.table)
library(Seurat)
library(future);cat("future:",as.character(packageVersion("future")),"\n")
plan("multisession",workers = 20)
options(future.globals.maxSize = 10*1024^3)

#Split Parameters
raw <-"/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/ileum"
column <- "cell_type"
rds.file <- paste0(raw,'/AM_ileum.RDS')
EC <- readRDS(rds.file)
table(Idents(EC))

anno.file <- '/data2/liuhuiling/gut_snRNA/10SCENIC/02ileum/SCENIC_annotation.txt'
anno.df <- read.delim(anno.file,header = F)
#select aim cell, get check
cell_type <- as.character(anno.df$V1)#For data extraction without N
check <- as.numeric(sapply(EC@meta.data[,column], function(x){x %in% c(cell_type)}))
EC@meta.data$check <- factor(check)


ident <- EC@meta.data$check
names(ident) <- row.names(EC@meta.data)
EC@active.ident <- ident

#Subset
EC.split <- subset(x = EC, idents = 1)
saveRDS(EC.split,'/data2/liuhuiling/gut_snRNA/10SCENIC/02ileum/ileum_scenic.RDS')
write.table(EC.split@meta.data,'/data2/liuhuiling/gut_snRNA/10SCENIC/02ileum/ileum_scenic.tsv',
            row.names = F,quote = F,sep = '\t')
##counts
EC <- JoinLayers(EC)
counts <- GetAssayData(object = EC,assay = "RNA",layer  = "counts")
write.csv(t(as.matrix(counts)),'/data2/liuhuiling/gut_snRNA/10SCENIC/02ileum/ileum_scenic_counts_t.csv',row.names = T)
```

# loom file conversion
```python
import os
import loompy as lp
import numpy as np
import scanpy as sc

# Check the current working directory
print("Current Working Directory:", os.getcwd())
print("Files in Directory:", os.listdir(os.getcwd()))

# Reading CSV Files
x = sc.read_csv("ileum_scenic_counts_t.csv")  

# 创建 row_attrs 和 col_attrs
row_attrs = {"Gene": np.array(x.var_names)}  # Genes as rows
col_attrs = {"CellID": np.array(x.obs_names)}  # Cells as columns

# Transpose and create loom file
lp.create("ileum_scenic.loom", x.X.transpose(), row_attrs, col_attrs)

#Check data correctness
import loompy
loom_file = loompy.connect('ileum_scenic.loom')
genes = loom_file.ra.Gene  
cells = loom_file.ca.CellID  
loom_file.close()
print("Genes in Loom file:", genes)
```
