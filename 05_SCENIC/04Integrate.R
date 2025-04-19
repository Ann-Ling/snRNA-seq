#library
library(Seurat)  
library(SCopeLoomR)  
library(AUCell)  
library(SCENIC)  
library(dplyr)  
library(KernSmooth)  
library(RColorBrewer)  
library(plotly)  
library(BiocParallel)  
library(grid)  
library(ComplexHeatmap)  
library(data.table)  
library(patchwork)  
library(ggplot2)  
library(stringr)  
library(circlize)  

#input
setwd("/data2/liuhuiling/gut_snRNA/10SCENIC/02ileum") 
setwd("/data2/liuhuiling/gut_snRNA/10SCENIC/02ileum")  

#load data  
loom <- open_loom('aucell.loom')  
regulons_incidMat <- SCopeLoomR::get_regulons(loom, column.attr.name="Regulons")  
regulons_incidMat[1:4,1:4]
# LOC408891 LOC410725 LOC725960 LOC409036
#Antp(+)         0         0         0         0
#Dl(+)           0         0         0         0
#Dll(+)          0         0         0         0
#E74(+)          0         0         0         0

regulons <- SCENIC::regulonsToGeneLists(regulons_incidMat)  
#regulunsToGeneLists function returns a list of target genes regulated by each transcription factor
head(regulons,3)
regulonAUC <- SCopeLoomR::get_regulons_AUC(loom,column.attr.name='RegulonsAUC')  
regulonAUC
regulonAucThresholds <- SCopeLoomR::get_regulon_thresholds(loom)  
head(regulonAucThresholds,3)
embeddings <- SCopeLoomR::get_embeddings(loom)  
embeddings
close_loom(loom)

#Transcription factor activity scores were integrated into the cell metadata.
#precoess  
sub_regulonAUC <- regulonAUC[,match(colnames(EC),colnames(regulonAUC))] 
# check  
identical(colnames(sub_regulonAUC), colnames(EC))  
# combine  
auc_data <- getAUC(sub_regulonAUC)
auc_data <- auc_data[, match(colnames(EC), colnames(auc_data))]
EC@meta.data <- cbind(EC@meta.data, t(auc_data))
#check
head(EC)

#Save
saveRDS(EC,"ileum_scenic_final.RDS")
