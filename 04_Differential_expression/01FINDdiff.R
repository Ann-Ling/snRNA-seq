#Package
library(Seurat);cat("Seurat:",as.character(packageVersion("Seurat")),"\n")
library(tidyverse)
library(harmony)
library(ggrepel)
library(patchwork)
library(future);cat("future:",as.character(packageVersion("future")),"\n")
plan("multisession",workers = 15)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 20*1024^3)

#Environment
rm(list=ls())

#input
testseu=readRDS("/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/AM.RDS")
testseu <- JoinLayers(testseu)
head(testseu)

#Find different expression genes
testseu@meta.data$celltype_condition=paste(testseu@meta.data$cell_type,testseu@meta.data$orig.ident,sep = "_")

marker_condition=data.frame()
Idents(testseu)="celltype_condition"
for ( ci in sort(as.character(unique(testseu@meta.data$cell_type))) ) {
  tmp.marker <- FindMarkers(
    testseu, logfc.threshold = 0, min.pct = 0.1,
    only.pos = F, test.use = "wilcox",
    ident.1=paste0(ci,"_CV"),ident.2=paste0(ci,"_GF")
  )
  
  tmp.marker$gene=rownames(tmp.marker)
  tmp.marker$condition=ifelse(tmp.marker$avg_log2FC > 0,paste0(ci,"_CV"),paste0(ci,"_GF"))
  tmp.marker$cluster=ci
  
  #tmp.marker=tmp.marker%>%filter(p_val_adj < 0.01)
  tmp.marker=as.data.frame(tmp.marker)
  tmp.marker=tmp.marker%>%arrange(desc(avg_log2FC))
  
  marker_condition=marker_condition%>%rbind(tmp.marker)
}
path <-'/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/DEG'
setwd(path)
write.table(marker_condition,file = "markers.BasedOncondition.txt",quote = F,sep = "\t",
            row.names = F,col.names = T)

