#package
library(Seurat);cat("Seurat:",as.character(packageVersion("Seurat")),"\n")
library(data.table);cat("data.table:",as.character(packageVersion("data.table")),"\n")
library(tidyr);cat("tidyr:",as.character(packageVersion("tidyr")),"\n")
library(dplyr);cat("dplyr:",as.character(packageVersion("dplyr")),"\n")
library(patchwork);cat("patchwork:",as.character(packageVersion("patchwork")),"\n")
library(ggpubr);cat("ggpubr:",as.character(packageVersion("ggpubr")),"\n")
library(future);cat("future:",as.character(packageVersion("future")),"\n")

#environment
rm(list = ls())
options(stringsAsFactors = F)

#parameter
res.usage.i <- 0.3
thread <- 10
ram <- 10
plan("multisession",workers = thread)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = ram*1024^3)

#workdir
path <- "/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/ileum/DEG/cluster_DEG_top30"
path.res <- paste0(path,'/','res_',res.usage.i)
if (!dir.exists(path.res)) {
  dir.create(path.res)
}
setwd(path.res)

##input RDS
EC <- readRDS("/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/ileum/AM_ileum.RDS")
DefaultAssay(EC) <- 'RNA'

#print parameter
cat('resolution:',mode(res.usage.i),';',res.usage.i,"\n")
cat('out directory:',mode(path.res),';',path.res,"\n")
cat('thread:',mode(thread),';',thread,"\n")
cat('ram:',mode(ram),';',ram,"\n")

#Diff analysis:find marker
EC.assay <- names(EC@assays)
if ('integrated' %in% EC.assay) {
  assay.m <- 'integrated'
}else if ('SCT' %in% EC.assay) {
  assay.m <- 'SCT'
}else{
  assay.m <- 'RNA'
}

res <- paste0(assay.m,'_snn_res.',res.usage.i)
ident <- EC@meta.data[,res]
names(ident) <- row.names(EC@meta.data)
EC@active.ident <- ident
EC <- JoinLayers(EC)
EC.markers <- FindAllMarkers(EC, test.use = "wilcox",only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.1)
markers <- select(EC.markers,c('cluster','gene','p_val_adj','p_val','avg_log2FC'))
write.csv(markers,'./02_diff_marker.csv',row.names = F,quote = F)

#top30_marker_plot
top_gene_df <- data.frame(cluster=character(0),gene=character(0))
top_gene <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC) %>% select(c('cluster','gene','p_val_adj','p_val','avg_log2FC'))
cluster_i <- as.character(unique(top_gene$cluster))
for (i in seq(1,length(cluster_i))) {
  top_gene_i <- filter(top_gene,cluster==cluster_i[i])
  gene_i <- top_gene_i$gene
  gene_i.1 <- paste(gene_i[seq(1,10)],collapse = ',')
  gene_i.2 <- paste(gene_i[seq(11,20)],collapse = ',')
  gene_i.3 <- paste(gene_i[seq(21,30)],collapse = ',')
  gene_i <- paste(gene_i.1,gene_i.2,gene_i.3,sep = '\n')
  top_gene_df[i,1] <- cluster_i[i]
  top_gene_df[i,2] <- gene_i
}

EC@meta.data$celltype <- EC@active.ident
p <- ggtexttable(top_gene_df, rows = NULL)
p
write.csv(top_gene,'./merge_res_0.3 top 30 DEG.csv',row.names = F,quote = F)
