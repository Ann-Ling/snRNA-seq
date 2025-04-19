#Package
library(devtools)
library(monocle)
library(Seurat)
library(future);cat("future:",as.character(packageVersion("future")),"\n")
plan("multisession",workers = 20)
options(future.globals.maxSize = 10*1024^3)
library(ggplot2)

#Environment
rm(list=ls())

#Input
EC<-readRDS("/data2/liuhuiling/gut_snRNA/09pseudotime/ileum/ileum_pseudotime.RDS")
#check
head(EC)
head(EC@active.ident)
EC <- JoinLayers(EC)
mat <-  GetAssayData(object = EC, layer = "counts")
data<-as(as.matrix(mat),'dgCMatrix')
pd<-new("AnnotatedDataFrame",data=EC@meta.data)
fData<-data.frame(gene_short_name=row.names(data),row.names = row.names(data))
fd<-new("AnnotatedDataFrame",data=fData)
mycds <- newCellDataSet(data,
                        phenoData=pd,
                        featureData=fd,
                        expressionFamily=negbinomial.size())
dim(mycds)

mycds<-estimateSizeFactors(mycds)
mycds<-estimateDispersions(mycds,cores=4,relative_expr=TRUE)
disp_table <- dispersionTable(mycds)
View(disp_table)
disp.genes<- subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 *dispersion_fit)$gene_id
mycds<-setOrderingFilter(mycds,disp.genes)
p <- plot_ordering_genes(mycds)
p
pdf("/data2/liuhuiling/gut_snRNA/09pseudotime/ileum/ileum_geneFilter.pdf",width = 5,height=5)
print(p)
dev.off()

#Dimensionality reduction
mycds<-reduceDimension(mycds,max_components = 2,method = 'DDRTree')

#Set ISC/EB as the root node of the trajectory
GM_state <- function(mycds){
  if (length(unique(pData(mycds)$State)) > 1){
    T0_counts <- table(pData(mycds)$State, pData(mycds)$cell_type)[,"ISC/EB"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
mycds <- orderCells(mycds, root_state = GM_state(mycds))

#Save
saveRDS(mycds, file = "/data2/liuhuiling/gut_snRNA/09pseudotime/ileum/ileum_monocle.RDS")

#Visualization
#Plot by Pseudotime
p2 <- plot_cell_trajectory(mycds,color_by="Pseudotime",show_backbone=TRUE)
p2
pdf("/data2/liuhuiling/gut_snRNA/09pseudotime/ileum/ileum_pseudotime_monocle.pdf",width = 10,height=10)
print(p2)
dev.off()

#Plot by cell type
colour <- c("ISC/EB" = "#248067","EE" = "#61649f","EC" = "#f28e16")
p3 <- plot_cell_trajectory(mycds,color_by="sub_cell_type",show_backbone=TRUE)+
  scale_color_manual(values=colour)+
  theme(legend.position = "right")
p3
pdf("/data2/liuhuiling/gut_snRNA/09pseudotime/ileum/ileum_celltype_monocle.pdf",width = 10,height=8)
print(p3)
dev.off()
