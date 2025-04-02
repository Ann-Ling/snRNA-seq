#package
library(Seurat);cat("Seurat:",as.character(packageVersion("Seurat")),"\n")
library(clustree);cat("clustree:",as.character(packageVersion("clustree")),"\n")
library(data.table);cat("data.table:",as.character(packageVersion("data.table")),"\n")
library(future);cat("future:",as.character(packageVersion("future")),"\n")

#environment
rm(list = ls())
options(stringsAsFactors = F)

#multiprocess
thread <- 15
ram <- 20
plan("multisession",workers = thread)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = ram*1024^3)

##Get the parameters
#parser = argparse::ArgumentParser(description="script to Cluster scRNA data")
#parser$add_argument('-d','--dim',help='dim usage')
#parser$add_argument('-k','--knn',help='defines k for the k-nearest neighbor algorithm')
#parser$add_argument('-r','--res_file',help='resolution usage')
#parser$add_argument('-i','--input', help='input merge seurat object rds')
#parser$add_argument('-o','--out',help='out directory')
#parser$add_argument('-t','--thread', help='thread used')
#parser$add_argument('-g','--ram', help='ram used(Gb)')
#args = parser$parse_args()

#print parameter
cat('dim usage:',mode(dim.usage),';',dim.usage,"\n")
cat('k for the k-nearest neighbor algorithm:',mode(k.usage),';',k.usage,"\n")
cat('resolution:',mode(res.usage.seq),';',res.usage.seq,"\n")
cat('input seurat rds:',mode(input),';',input,"\n")
cat('out directory:',mode(path),';',path,"\n")
cat('thread:',mode(thread),';',thread,"\n")
cat('ram:',mode(ram),';',ram,"\n")

#dim and cluster
dim.usage <- 30
k.usage <- 20
res.usage.seq <-as.vector(as.matrix(read.table("/data2/liuhuiling/gut_snRNA/03cluster_resolution/code/resolution.txt",
                                               header=F,stringsAsFactors=F)))
##data input
input <- "/data2/liuhuiling/gut_snRNA/02merge_RDS/merge/merge_final.RDS"
#workdir
path <- "/data2/liuhuiling/gut_snRNA/03cluster_resolution/merge"
if (!dir.exists(path)) {
  dir.create(path)
}
setwd(path)
 
#read data
EC <- readRDS(input)
EC.assay <- names(EC@assays)
if ('integrated' %in% EC.assay) {
  assay.m <- 'integrated'
}else if ('SCT' %in% EC.assay) {
  assay.m <- 'SCT'
}else{
  assay.m <- 'RNA'
}
DefaultAssay(EC) <- assay.m

#standard analysis
EC <- RunTSNE(EC, dims = 1:dim.usage)
#cluster
EC <- FindNeighbors(EC, reduction = "pca", dims = 1:dim.usage, k.param = k.usage, annoy.metric = 'euclidean')
EC <- FindClusters(EC, resolution = res.usage.seq, algorithm = 1, random.seed = 0)

#RNA assay all genes scale
if ('SCT' %in% EC.assay) {
  assay <- 'RNA'
}else{
  assay <- 'RNA'
}
DefaultAssay(EC) <- assay
all.genes <- rownames(EC)
EC <- ScaleData(EC, features = all.genes)
a <- row.names(EC@meta.data)
b <- sub("^[^_]*_([^_]*)_.*$", "\\1", a)
EC@meta.data$region <- b
head(EC)
saveRDS(EC,paste(path,'/',"01_cluster.RDS",sep=""))

#clustree_plot
p <- clustree(EC@meta.data,prefix = paste0(assay.m,'_snn_res.'),node_size = 10)
p
pdf(paste(path,"/01_clustree.pdf",sep=""),
    height = length(res.usage.seq)*3,width = length(unique(EC@active.ident))*5)
print(p)
dev.off()

