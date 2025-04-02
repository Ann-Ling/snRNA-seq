#package
library(future);cat("future:",as.character(packageVersion("future")),"\n")
plan("multisession",workers = 10)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 20*1024^3)
library(Seurat);cat("Seurat:",as.character(packageVersion("Seurat")),"\n")
library(tidyr);cat("tidyr:",as.character(packageVersion("tidyr")),"\n")
rm(list=ls())

#Get the parameters
dim.usage <- 30
pc.usage <- 50
seed.usage <- 0
maxdim.usage <- "2L"

#FindIntegrationAnchors
reduction.method <- 'rpca' #c("cca", "rpca", "rlsi")
norm.method <- 'LogNormalize' #c("LogNormalize", "SCT")

#Input RDS files
out <- '/data2/liuhuiling/gut_snRNA/02merge_RDS/merge'
#Output file
input <- '/data2/liuhuiling/gut_snRNA/02merge_RDS/code/merge_rds_list.txt'
files <- as.vector(as.matrix(read.table(input,header=F,stringsAsFactors=F)))

#Read data
objectlist <- list()
for(i in 1:length(files)){
  objectlist[[i]] <- readRDS(files[i])
}
cat("The number of RDS inputed:",length(files),"\n")

print('Objectlist creat begin!')
objectlist <- lapply(X = objectlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x)
  x <- RunPCA(x)
  x <- RunUMAP(x,dims = 1:dim.usage)
})
saveRDS(objectlist,paste0(out,"/objectlist.RDS"))
print('Objectlist done!')
cat('\n')

### Seurat Integrate
print('Integrate begin!')
features <- SelectIntegrationFeatures(object.list = objectlist)
#merge the data
object.anchors <- FindIntegrationAnchors(object.list = objectlist,
                                         dims = 1:dim.usage,
                                         reduction=reduction.method,
                                         normalization.method=norm.method)
object.combined <- IntegrateData(anchorset = object.anchors, dims = 1:dim.usage)
saveRDS(object.combined,paste0(out,"/merge.RDS"))
print('Integrate done!')
cat('\n')

#Standardization, normalization, dimensionality reduction
print('Standard analysis begin!')
DefaultAssay(object.combined) <- "integrated"
object.combined <- ScaleData(object.combined, verbose = FALSE)
object.combined <- RunPCA(object.combined, npcs = pc.usage, seed.use = seed.usage,verbose = FALSE)
object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:dim.usage, seed.use = seed.usage, max.dim = maxdim.usage)

#Save the file
saveRDS(object.combined,paste0(out,"/merge_final.RDS"))
print('Standard analysis done!')

