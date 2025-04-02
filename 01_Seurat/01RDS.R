#gzip -d barcodes.tsv.gz
#gzip -d features.tsv.gz
#gzip -d matrix.mtx.gz
#package
library(future);cat("future:",as.character(packageVersion("future")),"\n")
plan("multisession",workers = 10)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 10*1024^3)
### package
library(Seurat);cat("Seurat:",as.character(packageVersion("Seurat")),"\n")
rm(list=ls())

#Raw data address
input='/data2/liuhuiling/gut_snRNA/00Raw_data/rectum/GF_rectum/GF_rectum_04' 
#Output dddress
out="/data2/liuhuiling/gut_snRNA/01RDS/merge/GF"
batch='GF_rectum_04'

#Reading Data
input=input
EC.data <- Read10X(data.dir = input,gene.column = 1)
EC.data <- EC.data[,-1]
batch <- batch
cell_name <- paste0(batch,'_',seq(1,ncol(EC.data)))
colnames(EC.data) <- cell_name
EC <- CreateSeuratObject(EC.data,project = batch)
rm(EC.data)

#Standard Process
EC <- NormalizeData(EC)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000)
EC <- ScaleData(EC)
EC <- RunPCA(EC)
dim.usage <- 30
EC <- RunUMAP(EC, dims = 1:dim.usage) 

#Output File
out=out
saveRDS(EC,paste(out,"/",batch,"_QCobject.RDS",sep=""))

