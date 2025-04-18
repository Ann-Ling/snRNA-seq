#Package
library("clusterProfiler")
library('grid')
library('readr')
vplayout <- function(x,y){
  viewport(layout.pos.row=x,layout.pos.col=y)}

#Environment
rm(list=ls())
path <- "/data2/liuhuiling/gut_snRNA/08KEGG/code"
setwd(path)

#input
#KEGG annotation file
name_term_gene <- read_csv("./name_term_gene_3.csv")
term2name <- name_term_gene[,c(2,1)]
term2gene <- name_term_gene[,c(2,3)]
gene1 <- read.csv('/data2/liuhuiling/gut_snRNA/04DEG_find/merge/res_0.3/DEG.csv')
#There is a 'cluster' column in gene1 that represents cell type information
gene1$cluster <- as.character(gene1$cluster)
unique_groups <- unique(gene1$cluster)
all_results <- data.frame()
for (group in unique_groups) {
  genes <- gene1$gene[gene1$cluster == group]
  print(paste("Analyzing group:", group))
  print(genes)
  
  # Enrichment
  enrichment_result <- enricher(genes, TERM2GENE = term2gene, TERM2NAME = term2name,
                                pvalueCutoff = 1, pAdjustMethod = 'BH', qvalueCutoff = 0.5,
                                minGSSize = 0, maxGSSize = 10000)
  
  if (!is.null(enrichment_result) && nrow(as.data.frame(enrichment_result)) > 0) {
    enrichment_df <- as.data.frame(enrichment_result)
    enrichment_df$cluster <- group
    all_results <- rbind(all_results, enrichment_df)
  } else {
    print(paste("No enrichment results for group:", group))
  }
}

#Output
write.csv(all_results, "/data2/liuhuiling/gut_snRNA/08KEGG/merge/res_0.3/DEG_KEGG.csv", 
          row.names = FALSE, quote = FALSE)
