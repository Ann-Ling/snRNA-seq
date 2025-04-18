**Tips**
1. In `01FINDdiff.R`, the `FindMarkers` function is used to compute cell type-specific differentially expressed genes (DEGs) betwwen MF and CV groups.
2. In `02Volcano.R`, DEGs between MF and CV groups are visualized using the **EnhancedVolcano** and **ggplot2** packages.
3. In `03KEGG.R`
     * KEGG annotation file should be determined in advance based on your specific data, and the code should be adapted accordingly.
     * Genes listed in `DEG.CSV` met the criteria of |log2 fold change| >1 and p < 0.05.
