**Tips**:
1. `01PrepareSce.R`: this script converts the Seurat object into a SingleCellExperiment object, enabling compatibility with downstream statistical analyses.
2. `02CellRatioShifts.R`: this script calculates the frequency of each cell type across samples. Differential abundance analysis between groups is performed using edgeR, primarily through the glmFit() and glmLRT() functions.
3. `03Visualization.R`: this script generates visualizations of cell type proportions across treatment groups or intestinal compartments.
4. These scripts are also applicable for the analysis and comparison of EC and EE subtypes.
