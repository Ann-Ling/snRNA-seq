**Tips**:
1. For `01PrepareSce.R`, this script converts the Seurat object into a SingleCellExperiment object, enabling compatibility with downstream statistical analyses.
2. For `02CellRatioShifts.R`, this script calculates the frequency of each cell type across samples. Differential abundance analysis between groups is performed using edgeR, primarily through the glmFit() and glmLRT() functions.
3. For `03Visualization.R`, This script generates visualizations of cell type proportions across treatment groups or intestinal compartments.
4. These scripts are equally applicable for analyzing and comparing EC or EE subtypes.
