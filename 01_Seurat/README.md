Tips
1. For **02merge_RDS.R**, you should create a **merge_rds_list.txt**, which contains the addresses of all files output by **01RDS.R**
2. For **03cluster_tree.R**, you should create a **resolution.txt**, which contains different resolution parameter (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.2), followed by clustree analysis.
3. For **04DEG_find.R** Markers was defined utilizing FindAllMarkers function in this script. Significantly overexpressed genes were defined based on the Wilcoxon rank-sum test with a P value < 0.01 and log2 FC > 0.1
