#package
rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)

#Input
EC <- readRDS("/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/01_cluster_anno.RDS")
EC <- JoinLayers(EC)
head(EC)

#Classical marker genes
selected_ISCEB <- c("LOC406097", "LOC100578939","LOC100577397","LOC726514","LOC727284")
selected_EC <- c("LOC552299","LOC412810","LOC724308","LOC410530")
selected_EE <- c("LOC408960","brp","LOC406073")
selected_VM <- c("LOC725878","LOC409805","LOC409881","LOC725404")
selected_EN <- c("LOC726184","LOC552792","LOC551571","LOCPsq")
selected_others <- c("LOC412976")

EC$integrated_snn_res.0.3 <- factor(EC$integrated_snn_res.0.3,
                                    levels = c("10","14","18","19","20",
                                               "5","22",
                                               "0","2","3","6","7","9","11","17","21",
                                               "8","12","28","29",
                                               "16",
                                               "1","4","13","15","23","24","25","26","27"))

p <- DotPlot(EC,features = c(selected_ISCEB, selected_EE, selected_EC, selected_VM, selected_EN, selected_others),
             group.by = 'integrated_snn_res.0.3') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 45),
        plot.margin = margin(t = 10, r = 10, b = 50, l = 10, unit = "pt")) + 
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3))+ 
  coord_flip()
p

#Save the plot to a PDF file
pdf('/data2/liuhuiling/gut_snRNA/05cluster_marker_bar_umi/merge/res_0.3/marker.pdf',
    width = 10, height = 6)
print(p)
dev.off()




