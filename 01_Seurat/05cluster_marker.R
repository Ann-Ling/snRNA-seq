rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
path <-"/data2/liuhuiling/gut_snRNA/05cluster_marker_bar_umi/merge/res_0.3"
EC <- readRDS("/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/cell_type_RDS/EE/01_cluster_anno.RDS")
EC <- JoinLayers(EC)
head(EC)

#cluster在纵坐标(经典marker)
selected_ISCEB <- c("LOC410351","LOC725617","LOC726514","LOC100577397","LOC406097","LOC727284")
selected_EE <- c("Apime-ASTA","brp","LOC100576163","LOC406073","LOC408960")
selected_EC <- c("LOC552299","LOC551687","LOC724259","LOC724308","LOC406114")
selected_VM <- c("LOC413946","LOC725878","LOC725404","LOC409805")
selected_others <- c("LOC412976")

EC_submarker <- c("LOC724308","LOC406114","LOC412810","LOC724259","LOC552299")

EE_submarker <- c("LOC113218809","Apime-ASTA","LOC406073","brp","LOC408960")
#cluster在纵坐标(所有marker)
#selected_ISCEB <- c('arm','LOC100577348',"LOC409650","LOC410351",'LOC552697',"LOC552742","LOC725617","LOC726514",
#                    "LOC413261","LOC413870","LOC100577692","LOC102654776","LOC100577397","LOC100578939","LOC107964400",
#                    "LOC406097","LOC410464","LOC413742","LOC725336","LOC409027","LOC412092")
#selected_EE <- c("LOC726150","LOC413289","LOC409822","LOC727284","Apime-ASTA","LOC100577681","brp","LOC100576163",
#                 "LOC113218809","LOC406073","LOC412839","LOC408960","LOC551717","LOC725573","LOC724270","LOC724113",
#                 "Burs","Tk")
#selected_EC <- c("LOC552299","LOC727092","LOC551687","LOC724259","LOC724308","LOC724422","LOC100576838","LOC113218562",
#                 "LOC406114","LOC410530","LOC412810","Pgrp-s2")
#selected_VM <- c("LOC413946","LOC725878",'LOC413502',"LOC725404","LOC409805")
#selected_others <- c("LOC412976",'LOC100578215',"LOC102656552")


#####cluster在纵坐标
#EC$integrated_snn_res.0.1 <- factor(EC$integrated_snn_res.0.1,
#                                    levels = c("8", "13", "14",
#                                               "0", "2", "3", "4", "6", "7","11", "12","16",
#                                               "17", "18",
#                                               "10",
#                                               "1","5","9","15"))

#EC$integrated_snn_res.0.3 <- factor(EC$integrated_snn_res.0.3,
#                                    levels = c("10","14","18","19","20",
#                                               "5","22",
#                                               "0","2","3","6","7","9","11","17","21",
#                                               "8","12","28","29",
#                                               "1","4","13","15","16","23","24","25","26","27"))
EC$integrated_snn_res.0.3 <- factor(EC$integrated_snn_res.0.3,
                                    levels = c("7","9","11","3","6","2","0","17","21"))

EC$integrated_snn_res.0.3 <- factor(EC$integrated_snn_res.0.3,
                                    levels = c("5","22"))

##cluster在纵坐标
#p <- DotPlot(
#  EC,
#  features = c(
#    selected_ISC.EB, selected_EE, selected_EC, selected_VM, selected_others
#  ),
#  group.by = 'integrated_snn_res.0.1',
#  assay = 'RNA'
#)
#theme_bw()+
#  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5,angle = 45))+
#  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
# Save the plot to a PDF file
#pdf('/data2/liuhuiling/gut_snRNA/05cluster_marker_bar_umi/merge/res_0.1/Dotplot_0.1.pdf', width = 15, height = 6)
#print(p)
#dev.off()
##cluster在纵坐标--横坐标倾斜
#p <- DotPlot(EC,
#             features = c(selected_ISCEB, selected_EE, selected_EC, selected_VM, selected_others),
#             group.by = 'integrated_snn_res.0.3') +
#  theme_bw() +
#  theme(panel.grid = element_blank(),
#        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 45),
#        plot.margin = margin(t = 10, r = 10, b = 50, l = 10, unit = "pt")) + # 调整绘图区域边距
#  labs(x = NULL, y = NULL) +
#  guides(size = guide_legend(order = 3))+ 
#  coord_flip()
#p

p <- DotPlot(EC,
             features = EC_submarker,
             group.by = 'integrated_snn_res.0.3') +
  scale_color_gradientn(colors = c("blue", "white", "red"))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 45),
        plot.margin = margin(t = 10, r = 10, b = 50, l = 10, unit = "pt")) + # 调整绘图区域边距
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3))+ 
  coord_flip()
p
# 打印绘图
pdf('/data2/liuhuiling/gut_snRNA/05cluster_marker_bar_umi/merge/res_0.3/EE/EE_submarker.pdf',
    width = 5.5, height = 5)
print(p)
dev.off()




