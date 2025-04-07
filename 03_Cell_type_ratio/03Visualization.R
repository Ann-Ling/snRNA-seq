# Package and environment
rm(list = ls())
library(tidyr)
library(dplyr)
library(ggplot2)

#Input
data0 <- read.table("/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/03pb_da/03_1_cell_cnt_freq.txt",
                 header = T)
data0$group_id <- as.factor(data0$group_id)
data0$cluster_id <- as.factor(data0$cluster_id)

# Output
library(ggplot2)
data0$group_id <- factor(data0$group_id, 
                         levels = c("CV", "GF"), ordered = TRUE)
data$cluster_id <- as.factor(data0$cluster_id)
p <- ggplot(data0, aes(x = group_id, y = frequency, fill = group_id)) +
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "#000000")) +  
  geom_boxplot(size = 1, outlier.fill = "white", outlier.color = "white") +  
  scale_fill_manual(values = c("#0B60B0", "#95D2B3")) +  
  geom_jitter(aes(fill = group_id), width = 0.2, shape = 21, size = 3) +  
  scale_color_manual(values = c("#000000", "#000000", "#000000")) +  
  ggtitle("region and cluster") + 
  theme_bw() +  
  theme(
    legend.position = "none",  
    axis.text.x = element_text(colour = "#000000", family = "Arial", size = 14),  
    axis.text.y = element_text(family = "Arial", size = 14, face = "plain"), 
    axis.title.y = element_text(family = "Arial", size = 14, face = "plain"),  
    axis.title.x = element_text(family = "Arial", size = 14, face = "plain"), 
    plot.title = element_text(family = "Arial", size = 15, face = "bold", hjust = 0.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  ) +
  ylab("Frequency") +  
  facet_wrap(~ cluster_id)
p
pdf('./cell_ratio.pdf')
print(p)
dev.off()
