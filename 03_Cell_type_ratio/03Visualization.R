# Package and environment
rm(list = ls())
library(tidyr)
library(dplyr)
library(ggplot2)
df <- read.table("/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/03pb_da/03_1_cell_cnt_freq.txt",
                 header = T)
# 计算每个 cluster 在 CV 和 GF 组的平均值
averaged_data <- df %>%
  group_by(cluster_id, group_id) %>%
  summarize(mean_value = mean(frequency), .groups = 'drop') %>%
  pivot_wider(names_from = group_id, values_from = mean_value, names_prefix = "mean_")

# 绘制散点图
# 假设 cluster_colors 是一个命名向量，包含 cluster_id 和对应颜色
cluster_colors <- c("0" = "#f9cb8b","1"="#eef7f2", "2" = "#f28e16", "3" = "#f58825","4" = "#d8e3e7","5" = "#61649f",
                    "6" = "#f27635","7" = "#fc6315","8" = "#1677b3",
                    "9" = "#d85916", "10" = "#2c9678","11" = "#e46828", "12" = "#1661ab","13" = "#ede3e7","14" = "#248067",
                    "15" = "#d1c2d3","17" = "#b7511d","18" = "#83a78d","19" = "#6e8b74","20" = "#314a43","21" = "#a6522c",
                    "23" = "#dddddd","24" = "#ccccd6","25" = "#c0c4c3","26" = "#74787a","27" = "#5e616d")  # 根据实际情况修改

# 绘制散点图，使用指定的颜色
# 计算每个点与45°线的距离
averaged_data <- averaged_data %>%
  mutate(distance_to_line = abs(mean_CV - mean_GF))
# 定义一个阈值，决定何为"靠近" 45° 线
threshold <- 0.02  # 根据需要调整阈值
# 重新建立一个point_color列，方便赋色
averaged_data <- averaged_data %>%
  mutate(cluster_id = factor(cluster_id),  # 确保 cluster_id 为因子
         point_color = ifelse(distance_to_line < threshold, "close", 
                              as.character(cluster_id)))  # 使用 cluster_id 直接赋值
# 画图
p <- ggplot(averaged_data, aes(x = mean_CV, y = mean_GF)) +
  geom_point(aes(color = point_color), size = 4) +  # 使用新的颜色列
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +  # 45° 分界线
  geom_text(aes(label = ifelse(distance_to_line >= threshold, point_color, "")), 
            vjust = -1, hjust = 0.5, size = 4) +  # 仅在不靠近线的点上添加标签
  scale_color_manual(values = c("close" = "grey", cluster_colors), name = "Cluster") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(size = 0.3),
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.line.y = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black")
  ) +
  coord_cartesian(clip = "off")

p
pdf('/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/ileum/ileum_cluster/03pb_da/cluster_ratio.pdf',width=7,height=5)
print(p)
dev.off()
write.csv(averaged_data,
          "/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/ileum/ileum_cluster/03pb_da/cluster_ratio data.csv",
          row.names = FALSE)



