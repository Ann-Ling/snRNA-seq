#Calculate the proportion of different branches
#Package
library(dplyr)

# Extracting branch and group information of cells
cell_data <- pData(mycds)

# Count cells by branch and group
branch_group_counts <- cell_data %>%
  group_by(State, orig.ident) %>%
  summarise(count = n(), .groups = "drop") 

# Calculate the proportion of groups in each branch
branch_group_ratio <- branch_group_counts %>%
  group_by(State) %>%
  mutate(ratio = count / sum(count))
print(branch_group_ratio)
write.table(branch_group_ratio,"/data2/liuhuiling/gut_snRNA/09pseudotime/ileum/State ratio.txt",sep="\t",
            quote = F,row.names = F)

#Visualization
p <- ggplot(branch_group_ratio, aes(x = "", y = ratio, fill = orig.ident)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~sub_cell_type) +
  labs(fill = "Group", y = "Proportion", x = NULL) +
  theme_void() +
  theme(legend.position = "right") +
  ggtitle("Proportion of CV and GF in Each Branch (State)")+
  scale_fill_manual(values = c("CV" = "#0B60B0", "GF" = "#95D2B3"))
p
pdf("/data2/liuhuiling/gut_snRNA/09pseudotime/ileum/subtype ratio.pdf",width = 20,height=10)
print(p)
dev.off()

##########################
branch_group_counts <- cell_data %>%
  group_by(cell_type, orig.ident, State) %>%
  summarise(count = n(), .groups = "drop")

# Calculate the proportion of each cell type in the group
branch_group_ratio <- branch_group_counts %>%
  group_by(State) %>%
  mutate(ratio = count / sum(count))
print(branch_group_ratio)

#Save
write.table(branch_group_ratio,"/data2/liuhuiling/gut_snRNA/09pseudotime/ileum/cell type and state ratio(by state).txt",sep="\t",
            quote = F,row.names = F)

