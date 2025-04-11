#Package
library(ggplot2)
library(tidyverse)
library(ggrepel)

#Environment
rm(list=ls())

#Input
setwd("/data2/liuhuiling/gut_snRNA/07cell_type_group/AM_merge/ileum/DEG")
df <- read.table("markers.BasedOncondition.txt",header = T)
head(df)

logFC_t <- 1
df$label <- ifelse(df$p_val>0.05,'stable',
                   ifelse( df$avg_log2FC > logFC_t,'CV',
                           ifelse( df$avg_log2FC< -logFC_t,'GF','stable')))
head(df)

df$cluster <- factor(df$cluster,levels=c("ISC/EB","EE","EC",'VM',"N"))
dfbar<-data.frame(x=c("ISC/EB","EE","EC",'VM',"N"),
                  y=c(1,1,1,1,1))
dfbar1<-data.frame(x=c("ISC/EB","EE","EC",'VM',"N"),
                   y=c(-4.8,-4.8,-4.8,-4.85,-4.8))
dfbar$x <- factor(dfbar$x,levels=c("ISC/EB","EE","EC",'VM',"N"))
dfbar1$x <- factor(dfbar1$x,levels=c("ISC/EB","EE","EC",'VM',"N"))

#Plot
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = df,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)
p2

dfcol<-data.frame(x=c("ISC/EB","EE","EC",'VM',"N"),
                  y=c(0,0,0,0,0),
                  label=c("ISC/EB","EE","EC",'VM',"N"))
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.35,
                     color = "black",
                     alpha = 1,
                     show.legend = F)
p3

p4 <- p3 +
  scale_color_manual(name=NULL,
                     values = c("#B33030","#3271ae","#6B728E"))
p4

p5 <- p4+
  labs(x="cell type",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="black")
p5

p6 <- p5+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               linewidth = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p6

#Output
pdf('./DEG_volcano.pdf',width = 10,height = 8)
print(p6)
dev.off()

