require(ggplot2)
require(RColorBrewer)
require(reshape2)
require(cowplot)
require(stringr)
setwd("/data2/liuhuiling/gut_snRNA/06cross_species/merge/res_0.3")
source("2017-08-28-runMN-US.R")
source("2017-08-28-runMN-US.pearson.R")

data<-read.table("merged.exp.txt")
rownames(data) <- str_split_fixed(rownames(data),'\\(',2)[,1]
pheno<-read.table("merged.pheno.cell.txt",header=T)
celltypes <- unique(pheno$Celltype)
var.genes<-read.delim("marker.list",header=F)$V1
var.genes=get_variable_genes(data, pheno)
AUROC.matrix.s=run_MetaNeighbor_s(var.genes, data, celltypes, pheno)
AUROC.matrix.p=run_MetaNeighbor_p(var.genes, data, celltypes, pheno)
AUROC.data.s<-melt(AUROC.matrix.s,value.name = "AUROCs")
AUROC.data.p<-melt(AUROC.matrix.p,value.name = "AUROCp")
res<-merge(AUROC.data.p,AUROC.data.s,by=c("Var1","Var2"))
write.table(res,"AUROC.z.txt",sep="\t",quote=F,row.names=F)


data<-read.table("AUROC.z.txt",header=T)

plot_theme<-theme(panel.background = element_blank(),
axis.line = element_line(size=0.1),axis.ticks = element_line(size=0.1),
axis.text = element_text(size=2.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
axis.ticks.length = unit(1, "pt"),plot.title = element_text(size = 3, face = "bold",margin=margin(0,0,4,0)),
plot.margin = unit(c(0,0,0,0),"cm"),legend.key.size = unit(0.2, 'cm'),legend.key.height = unit(0.2, 'cm'),
legend.key.width = unit(0.2, 'cm'),legend.title = element_text(size=4),legend.text = element_text(size=4))


MphaDmel<-subset(data, grepl("Mpha",Var1) & grepl("Amel",Var2))

MphaDmel$AUROCp2<-ifelse(MphaDmel$AUROCp<=0.5, "", round(MphaDmel$AUROCp,2))

MphaDmel$Var1<-factor(MphaDmel$Var1,levels=c("Dros_cardia","Dros_ISC.EB","Dros_AstA.EE","Dros_NPF.EE","Dros_AstC.EE","Dros_dEC","Dros_aEC1",
"Dros_aEC2","Dros_aEC3","Dros_aEC4","Dros_mEC","Dros_copper.and.iron.cell","Dros_LFC","Dros_pEC1","Dros_pEC2","Dros_pEC3",
"Dros_EC.like1","Dros_EC.like2","Dros_EC.like3","Dros_others","Dros_unk1","Dros_unk2"))
MphaDmel$Var2<-factor(MphaDmel$Var2,levels=c("AM_c10","AM_c14","AM_c18","AM_c19","AM_c20",
                                             "AM_c5","AM_c22",
                                             "AM_c0","AM_c2","AM_c3","AM_c6","AM_c7","AM_c9","AM_c11",
                                             "AM_c17","AM_c21",
                                             "AM_c8","AM_c12","AM_c28","AM_c29",
                                             "AM_c1","AM_c4","AM_c13","AM_c15","AM_c16",
                                             "AM_c23","AM_c24","AM_c25","AM_c26","AM_c27"))


p1<-ggplot(MphaDmel, aes(Var1,Var2, fill = AUROCp))+ geom_tile(color = "white")+scale_fill_distiller(palette = "RdBu", limits = c(0,1))+geom_text(aes(label=AUROCp2),size=0.9)+coord_fixed()+labs(x="",y="")+plot_theme
p2<-ggplot(MphaDmel, aes(Var1,Var2, fill = AUROCs))+ geom_tile(color = "white")+scale_fill_distiller(palette = "RdBu", limits = c(0,1))+geom_text(aes(label=round(AUROCs,3)),size=1)+coord_fixed()+labs(x="",y="")+plot_theme
ggsave(filename="Dros_Amel.z.auroc.pdf",plot=plot_grid(p1,p2,nrow = 1),height=6,width=12)
ggsave(filename="AUROCs.pdf",p2,height=6,width=6)
