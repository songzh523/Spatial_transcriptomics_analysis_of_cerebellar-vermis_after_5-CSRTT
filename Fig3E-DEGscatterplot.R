suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library(clusterProfiler))
library(VennDiagram)
library(grid)
library(stringr)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)

#load data####
celltype="Granular_cell" 
celltypefix="GC" 

###scatter of gene between ctrl and 5csrt[lobus; posterior]####
firstgene#DEGs in anterior
secondgene#DEGs in posterior
firsts<-firstgene[firstgene$cluster=="5CSRT",c("gene","avg_log2FC")]
firsts$type<-"lobus"
seconds<-secondgene[firstgene$cluster=="5CSRT",c("gene","avg_log2FC")]
seconds$type<-"posterior"
genesexp<-merge(firsts,seconds,by="gene")
p0<-ggpubr::ggscatter(genesexp,"avg_log2FC.y","avg_log2FC.x",xlim=c(-1,1),ylim=c(-1,1),size=0.6,xlab = "posterior Log2FC (5CSRT vs Ctrl)",ylab = "anterior Log2FC (5CSRT vs Ctrl)")


selected<-c("Itpr1","Unc80","Dgkz","Rgs8","Gad1","Dctn1","Atp1a1","Bin1","Clstn3","Tmem108","Grin1","Vamp1","Slc6a17","Plp1","Mobp","Coro1a","Calb2","Aopep","Ywhae","Jph3")

genesexp$labels<-ifelse(genesexp$gene %in% selected,genesexp$gene,"")
genesexp$show<-ifelse((genesexp$gene %in% firstsigs  & !genesexp$gene %in% secondsigs),"lobus",ifelse((genesexp$gene %in% secondsigs & !genesexp$gene %in% firstsigs),"posterior",ifelse((genesexp$gene %in% secondsigs & genesexp$gene %in% firstsigs),"overlap","zzz")))
genesexp$bond<-ifelse(genesexp$gene %in% selected,"show","none")

p<-ggpubr::ggscatter(genesexp,"avg_log2FC.y","avg_log2FC.x",fill = "bond",color="show",xlim=c(-1,1),ylim=c(-1,1),xlab = "PLGCs Log2FC (5-CSRTT vs Naive)",ylab = "ALGCs Log2FC (5-CSRTT vs Naive)")
p+geom_text_repel(data = genesexp, aes(label = labels),size = 3,max.overlaps=1000)+scale_fill_manual(values = c("black","gray"))+scale_color_manual(values = c("#CA6E82","black","#638ACD","gray"))+theme(legend.position='none')+geom_abline(intercept = 0, slope = 1,linetype=2)

