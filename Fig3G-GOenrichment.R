library(ggpubr)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)

resgo1t<-#enriched GO terms in GCs of anterior in 5CSRTT treatment
resgo1c<-#enriched GO terms in GCs of anterior in ctrl treatment 
resgo2t<-#enriched GO terms in GCs of posterior in 5CSRTT treatment
resgo2c<-#enriched GO terms in GCs of posterior in ctrl treatment
resgo1<-rbind(resgo1t[1:10,],resgo1c[1:10,])
resgo1$logp<-(-log10(as.numeric(resgo1$p.adjust)))
resgo1$type<-c(rep("5CSRT",10),rep("Ctrl",10))
resgo1$Count<-as.numeric(resgo1$Count)
resgo1$type<-factor(resgo1$type,levels = c("Ctrl","5CSRT"))

ggpubr::ggscatter(resgo1,"Description","type",color = "logp",size = "Count",sort.by.groups = T,sort.val="desc",orientation="horizontal",ylab="",xlab = "",title = paste("lobus-",celltypefix,sep=""),legend = "bottom",legend.title="-log10(FDR)")+ scale_color_gsea()#+ scale_color_gsea()

resgo2<-rbind(resgo2t[1:10,],resgo2c[1:10,])
resgo2$logp<-(-log10(as.numeric(resgo2$p.adjust)))
resgo2$type<-c(rep("5CSRT",10),rep("Ctrl",10))
resgo2$Count<-as.numeric(resgo2$Count)
resgo2$type<-factor(resgo2$type,levels = c("Ctrl","5CSRT"))
resgo3<-resgo2[!duplicated(resgo2$Description),]

ggpubr::ggscatter(resgo3,"Description","type",color = "logp",size = "Count",sort.by.groups = T,sort.val="desc",orientation="horizontal",ylab="",xlab = "",title = paste("posterior-",celltypefix,sep=""),legend = "bottom",legend.title="-log10(FDR)")+ scale_color_gsea()#+ scale_color_gsea()

