#library R4.1.3 Seurat; R4.5.1 RRHO ####
suppressMessages(library(RColorBrewer))
library(Seurat) #version 4.3.0
library(dplyr)
library(ggplot2)
library(RRHO)
library(lattice)
library(ggplot2)
library(lattice)
library("circlize")
library("RColorBrewer")
library(scales)
library(ggpubr)
library(ggsci)
#data location####
resloc=yourpath
##fig1.e rank-rank hypergeometric overlap####
type=c("5CSRTvsCtrl")

  files<-list.files(path = resloc,pattern = "*-sigmarkers-all.txt")
  filesnam<-sapply(strsplit(files,"-",fixed=T), "[[",1)

  for(i in 1:(length(files)-1)){
    for(j in (i+1):length(files)){
      
      dat1t<-read.table(paste(resloc1,files[i],sep=""))
      dat2t<-read.table(paste(resloc1,files[j],sep=""))
      dat1<-dat1t[dat1t$avg_log2FC>0 & dat1t$p_val<0.05,] 
      dat2<-dat2t[dat2t$avg_log2FC>0 & dat2t$p_val<0.05,]
      dat1a<-dat1[dat1$gene %in% dat2$gene,]
      dat2a<-dat2[dat2$gene %in% dat1$gene,]
      
      list1<- data.frame(GeneIdentifier=dat1a$gene,
                         RankingVal=dat1a$avg_log2FC)
      list2<- data.frame(GeneIdentifier=dat2a$gene,
                         RankingVal=dat2a$avg_log2FC)
      
      RRHO.example <-  RRHO(list1,list2,BY=TRUE, alternative='two.sided')
      inp<-RRHO.example$hypermat
      inp[is.infinite(inp)]=range(inp,finite = T)[2]
      range(inp)
      print(paste(filesnam[i],filesnam[j],range(inp)[2]))
      p<-levelplot(inp,scales=list(draw=F),col.regions = colorRampPalette(c("blue","white","red")),at = c(-Inf,seq(0, 800, by = 200),Inf),xlab=filesnam[i],ylab=filesnam[j])
      print(p) 

    }
  }
  dev.off()
 
  