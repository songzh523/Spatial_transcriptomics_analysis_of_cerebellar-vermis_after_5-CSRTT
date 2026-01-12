#R(4.5.1)-load library####
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library(clusterProfiler))
library(Hmisc)
library(liger)
library(ggplot2)
library(ggrepel)
#load data location####
datloc=youpath
gmtloc=gmtpath
type="5CSRTvsContrl"

#load data
files0<-list.files(path = datloc,pattern = "*-sigmarkers.txt") 
filesnam<-sapply(strsplit(files,"-",fixed=T), "[[",1)
#Fig3F:significant genes among each leaf####
for(i in 1:length(files)){
  if(i==1){
    dat<-read.table(paste(datloc,files[i],sep=""))
    dat$region<-filesnam[i]
    dat$updown<-ifelse(dat$p_val<0.05 & dat$avg_log2FC>0.25,"UP",ifelse(dat$p_val<0.05 & dat$avg_log2FC<(-0.25),"DOWN","other"))
  }else{
    dat0<-read.table(paste(datloc,files[i],sep=""))
    dat0$region<-filesnam[i]
    dat0$updown<-ifelse(dat0$p_val<0.05 & dat0$avg_log2FC>0.25,"UP",ifelse(dat0$p_val<0.05 & dat0$avg_log2FC<(-0.25),"DOWN","other"))
    dat<-rbind(dat,dat0)
  }
}

inp0<-as.data.frame(table(dat$region,dat$updown))
inp0<-inp0[inp0$Freq!=0,]
inp1<-inp0[inp0$Var2!="other" & inp0$Var2!="DOWN",]


firstCap<-function(x){
  capitalize(tolower(x))
}

gmtfiles<-list.files(path=gmtloc,pattern = "*.txt")
gmtfilesnam<-sub("HP_","",sub("\\..*","",gmtfiles))
gs.lists<-list()
for(f in 1:length(gmtfiles)){
gene_sets <- read.table(paste0(gmtloc,gmtfiles[f]))
gene_names <- firstCap(unlist(gene_sets))
gs.lists[[f]]=gene_names
names(gs.lists)[f]=gmtfilesnam[f]
}
files0<-list.files(path = datloc,pattern = "*-sigmarkers.txt") 
filesnam<-sapply(strsplit(files,"-",fixed=T), "[[",1)

for(i in 1:length(files)){
    dat<-read.table(paste(datloc,files[i],sep=""))
    dat$region<-filesnam[i]
    dat$updown<-ifelse(dat$p_val<0.05 & dat$avg_log2FC>0.25,"UP",ifelse(dat$p_val<0.05 & dat$avg_log2FC<(-0.25),"DOWN","other"))
    geneList0<-dat$avg_log2FC[dat$updown=="UP"]
    names(geneList0)<-dat$gene
    geneList<-sort(geneList0,decreasing = T)
    gseavals<-iterative.bulk.gsea(values=geneList,set.list = gs.lists) #total or #select one [c(6,17,13,10,4,12)]
    if(i==1){
      dat<-read.table(paste0(figloc,"gsearesult/",files[i]))
      dat$region<-filesnam[i]
    }else{
      tmp<-read.table(paste0(figloc,"gsearesult/",files[i]))
      tmp$region<-filesnam[i]
      dat<-rbind(dat,tmp)  
    }
}

dat$Path<-gsub("[1-9]$","",rownames(dat))
library(pheatmap)
library(reshape2)
inp<-dat[,c("region","Path")]
inp$logp<-(-log10(dat$p.val))
inp1<-melt(inp)
inpheat0<-reshape2::dcast(inp1,region~Path,value.var="value")
inpheattotal<-inpheat0
rownames(inpheattotal)<-inpheat0

pheatmap(inpheattotal, display_numbers = matrix(ifelse(inpheattotal > 2, "**", ifelse(inpheattotal > 1.3,"*","")),nrow(inpheattotal)),cluster_row=F,treeheight_col =0,color=colorRampPalette(c("#3a97c0", "#FFFFFF", "#ec9594"))(50),border="white")

###Fig3F: DEGs number plot####
ggplot(data = inp1, mapping = aes(reorder(Var1,Freq), y = Freq,fill=Var2)) + geom_bar(stat = 'identity', position = 'identity') +theme_classic()+ scale_fill_manual(values = c('blue','red'), guide = FALSE) + ylab('Number')+xlab('')+geom_text(data = inp1, aes(x=Var1,  y=Freq, label=Freq))+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.line.x = element_blank(),axis.ticks.y  = element_blank(),axis.ticks.x  = element_blank())+scale_fill_manual(values=c("#b1d4e5","#f47b7b"))+labs(fill="",title = type)



