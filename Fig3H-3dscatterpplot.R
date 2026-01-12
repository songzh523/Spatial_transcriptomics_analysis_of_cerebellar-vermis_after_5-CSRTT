#R4.1.3 
library(Seurat)
library(scatterplot3d)
library(ggrepel)
library(rgl)
options(rgl.printRglwidget = TRUE)

###CV and Tau gene expression specific####
datloc=yourpath
celltype="Granular_cell" 
celltypefix="GC" 
cols<-c("#C57084","#648BD0")

mergedcom #seurat object

selected<-c("Itpr1","Unc80","Dgkz","Rgs8","Gad1","Dctn1","Atp1a1","Bin1","Clstn3","Tmem108","Grin1","Vamp1","Slc6a17","Plp1","Mobp","Coro1a","Calb2","Aopep","Ywhae","Jph3")
leafs<-unique(mergedcom@meta.data$leaf)
samples<-unique(as.character(mergedcom@meta.data$sample))
for(l in 1:length(leafs)){
for(s in 1:length(samples)){
a<-mergedcom@assays$Spatial@counts[,mergedcom@meta.data$leaf==leafs[l] & mergedcom@meta.data$sample==samples[s]]
tau<-apply(a,1, function(x) sum(1-(x/max(x)))/(ncol(a)-1))
cvs<-apply(a, 1, function(x) sd(x)/mean(x))
if(l==1 & s==1){
  outs<-as.data.frame(cvs)
  colnames(outs)<-paste(leafs[l],samples[s],sep="-")
  outs1<-as.data.frame(tau)
  colnames(outs1)<-paste(leafs[l],samples[s],sep="-")
}else{
  outs<-cbind(outs,as.data.frame(cvs))
  colnames(outs)[ncol(outs)]<-paste(leafs[l],samples[s],sep="-")
  outs1<-cbind(outs1,as.data.frame(tau))
  colnames(outs1)[ncol(outs1)]<-paste(leafs[l],samples[s],sep="-")
}
}
}


b1<-cbind(outs[firstuniq,1:2],outs1[firstuniq,1:2]) ##lobus Tau+CV
colnames(b1)<-paste(sub(".*-","",colnames(b1)),rep(c("CV","Tau"),each=2),sep="_")
b2<-cbind(outs[seconduniq,3:4],outs1[seconduniq,3:4]) ##posterior Tau+CV
colnames(b2)<-paste(sub(".*-","",colnames(b2)),rep(c("CV","Tau"),each=2),sep="_")

b1$lab<-rownames(b1)
b1$testcv<-(b1$`5CSRT_CV`-b1$`Ctrl_CV`)/b1$`Ctrl_CV`
b1$testtau<-(b1$`5CSRT_Tau`-b1$`Ctrl_Tau`)/b1$`Ctrl_Tau`
b1$gonum<-gonumbercount(sigresgo1,firstuniq)
b1$gonum[is.na(b1$gonum)]=0
b1$s<-"AL"


b2$lab<-rownames(b2)
b2$testcv<-(b2$`5CSRT_CV`-b2$`Ctrl_CV`)/b2$`Ctrl_CV`
b2$testtau<-(b2$`5CSRT_Tau`-b2$`Ctrl_Tau`)/b2$`Ctrl_Tau`
b2$gonum<-gonumbercount(sigresgo2,seconduniq)
b2$gonum[is.na(b2$gonum)]=0
b2$s<-"PL"


b1$type<-"AL"
b2$type<-"PL"
bb<-rbind(b1,b2)
bb$s<-factor(bb$s,levels = c("AL","PL"))


cols<-cols[as.numeric(as.factor(bb$s))]

s3d<-scatterplot3d(bb[,c(8,7,6)],pch=16,angle = 60,color = cols,box = F,xlab="GO term number",ylab="Tau variation",zlab="CV variation",main = "ALandPL")
text(s3d$xyz.convert(bb[,c(8,7,6)]),labels = bb$lab,offset = 0.6,cex = 0.6)
legend(s3d$xyz.convert(8,7,6), legend = levels(bb$s),col =  c("#C57084","#648BD0"), pch = 16)

s3d<-scatterplot3d(bb[,c(8,7,6)],pch=16,angle = 60,color = cols,xlab="GO term number",ylab="Tau variation",zlab="CV variation",main = "ALandPL")
text(s3d$xyz.convert(bb[,c(8,7,6)]),labels = bb$lab,offset = 0.6,cex = 0.6)
legend(s3d$xyz.convert(8,7,6), legend = levels(bb$s),col =  c("#C57084","#648BD0"), pch = 16)

