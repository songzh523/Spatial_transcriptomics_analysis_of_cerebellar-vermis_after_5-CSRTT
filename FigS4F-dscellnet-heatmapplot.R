datloc=youpath

inps #compared network
inps1<-reshape2::dcast(inps,from~to,value.var="edge")
inps2<-inps1[,-1]
rownames(inps2)<-inps1[,1]

inps2[is.na(inps2)]=0
View(inps2)
pheatmap::pheatmap(inps2,scale = "none",treeheight_col = 0,treeheight_row = 0) 
bk <- c(0,seq(0.9,7,by=0.001))
library(RColorBrewer)
pheatmap::pheatmap(inps2,scale = "none",treeheight_col = 0,treeheight_row = 0,color = c("white",colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(bk)-1)),legend_breaks=seq(1,7,2),breaks=bk)#+ggplot2::labs(x= "Ligands",y = "Receivers")

