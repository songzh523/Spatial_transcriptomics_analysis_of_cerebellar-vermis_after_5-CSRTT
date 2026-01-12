###R4.3.0####
#https://cloud.tencent.com/developer/article/1865540
#load library####
library(Seurat)
library(CellChat)
library(patchwork)
#set location####
datloc=datapath
resloc=resultpath
#load data####
kepsample="total" 
keptreat="Treat" 

cellchat <- createCellChat(object = subdata.input, meta = meta, group.by = "labels.Type")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels.Type")
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.mouse 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse) 

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 


netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
print(gg1 + gg2)
gg3<-netAnalysis_signalingRole_heatmap(cellchat,pattern = "incoming")+
  netAnalysis_signalingRole_heatmap(cellchat,pattern = "outgoing")
print(gg3)


save(cellchat,file = paste0(resloc,kepsample,"-",keptreat,"-cellchatnet.rds"))

##compare communication between "Treat" and "Ctrl"####
load(paste0(resloc,"total-Ctrl-cellchatnet.rds"))
cellchat.ctrl<-cellchat
load(paste0(resloc,"total-Treat-cellchatnet.rds"))
cellchat.treat<-cellchat

ctrl.net <- subsetCommunication(cellchat.ctrl)
treat.net <- subsetCommunication(cellchat.treat)
treat.net$source<-as.character(treat.net$source)
treat.net$target<-as.character(treat.net$target)
treat.net$ligand<-as.character(treat.net$ligand)
treat.net$receptor<-as.character(treat.net$receptor)
in.net<-treat.net[!(treat.net$source%in%ctrl.net$source & treat.net$target%in%ctrl.net$target & treat.net$ligand%in%ctrl.net$ligand & treat.net$receptor%in%ctrl.net$receptor),]
de.net<-ctrl.net[!(ctrl.net$source%in%treat.net$source & ctrl.net$target%in%treat.net$target & ctrl.net$ligand%in%treat.net$ligand & ctrl.net$receptor%in%treat.net$receptor),]
in.net$edge<-1
de.net$edge<-1

typs<-c("Cell-Cell Contact","ECM-Receptor","Secreted Signaling")
for(i in 1:length(typs)){
subin.net<-in.net[in.net$annotation==typs[i],]  
if(nrow(subin.net)>0){
subin.netmerge<-aggregate(subin.net$edge,by=list(subin.net$source,subin.net$target),FUN=sum)
colnames(subin.netmerge)<-c("source","target","edge")
}
subde.net<-de.net[de.net$annotation==typs[i],]
if(nrow(subde.net)>0){
subde.netmerge<-aggregate(subde.net$edge,by=list(subde.net$source,subde.net$target),FUN=sum)
colnames(subde.netmerge)<-c("source","target","edge")
}
}

in.netmerge<-aggregate(in.net$edge,by=list(in.net$source,in.net$target),FUN=sum)
colnames(in.netmerge)<-c("source","target","edge")
de.netmerge<-aggregate(de.net$edge,by=list(de.net$source,de.net$target),FUN=sum)
colnames(de.netmerge)<-c("source","target","edge")
####heatmap####
allcelltypes<-c("Granular_cell","Purkinje_cell","Golgi_cell","Basket_cell","Stellate_cell","Candelabrum_cell","Unipolar_brush_cell","Microglia","Astrocytes","Bergmann_glia_cell","Oligodendrocyte","Oligodendrocyte_precursor_cell","Choroid","Choroid_plexus_cell","Endo_stalk","Endothelial_cell","MLI")

inps1<-reshape2::dcast(in.netmerge,source~target,value.var="edge")
inps2<-inps1[,-1]
rownames(inps2)<-inps1[,1]
inps2[is.na(inps2)]=0
View(inps2)
pheatmap::pheatmap(inps2,scale = "none",treeheight_col = 0,treeheight_row = 0) #,color = c("white",colorRampPalette(colors = c("#e5c5bb","firebrick3"))(length(bk)-1)),legend_breaks=seq(1,7,2),breaks=bk)

result <- as.data.frame(matrix(0, nrow = length(allcelltypes), ncol = length(allcelltypes)))
rownames(result) <- allcelltypes
colnames(result) <- allcelltypes
existing_rows <- allcelltypes %in% rownames(inps2)
existing_cols <- allcelltypes %in% colnames(inps2)
result[existing_rows, existing_cols] <- inps2[match(allcelltypes[existing_rows], rownames(inps2)),match(allcelltypes[existing_cols], colnames(inps2))]
library(RColorBrewer)
bk <- c(0,seq(0.9,max(result),by=0.001))

pheatmap::pheatmap(result,scale = "none",cluster_rows = T,cluster_cols = T,treeheight_col = 0,treeheight_row = 0,color = c("white",colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(bk)-1)),legend_breaks=seq(1,max(inps2),2),breaks=bk) #,color = c("white",colorRampPalette(colors = c("#e5c5bb","firebrick3"))(length(bk)-1)),legend_breaks=seq(1,7,2),breaks=bk)



bk <- c(0,seq(0.9,max(inps2),by=0.001))
library(RColorBrewer)
xord<-c("Granular_cell","Choroid","Astrocytes","Choroid_plexus_cell","Basket_cell","Oligodendrocyte","MLI","Stellate_cell","Purkinje_cell","Endothelial_cell")
yord<-c("Purkinje_cell","Stellate_cell","Choroid","Oligodendrocyte","Endothelial_cell","Basket_cell","Astrocytes","MLI","Granular_cell","Choroid_plexus_cell")
inps3<-inps2[yord,xord]
pheatmap::pheatmap(inps2,scale = "none",cluster_rows = T,cluster_cols = T,treeheight_col = 0,treeheight_row = 0,color = c("white",colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(bk)-1)),legend_breaks=seq(1,max(inps2),2),breaks=bk)#+ggplot2::labs(x= "Ligands",y = "Receivers")



