##Figure 2

# load the following objects from https://osf.io/wxpgn/
##Figure 2A
Tert  <- readRDS("Tert_Subtypes.rds")
# UMAP-plot of seurat object based on clusters
DimPlot(Tert, reduction="umap", group.by="Subtype",label = TRUE, pt.size = .1)+ NoLegend() 

##Figure 2B
# UMAP-plot of seurat object based on cell lines
DimPlot(Tert, reduction="umap", group.by="Dataset",label = TRUE, pt.size = .1)+ NoLegend() 

##Figure 2C
#Heatmap for the top 10 enriched genes per each cluster
Idents(object = Tert) <- "RNA_snn_res.0.15"
Tert.markers <- FindAllMarkers(Tert,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Tert.markers <- Tert.markers[Tert.markers$p_val_adj < 0.05,]

top10 <- Tert.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Tert, features = top10$gene)

#Number of cells in each cluster
table(Tert@meta.data$RNA_snn_res.0.15)

##Figure 2D
#Overlap of genesfor clusters and cell lines (Jaccard index)
#Load the following objects 
 Tert_markers <- read.delim("Tert.markers.txt", h=T) 
 Tert_markers_Dataset <- read.delim("Tert.markers_Dataset.txt", h=T) 
 
 Tert<- merge(Tert_markers, Tert_markers_Dataset, all.x=TRUE, all.y=TRUE)

 Cluster0 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==0 ,c('gene')])
 Cluster1 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==1 ,c('gene')])
 Cluster2 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==2 ,c('gene')])
 Cluster3 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==3 ,c('gene')])
 Cluster4 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==4 ,c('gene')])
 Cluster5 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==5 ,c('gene')])
 Cluster6 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==6 ,c('gene')])
 Cluster7 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==7 ,c('gene')])
 Cluster8 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==8 ,c('gene')])
 Cluster9 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==9 ,c('gene')])
 Cluster10 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster==10 ,c('gene')])
 
 AD10 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster=="AD10" ,c('gene')])
 DD8 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster=="DD8" ,c('gene')])
 CB4 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster=="CB" ,c('gene')])
 CD8 <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster=="CD" ,c('gene')])
 Msc <- data.frame(Tert[Tert$p_val_adj<0.05& Tert$cluster=="Msc" ,c('gene')])
 
  names(Cluster0)<-"Symbol"
  names(Cluster1)<-"Symbol"
  names(Cluster2)<-"Symbol"
  names(Cluster3)<-"Symbol"
  names(Cluster4)<-"Symbol"
  names(Cluster5)<-"Symbol"
  names(Cluster6)<-"Symbol"
  names(Cluster7)<-"Symbol"
  names(Cluster8)<-"Symbol"
  names(Cluster9)<-"Symbol"
  names(Cluster10)<-"Symbol"
  names(AD10)<-"Symbol"
  names(DD8)<-"Symbol"
  names(CB4)<-"Symbol"
  names(CD8)<-"Symbol"
  names(Msc)<-"Symbol"


Tert_markers <- Tert_markers[Tert_markers$p_val_adj < 0.05,]
df<- list()
df[[1]] <- Cluster0
df[[2]] <- Cluster1
df[[3]] <- Cluster2
df[[4]] <- Cluster3
df[[5]] <- Cluster4
df[[6]] <- Cluster5
df[[7]] <- Cluster6
df[[8]] <- Cluster7
df[[9]] <- Cluster8
df[[10]] <- Cluster9
df[[11]] <- Cluster10
df[[12]] <- AD10
df[[13]] <- DD8
df[[14]] <- CB4
df[[15]] <- CD8
df[[16]] <- Msc
names(df) <- c("Cluster0","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8","Cluster9","Cluster10","AD10","DD8","CB4","CD8","Msc") 


sq <- seq(max(length(Cluster0$Symbol),length(Cluster1$Symbol),length(Cluster2$Symbol), length(Cluster3$Symbol), length(Cluster4$Symbol),length(Cluster5$Symbol),length(Cluster6$Symbol), length(Cluster7$Symbol),length(Cluster8$Symbol),length(Cluster9$Symbol),length(Cluster10$Symbol), length(AD10$Symbol),length(DD8$Symbol),length(CB4$Symbol),length(CD8$Symbol), length(Msc$Symbol)))
df <- data.frame(Cluster0$Symbol[sq], Cluster1$Symbol[sq], Cluster2$Symbol[sq],Cluster3$Symbol[sq],Cluster4$Symbol[sq], Cluster5$Symbol[sq], Cluster6$Symbol[sq],Cluster7$Symbol[sq],Cluster8$Symbol[sq], Cluster9$Symbol[sq], Cluster10$Symbol[sq],AD10$Symbol[sq],DD8$Symbol[sq], CB4$Symbol[sq], CD8$Symbol[sq],Msc$Symbol[sq])


input.variables =df

m = matrix(data = 0, nrow = length(input.variables), ncol = length(input.variables))
for (r in 1:length(input.variables)) {
        for (c in 1:length(input.variables)) {
                if (c == r) {
                        m[r,c] = 1
                } else if (c > r) {
                        m[r,c] =length(intersect(input.variables[,r], input.variables[,c]))/length(union(input.variables[,r],input.variables[,c]))
                        
                }
        }
}

variable.names = sapply(input.variables, attr, "label")
colnames(m) = colnames(df)
rownames(m) = colnames(df)   

jaccards = m


library(RColorBrewer)
library(dplyr)

col <- designer.colors(n=50, col=c("white","red", "black"))
col_breaks <- seq(0,1,length=51) 
heatmap.2(jaccards, scale = "none", col = col, breaks=col_breaks ,cexRow=1,cexCol = 0.9, 
          Colv = F, Rowv = F, margins = c(10,10),trace="none", dendrogram='none')
		  
#Number of marker genes in each cluster and cell line
length(unique(df$Cluster0.Symbol.sq.))
length(unique(df$Cluster1.Symbol.sq.))
length(unique(df$Cluster2.Symbol.sq.))
length(unique(df$Cluster3.Symbol.sq.))
length(unique(df$Cluster4.Symbol.sq.))
length(unique(df$Cluster5.Symbol.sq.))
length(unique(df$Cluster6.Symbol.sq.))
length(unique(df$Cluster7.Symbol.sq.))
length(unique(df$Cluster8.Symbol.sq.))
length(unique(df$Cluster9.Symbol.sq.))
length(unique(df$Cluster10.Symbol.sq.))
length(unique(df$AD10.Symbol.sq.))
length(unique(df$DD8.Symbol.sq.))
length(unique(df$CB4.Symbol.sq.))
length(unique(df$CD8.Symbol.sq.))
length(unique(df$Msc.Symbol.sq.))

## Figure 2E
#Load the following objects 
Counts_Focal Adhesion <- read.delim("Counts_Focal Adhesion.txt", h=T) 
Counts_HIPPO <- read.delim("Counts_HIPPO.txt", h=T)
Counts_IGF1 <- read.delim("Counts_IGF1.txt", h=T)
Counts_Proliferative <- read.delim("Counts_Proliferative.txt", h=T)
Counts_Glyco <- read.delim("Counts_Glyco.txt", h=T)
Counts_G protein signaling <- read.delim("Counts_G protein signaling.txt", h=T)
Counts_Parkin-Ubiquitin <- read.delim("Counts_Parkin-Ubiquitin.txt", h=T)
Counts_OP <- read.delim("Counts_OP.txt", h=T)
Counts_BMP <- read.delim("Counts_BMP.txt", h=T)
Counts_TF <- read.delim("Counts_TF.txt", h=T)
Counts_WNT <- read.delim("Counts_WNT.txt", h=T)


AD10_markers <- read.delim("AD10.markers.txt", h=T) 
DD8_markers <- read.delim("DD8.markers.txt", h=T) 
CB4_markers <- read.delim("CB4.markers.txt", h=T) 
CD8_markers <- read.delim("CD8.markers.txt", h=T) 
MSC_markers <- read.delim("Msc.markers.txt", h=T)

x<-subset(AD10_markers, gene %in% Counts_Focal Adhesion$Gene)
x0<-x[x$cluster== 'Focal adhesion signaling  ',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"

x<-subset(DD8_markers, gene %in% Counts_Focal Adhesion$Gene)
x1<-x[x$cluster== 'Focal adhesion signaling',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_Focal Adhesion$Gene)
x2<-x[x$cluster== 'Focal adhesion signaling ',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_Focal Adhesion$Gene)
x3<-x[x$cluster== 'Focal adhesion signaling ',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"

x<-subset(MSC_markers, gene %in% Counts_Focal Adhesion$Gene)
x4<-x[x$cluster== 'Focal adhesion signaling',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_Focal Adhesion$Gene)
x5<-x[x$cluster== '0',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

Focal Adhesion_0<- merge(Counts_Focal Adhesion,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
Focal Adhesion_1<- merge(Focal Adhesion_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
Focal Adhesion_2<- merge(Focal Adhesion_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
Focal Adhesion_3<- merge(Focal Adhesion_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
Focal Adhesion_4<- merge(Focal Adhesion_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
Focal Adhesion_5<- merge(Focal Adhesion_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


Focal Adhesion<- Focal Adhesion_5[,c(1,2,4,10,16,22,28,34)]

Focal Adhesion[is.na(Focal Adhesion)] <- 0

x<-subset(AD10_markers, gene %in% Counts_HIPPO$Gene)
x0<-x[x$cluster== 'Pathways Regulating Hippo Signaling',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_HIPPO$Gene)
x1<-x[x$cluster== 'Pathways Regulating Hippo Signaling',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_HIPPO$Gene)
x2<-x[x$cluster== 'Pathways Regulating Hippo Signaling',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_HIPPO$Gene)
x3<-x[x$cluster== 'Pathways Regulating Hippo Signaling',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_HIPPO$Gene)
x4<-x[x$cluster== 'Pathways Regulating Hippo Signaling',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_HIPPO$Gene)
x5<-x[x$cluster== '1',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

HIPPO_0<- merge(Counts_HIPPO,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
HIPPO_1<- merge(HIPPO_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
HIPPO_2<- merge(HIPPO_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
HIPPO_3<- merge(HIPPO_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
HIPPO_4<- merge(HIPPO_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
HIPPO_5<- merge(HIPPO_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


HIPPO <- HIPPO_5[,c(1,2,4,10,16,22,28,34)]

HIPPO[is.na(HIPPO)] <- 0



x<-subset(AD10_markers, gene %in% Counts_IGF1$Gene)
x0<-x[x$cluster== 'IGF1−Akt signaling',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_IGF1$Gene)
x1<-x[x$cluster== 'IGF1−Akt signaling',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_IGF1$Gene)
x2<-x[x$cluster== 'IGF1−Akt signaling',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_IGF1$Gene)
x3<-x[x$cluster== 'IGF1−Akt signaling',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_IGF1$Gene)
x4<-x[x$cluster== 'IGF1−Akt signaling',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_IGF1$Gene)
x5<-x[x$cluster== '2',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

IGF1_0<- merge(Counts_IGF1,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
IGF1_1<- merge(IGF1_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
IGF1_2<- merge(IGF1_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
IGF1_3<- merge(IGF1_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
IGF1_4<- merge(IGF1_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
IGF1_5<- merge(IGF1_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


IGF1 <- IGF1_5[,c(1,2,4,10,16,22,28,34)]

IGF1[is.na(IGF1)] <- 0


x<-subset(AD10_markers, gene %in% Counts_Proliferative$Gene)
x0<-x[x$cluster== 'Proliferative',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_Proliferative$Gene)
x1<-x[x$cluster== 'Proliferative',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_Proliferative$Gene)
x2<-x[x$cluster== 'Proliferative',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_Proliferative$Gene)
x3<-x[x$cluster== 'Proliferative',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_Proliferative$Gene)
x4<-x[x$cluster== 'Proliferative',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_Proliferative$Gene)
x5<-x[x$cluster== '3',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

Proliferative_0<- merge(Counts_Proliferative,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
Proliferative_1<- merge(Proliferative_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
Proliferative_2<- merge(Proliferative_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
Proliferative_3<- merge(Proliferative_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
Proliferative_4<- merge(Proliferative_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
Proliferative_5<- merge(Proliferative_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


Proliferative <- Proliferative_5[,c(1,2,4,10,16,22,28,34)]

Proliferative[is.na(Proliferative)] <- 0


x<-subset(AD10_markers, gene %in% Counts_Glyco$Gene)
x0<-x[x$cluster== 'Glycolysis and Gluconeogenesis',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_Glyco$Gene)
x1<-x[x$cluster== 'Glycolysis and Gluconeogenesis',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_Glyco$Gene)
x2<-x[x$cluster== 'Glycolysis and Gluconeogenesis',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_Glyco$Gene)
x3<-x[x$cluster== 'Glycolysis and Gluconeogenesis',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_Glyco$Gene)
x4<-x[x$cluster== 'Glycolysis and Gluconeogenesis',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_Glyco$Gene)
x5<-x[x$cluster== '4',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

Glyco_0<- merge(Counts_Glyco,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
Glyco_1<- merge(Glyco_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
Glyco_2<- merge(Glyco_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
Glyco_3<- merge(Glyco_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
Glyco_4<- merge(Glyco_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
Glyco_5<- merge(Glyco_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


Glyco <- Glyco_5[,c(1,2,4,10,16,22,28,34)]

Glyco[is.na(Glyco)] <- 0


x<-subset(AD10_markers, gene %in% Counts_G protein signaling$Gene)
x0<-x[x$cluster== 'G protein signaling',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_G protein signaling$Gene)
x1<-x[x$cluster== 'G protein signaling',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_G protein signaling$Gene)
x2<-x[x$cluster== 'G protein signaling',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_G protein signaling$Gene)
x3<-x[x$cluster== 'G protein signaling',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_G protein signaling$Gene)
x4<-x[x$cluster== 'G protein signaling',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_G proteinn signaling$Gene)
x5<-x[x$cluster== '5',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

G protein signaling_0<- merge(Counts_G protein signaling,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
G protein signaling_1<- merge(G protein signaling_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
G protein signaling_2<- merge(G protein signaling_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
G protein signaling_3<- merge(G protein signaling_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
G protein signaling_4<- merge(G protein signaling_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
G protein signaling_5<- merge(G protein signaling_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


G protein signaling <- G protein signaling_5[,c(1,2,4,10,16,22,28,34)]

G protein signaling[is.na(G protein signaling)] <- 0


x<-subset(AD10_markers, gene %in% Counts_Parkin-Ubiquitin$Gene)
x0<-x[x$cluster== 'Parkin−Ubiquitin proteasomal system',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_Parkin-Ubiquitin$Gene)
x1<-x[x$cluster== 'Parkin−Ubiquitin proteasomal system',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_Parkin-Ubiquitin$Gene)
x2<-x[x$cluster== 'Parkin−Ubiquitin proteasomal system',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_Parkin-Ubiquitin$Gene)
x3<-x[x$cluster== 'Parkin−Ubiquitin proteasomal system',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_Parkin-Ubiquitin$Gene)
x4<-x[x$cluster== 'Parkin−Ubiquitin proteasomal system',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_Parkin-Ubiquitin$Gene)
x5<-x[x$cluster== '6',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

Parkin-Ubiquitin_0<- merge(Counts_Parkin-Ubiquitin,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
Parkin-Ubiquitin_1<- merge(Parkin-Ubiquitin_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
Parkin-Ubiquitin_2<- merge(Parkin-Ubiquitin_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
Parkin-Ubiquitin_3<- merge(Parkin-Ubiquitin_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
Parkin-Ubiquitin_4<- merge(Parkin-Ubiquitin_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
Parkin-Ubiquitin_5<- merge(Parkin-Ubiquitin_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


Parkin-Ubiquitin <- Parkin-Ubiquitin_5[,c(1,2,4,10,16,22,28,34)]

Parkin-Ubiquitin[is.na(Parkin-Ubiquitin)] <- 0


x<-subset(AD10_markers, gene %in% Counts_OP$Gene)
x0<-x[x$cluster== 'Oxidative phosphorylation',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_OP$Gene)
x1<-x[x$cluster== 'Oxidative phosphorylation',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_OP$Gene)
x2<-x[x$cluster== 'Oxidative phosphorylation',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_OP$Gene)
x3<-x[x$cluster== 'Oxidative phosphorylation',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_OP$Gene)
x4<-x[x$cluster== 'Oxidative phosphorylation',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_OP$Gene)
x5<-x[x$cluster== '7',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

OP_0<- merge(Counts_OP,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
OP_1<- merge(OP_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
OP_2<- merge(OP_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
OP_3<- merge(OP_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
OP_4<- merge(OP_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
OP_5<- merge(OP_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


OP <- OP_5[,c(1,2,4,10,16,22,28,34)]

OP[is.na(OP)] <- 0


x<-subset(AD10_markers, gene %in% Counts_BMP$Gene)
x0<-x[x$cluster== 'BMP signalling',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_BMP$Gene)
x1<-x[x$cluster== 'BMP signalling',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_BMP$Gene)
x2<-x[x$cluster== 'BMP signalling',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_BMP$Gene)
x3<-x[x$cluster== 'BMP signalling',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_BMP$Gene)
x4<-x[x$cluster== 'BMP signalling',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_BMP$Gene)
x5<-x[x$cluster== '8',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

BMP_0<- merge(Counts_BMP,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
BMP_1<- merge(BMP_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
BMP_2<- merge(BMP_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
BMP_3<- merge(BMP_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
BMP_4<- merge(BMP_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
BMP_5<- merge(BMP_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


BMP <- BMP_5[,c(1,2,4,10,16,22,28,34)]

BMP[is.na(BMP)] <- 0

x<-subset(AD10_markers, gene %in% Counts_TF$Gene)
x0<-x[x$cluster== 'Translation Factors',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_TF$Gene)
x1<-x[x$cluster== 'Translation Factors',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_TF$Gene)
x2<-x[x$cluster== 'Translation Factors',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_TF$Gene)
x3<-x[x$cluster== 'Translation Factors',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_TF$Gene)
x4<-x[x$cluster== 'Translation Factors',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_TF$Gene)
x5<-x[x$cluster== '9',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

TF_0<- merge(Counts_TF,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
TF_1<- merge(TF_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
TF_2<- merge(TF_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
TF_3<- merge(TF_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
TF_4<- merge(TF_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
TF_5<- merge(TF_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


TF <- TF_5[,c(1,2,4,10,16,22,28,34)]

TF[is.na(TF)] <- 0


x<-subset(AD10_markers, gene %in% Counts_WNT$Gene)
x0<-x[x$cluster== 'Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds',]
colnames(x0)[2] <- "logFC_AD10"
colnames(x0)[7] <- "Gene"


x<-subset(DD8_markers, gene %in% Counts_WNT$Gene)
x1<-x[x$cluster== 'Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds',]
colnames(x1)[2] <- "logFC_DD8"
colnames(x1)[7] <- "Gene"

x<-subset(CB4_markers, gene %in% Counts_WNT$Gene)
x2<-x[x$cluster== 'Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds',]
colnames(x2)[2] <- "logFC_CB4"
colnames(x2)[7] <- "Gene"

x<-subset(CD8_markers, gene %in% Counts_WNT$Gene)
x3<-x[x$cluster== 'Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds',]
colnames(x3)[2] <- "logFC_CD8"
colnames(x3)[7] <- "Gene"


x<-subset(MSC_markers, gene %in% Counts_WNT$Gene)
x4<-x[x$cluster== 'Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds',]
colnames(x4)[2] <- "logFC_Msc"
colnames(x4)[7] <- "Gene"

x<-subset(Tert_markers, gene %in% Counts_WNT$Gene)
x5<-x[x$cluster== '10',]
colnames(x5)[2] <- "logFC_ALL"
colnames(x5)[7] <- "Gene"

WNT_0<- merge(Counts_WNT,x0 ,by="Gene",all.x=TRUE, all.y=TRUE )
WNT_1<- merge(WNT_0,x1 ,by="Gene",all.x=TRUE, all.y=TRUE )
WNT_2<- merge(WNT_1,x2 ,by="Gene", all.x=TRUE, all.y=TRUE)
WNT_3<- merge(WNT_2,x3 ,by="Gene", all.x=TRUE, all.y=TRUE)
WNT_4<- merge(WNT_3, x4, by="Gene", all.x=TRUE, all.y=TRUE)
WNT_5<- merge(WNT_4, x5, by="Gene", all.x=TRUE, all.y=TRUE)


WNT <- WNT_5[,c(1,2,4,10,16,22,28,34)]

WNT[is.na(WNT)] <- 0

#Dot plot of marker genes across the cell lines
Focaladhesion_freq <- as.data.frame(table(Focal Adhesion$Count)/nrow(Focal Adhesion)*100)
Focaladhesion_freq$Cluster <- "Focal adhesion"
HIPPO_freq <- as.data.frame(table(HIPPO$Count)/nrow(HIPPO)*100)
HIPPO_freq$Cluster<- "HIPPO"
IGF1_freq <- as.data.frame(table(IGF1$Count)/nrow(IGF1)*100)
IGF1_freq$Cluster<- "IGF1"
Prolif_freq<- as.data.frame(table(Proliferative$Count)/nrow(Proliferative)*100)
Prolif_freq$Cluster <-"Proliferation"
Glyco_freq<- as.data.frame(table(Glyco$Count)/nrow(Glyco)*100)
Glyco_freq$Cluster<- "Glycogenesis"
Gprotein_freq<- as.data.frame(table(G protien signaling$Count)/nrow(G protien signaling)*100)
Gprotein_freq$Cluster<- "G protein signaling"
Ubiquitin_freq <- as.data.frame(table(Parkin-Ubiquitin$Count)/nrow(Parkin-Ubiquitin)*100)
Ubiquitin_freq$Cluster <-"Ubiquitin"
Oxhpos_freq <- as.data.frame(table(OP$Count)/nrow(OP)*100)
Oxhpos_freq$Cluster <-"Oxphos"
BMP_freq <- as.data.frame(table(BMP$Count)/nrow(BMP)*100)
BMP_freq$Cluster <-"BMP"
TF_freq<- as.data.frame(table(TF$Count)/nrow(TF)*100)
TF_freq$Cluster<-"Translation G protien signalingctor"
WNT_freq<- as.data.frame(table(WNT$Count)/nrow(WNT)*100)
WNT_freq$Cluster<- "WNT"
Freq_data_all <- rbind(Focaladhesion_freq, HIPPO_freq, IGF1_freq, Prolif_freq, Glyco_freq, Gprotein_freq,
                       Ubiquitin_freq, Oxhpos_freq, BMP_freq, TF_freq, WNT_freq)

Freq_data_all$Cluster <- G protien signalingctor(Freq_data_all$Cluster, levels = c("Focal adhesion","HIPPO","IGF1","Proliferation","Glycogenesis",
                                                          "G protein signaling","Ubiquitin","Oxphos","BMP","Translation G protien signalingctor",
                                                          "WNT"))


ggplot(Freq_data_all, aes(x = G protien signalingctor(Cluster),
               y = Var1,
               size = Freq)) +
  geom_point(shape = 21, fill = "pink", color = "black", alpha = 0.8) +
  scale_size(range = c(1, 15), name = "Freq (%)") +
  labs(x = "Cluster", y = "Gene count") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90")
  )

##Figure 2F
#Dotplot of marker genes per each cell line
Focal Adhesion <- Focal Adhesion[order(-Focal Adhesion$Count),]
Focal Adhesion_common<-Focal Adhesion[c(1:2),]
Focal Adhesion_common$Cluster<- "Focal Adhesion"

HIPPO <- HIPPO[order(-HIPPO$Count),]
HIPPO_common<-HIPPO[C(1:63),]
HIPPO_common$Cluster<- "HIPPO"

Common<- merge(Focal Adhesion_common, HIPPO_common, all.x=TRUE, all.y=TRUE)

IGF1 <- IGF1[order(-IGF1$Count),]
IGF1_common<-IGF1[C(1:2),]
IGF1_common$Cluster<- "IGF1"

Common<- merge(IGF1_common, Common, all.x=TRUE, all.y=TRUE)

Proliferative <- Proliferative[order(-Proliferative$Count),]
Proliferative_common<-Proliferative[c(1:14),]
Proliferative_common$Cluster<- "Proliferative"

Common<- merge(Proliferative_common, Common, all.x=TRUE, all.y=TRUE)

Glyco <- Glyco[order(-Glyco$Count),]
Glyco_common<-Glyco[c(1:13),]
Glyco_common$Cluster<- "Glyco"

Common<- merge(Glyco_common, Common, all.x=TRUE, all.y=TRUE)

G protien signaling <- G protien signaling[order(-G protien signaling$Count),]
G protien signaling_common<-G protien signaling[c(1:171),]
G protien signaling_common$Cluster<- "G protien signaling"

Common<- merge(G protien signaling_common, Common, all.x=TRUE, all.y=TRUE)

Parkin-Ubiquitin <- Parkin-Ubiquitin[order(-Parkin-Ubiquitin$Count),]
Parkin-Ubiquitin_common<-Parkin-Ubiquitin[c(1:6),]
Parkin-Ubiquitin_common$Cluster<- "Parkin-Ubiquitin"

Common<- merge(Parkin-Ubiquitin_common, Common, all.x=TRUE, all.y=TRUE)

OP <- OP[order(-OP$Count),]
OP_common<-OP[c(1:3),]
OP_common$Cluster<- "Oxphos"

Common<- merge(OP_common, Common, all.x=TRUE, all.y=TRUE)

BMP <- BMP[order(-BMP$Count),]
BMP_common<-BMP[c(1:8),]
BMP_common$Cluster<- "BMP"

Common<- merge(BMP_common, Common, all.x=TRUE, all.y=TRUE)

TF <- TF[order(-TF$Count),]
TF_common<-TF[c(1:238),]
TF_common$Cluster<- "Translation G protien signalingctors"

Common<- merge(TF_common, Common, all.x=TRUE, all.y=TRUE)

WNT <- WNT[order(-WNT$Count),]
WNT_common<-WNT[1,]
WNT_common$Cluster<- "WNT"

Common<- merge(WNT_common, Common, all.x=TRUE, all.y=TRUE)


library(dplyr)
library(ggplot2)

genes_of_interest <- c("COL1A1","ADH1B","ASAP1","CENPK","EEF1A1P12", "ACTA2","CD63","GDF5","BNC2","FLCN","AC013652.1")

dot_AD10 <- Common %>%
    filter(Gene %in% genes_of_interest) %>%        
    mutate(
   Gene = G protien signalingctor(Gene, levels = genes_of_interest)
      )

dot_AD10$Cluster <- G protien signalingctor(dot_AD10$Cluster, levels = rev(c("Focal Adhesion","HIPPO","IGF1","Proliferative","Glyco",
                                                    "G protien signaling","Parkin-Ubiquitin","Oxphos","BMP","Translation G protien signalingctors",
                                                    "WNT")))

library(ggplot2)
library(RColorBrewer)

ggplot(dot_AD10, aes(
  x = Gene,
  y = Cluster,
  color = logFC_AD10    
)) +
  geom_point(size = 8, alpha = 0.9) +
  scale_color_gradientn(
    colors = c(
      "white",  # very light pink
      "#fbc4e2",
      "#f768a1",
      "#dd3497",
      "#ae017e",
      "#7a0177",
      "#49006a"   # dark red/purple
    ),
    limits = c(0, 5),
    oob = scales::squish,
    name = "logFC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "AD10",
    x = "Gene",
    y = "Cluster"
  )


##DD8
dot_DD8 <- Common %>%
  filter(Gene %in% genes_of_interest) %>%        
  mutate(
    Gene = G protien signalingctor(Gene, levels = genes_of_interest)
  )

dot_DD8$Cluster <- G protien signalingctor(dot_DD8$Cluster, levels = rev(c("Focal Adhesion","HIPPO","IGF1","Proliferative","Glyco",
                                                            "G protien signaling","Parkin-Ubiquitin","Oxphos","BMP","Translation G protien signalingctors",
                                                            "WNT")))

library(ggplot2)
library(RColorBrewer)


ggplot(dot_DD8, aes(
  x = Gene,
  y = Cluster,
  color = logFC_DD8    
)) +
  geom_point(size = 8, alpha = 0.9) +
  scale_color_gradientn(
    colors = c(
      "white",  # very light pink
      "#fbc4e2",
      "#f768a1",
      "#dd3497",
      "#ae017e",
      "#7a0177",
      "#49006a"   # dark red/purple
    ),
    limits = c(0, 5),
    oob = scales::squish,
    name = "logFC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "DD8",
    x = "Gene",
    y = "Cluster"
  )

##CB4
dot_CB4 <- Common %>%
  filter(Gene %in% genes_of_interest) %>%        
  mutate(
    Gene = G protien signalingctor(Gene, levels = genes_of_interest)
  )

dot_CB4$Cluster <- G protien signalingctor(dot_CB4$Cluster, levels = rev(c("Focal Adhesion","HIPPO","IGF1","Proliferative","Glyco",
                                                           "G protien signaling","Parkin-Ubiquitin","Oxphos","BMP","Translation G protien signalingctors",
                                                           "WNT")))

library(ggplot2)
library(RColorBrewer)


ggplot(dot_CB4, aes(
  x = Gene,
  y = Cluster,
  color = logFC_CB4    
)) +
  geom_point(size = 8, alpha = 0.9) +
  scale_color_gradientn(
    colors = c(
      "white",  # very light pink
      "#fbc4e2",
      "#f768a1",
      "#dd3497",
      "#ae017e",
      "#7a0177",
      "#49006a"   # dark red/purple
    ),
    limits = c(0, 5),
    oob = scales::squish,
    name = "logFC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "CB4",
    x = "Gene",
    y = "Cluster"
  )

##CD8
dot_CD8 <- Common %>%
  filter(Gene %in% genes_of_interest) %>%        
  mutate(
    Gene = G protien signalingctor(Gene, levels = genes_of_interest)
  )

dot_CD8$Cluster <- G protien signalingctor(dot_CD8$Cluster, levels = rev(c("Focal Adhesion","HIPPO","IGF1","Proliferative","Glyco",
                                                           "G protien signaling","Parkin-Ubiquitin","Oxphos","BMP","Translation G protien signalingctors",
                                                           "WNT")))

library(ggplot2)
library(RColorBrewer)



ggplot(dot_CD8, aes(
  x = Gene,
  y = Cluster,
  color = logFC_CD8    
)) +
  geom_point(size = 8, alpha = 0.9) +
  scale_color_gradientn(
    colors = c(
      "white",  # very light pink
      "#fbc4e2",
      "#f768a1",
      "#dd3497",
      "#ae017e",
      "#7a0177",
      "#49006a"   # dark red/purple
    ),
    limits = c(0, 4),
    oob = scales::squish,
    name = "logFC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "CD8",
    x = "Gene",
    y = "Cluster"
  )

# Final aesthetics were done in Illustrator