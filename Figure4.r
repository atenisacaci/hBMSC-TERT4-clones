#Figure 4
#Load the following objects from https://osf.io/wxpgn/
All_genes_clones <- readRDS("All_genes_clones.rds")
Result <- readRDS("Clone_Cluster_markers.rds")

#Figure 4A
Heterogenity_data <- readRDS("MarkerGeneList.Rds")
AllGene_heterogenity <- readRDS("AllGene.Rds")


All_common <- All_genes_clones[All_genes_clones %in% AllGene_heterogenity]

Heterogenity_markers <- bind_rows(Heterogenity_data, .id = "column_label")
Heterogenity_markers <- Heterogenity_markers[Heterogenity_markers$Marker=="Enriched",]
Clusters <- c("Cl_0","Cl_1","Cl_2","Cl_3","Cl_4")

# Reduce both marker lists to the genes detected in both experiments
Heterogenity_markers <- Heterogenity_markers[Heterogenity_markers$Gene %in% All_common,]
Result <- Result[Result$Gene %in% All_common,]

Gene_groups_1 <- list()
for (i in 1: length(Clusters)){
  Gene_groups_1[[i]] <- Heterogenity_markers[Heterogenity_markers$column_label == Clusters[i],'Gene']
}
names(Gene_groups_1) <- Clusters

markers<-Result
Gene_groups <- list()
Gene_groups[[1]] <- markers[markers$column_label== "Cl_0" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[2]] <- markers[markers$column_label== "Cl_1" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[3]] <- markers[markers$column_label== "Cl_2" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[4]] <- markers[markers$column_label== "Cl_3" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[5]] <- markers[markers$column_label== "Cl_4" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[6]] <- markers[markers$column_label== "Cl_5" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[7]] <- markers[markers$column_label== "Cl_6" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[8]] <- markers[markers$column_label== "Cl_7" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[9]] <- markers[markers$column_label== "Cl_8" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[10]] <- markers[markers$column_label== "Cl_9" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[11]] <- markers[markers$column_label== "Cl_10" & markers$Marker =="Enriched",'Gene' ]
names(Gene_groups) <- c('Cl_0','Cl_1','Cl_2','Cl_3','Cl_4','Cl_5','Cl_6','Cl_7','Cl_8','Cl_9','Cl_10')

# Calculate enrichment based on a hypergeometric test

mat <- matrix(NA, ncol=length(Gene_groups_1), nrow=length(Gene_groups))
colnames(mat) <-names(Gene_groups_1)
rownames(mat) <-names(Gene_groups)

for(i in 1:length(Gene_groups_1)){
  for(k in 1:length(Gene_groups)){
    tmp_i <- Gene_groups_1[[i]]
    tmp_k <- Gene_groups[[k]]
    mat[k,i] <- phyper(length(tmp_i[tmp_i %in% tmp_k]), length(tmp_k), length(All_common[!All_common %in% tmp_k]), length(tmp_i), lower.tail = F)
  }
}

mat <- apply(mat,2,p.adjust,method = "BH")
mat <- -log10(mat)

#plot as heatmap
library(fields)
library(scales)
library(gplots)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(mat),length=51))
heatmap.2(mat,Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=rownames(mat),labCol=colnames(mat) )


#Figure 4B
###Relation of cluster specific genes to day7 adipocyte and osteoblast differentiation
Clones <- read.delim("Clones_result.txt",h=T)
Clones <- Clones[complete.cases(Clones), ]

# Alignment with expression data (https://doi.org/10.1038/s41588-019-0359-1), hypergeometric test to test overlap with osteogenic and adiopgenic induced genes 
TERT <- read.delim("final.TERT.txt",h=T)
TERT$SYMBOL <- sub("\\|.*", "", TERT$Annotation.Divergence)

library(dplyr)
library(fields)
library(gplots)

# Combine markers of single RNA-seq in one data frame and split enriched markers in a list
markers <- bind_rows(Result, .id = "column_label")
markers<-Result

All_genes_clones <- All_genes_clones[All_genes_clones %in% TERT$SYMBOL]
TERT <- TERT[TERT$SYMBOL %in% All_genes_clones,]

markers <- markers[markers$Gene %in% All_genes_clones,]


Gene_groups <- list()
Gene_groups[[1]] <- markers[markers$column_label== "Cl_0" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[2]] <- markers[markers$column_label== "Cl_1" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[3]] <- markers[markers$column_label== "Cl_2" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[4]] <- markers[markers$column_label== "Cl_3" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[5]] <- markers[markers$column_label== "Cl_4" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[6]] <- markers[markers$column_label== "Cl_5" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[7]] <- markers[markers$column_label== "Cl_6" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[8]] <- markers[markers$column_label== "Cl_7" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[9]] <- markers[markers$column_label== "Cl_8" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[10]] <- markers[markers$column_label== "Cl_9" & markers$Marker =="Enriched",'Gene' ]
Gene_groups[[11]] <- markers[markers$column_label== "Cl_10" & markers$Marker =="Enriched",'Gene' ]


names(Gene_groups) <- c('Cl_0','Cl_1','Cl_2','Cl_3','Cl_4','Cl_5','Cl_6','Cl_7','Cl_8','Cl_9','Cl_10')

# Group the osteoblast and adipocyte TERT genes in a list
Gene_groups_2 <- list()
Gene_groups_2[[1]] <- TERT[TERT$pVal_RNA_Ob7d_vs_Msc0h < 0.01 & TERT$log2FC_Ob7d < 0 ,'SYMBOL']
Gene_groups_2[[2]] <- TERT[TERT$pVal_RNA_Ob7d_vs_Msc0h < 0.01 & TERT$log2FC_Ob7d > 0 ,'SYMBOL']
Gene_groups_2[[3]] <- TERT[TERT$pVal_RNA_Ad7d_vs_Msc0h < 0.01 & TERT$log2FC_Ad7d < 0 ,'SYMBOL']
Gene_groups_2[[4]] <- TERT[TERT$pVal_RNA_Ad7d_vs_Msc0h < 0.01 & TERT$log2FC_Ad7d > 0 ,'SYMBOL']
names(Gene_groups_2) <- c('Ob_down','Ob_up','Ad_down','Ad_up')
# Test the overlap of both gene groups using a hypergeometric test
mat <- matrix(NA, ncol=length(Gene_groups),nrow=length(Gene_groups_2))
colnames(mat) <- names(Gene_groups)
rownames(mat) <- names(Gene_groups_2)

for (i in 1:length(Gene_groups_2)){
  for (k in 1:length(Gene_groups)){
    tmp_i <- Gene_groups_2[[i]]
    tmp_k <- Gene_groups[[k]]
    x <- length(tmp_k[tmp_k %in% tmp_i])
    m <- length(tmp_k)
    n <- length(Atenisa_all[!Atenisa_all %in% tmp_k])
    y <- length(tmp_i)
    mat[i,k] <- phyper(x,m,n,y,lower.tail = F)
  }
}
# transform to log scale
mat <- -log10(mat)
mat[mat=="Inf"]<-187
# Show enrichment in a heatmap
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(mat),length=51))
heatmap.2(mat,Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none')
rm(x,m,n,y,i,k,mat_col,mat_col_breaks,tmp_i,tmp_k,tmp_all,mat,markers, Result,Gene_groups,Gene_groups_2, All )


#Figure 4C
###Checking the enriched markers genes for eBMD associated SNPs in the vicinity of the TSS
# GWAS summary statistics from http://www.gefos.org/?q=content/data-release-2018 doi:10.1038/s41588-018-0302-x
GWAS_summary <- fread("Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",h=T)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(gplots)
library(fields)
library(S4Vectors)
# Define TSS of genes from RNA-seq on hg19 coordinates
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Get transcripts for hg19
GR <- data.frame(transcripts(txdb))
# Get symbols to UCSC IDs
GR2 <- select(Homo.sapiens,GR$tx_name, "SYMBOL","TXNAME")
GR <- merge(GR,GR2, by.x="tx_name", by.y="TXNAME")
rm(GR2)
GR <- GR[GR$seqnames %in% paste('chr',c(1:23,'X','Y','M'),sep=""),]
# Reduce to the gene expressed in the scRNAseq dataset
GR <- unique(GR[GR$SYMBOL %in% All_genes_clones,])
# Get lowest TSS for "+"-stranded
GR1 <- GR[GR$strand =="+",] 
GR1 <- GR1[order(GR1$SYMBOL,GR1$start),]
GR1$dup <- duplicated(GR1$SYMBOL)
GR1 <- GR1[GR1$dup =="FALSE",]
GR1$dup <- NULL
GR1 <- GR1[,c('SYMBOL','seqnames','start','strand')]
names(GR1) <- c("Symbol",'Chr',"TSS","Strand")
# Get highest TSS for "-"-stranded
GR2 <- GR[GR$strand =="-",] 
GR2 <- GR2[order(GR2$SYMBOL,-GR2$end),]
GR2$dup <- duplicated(GR2$SYMBOL)
GR2 <- GR2[GR2$dup =="FALSE",]
GR2$dup <- NULL
GR2 <- GR2[,c('SYMBOL','seqnames','end','strand')]
names(GR2) <- c("Symbol",'Chr',"TSS","Strand")
GR <- rbind(GR1,GR2)
rm(GR1,GR2, txdb)

# Reformat GWAS summary statistics
GWAS_summary <- GWAS_summary[,c('RSID','CHR','BP','P.NI','BETA')]
names(GWAS_summary)[1] <- 'SNP'
GWAS_summary$Chr <- paste("chr",GWAS_summary$CHR, sep="")
GWAS_summary$Start <- as.numeric(GWAS_summary$BP)
GWAS_summary <- GWAS_summary[complete.cases(GWAS_summary$Start),c('SNP','P.NI','Chr','Start','BETA')]
GWAS_summary <- data.frame(GWAS_summary)
GWAS_summary$Pval <- as.numeric(GWAS_summary$P.NI)

# Do overlap of SNPs in window of 5 Mb
grGenes <- with(unique(GR[,c("Symbol","Chr","TSS")]) , GRanges(Chr, IRanges(start=TSS - 5000000, end=TSS + 5000000, names=Symbol)))
grSummary <- with(unique(GWAS_summary[,c("SNP","Chr","Start")]) , GRanges(Chr, IRanges(start=Start, end=Start, names=SNP)))

hits = findOverlaps(grGenes,grSummary)
tmp2 <- cbind(data.frame(ranges(grGenes)[queryHits(hits)]),data.frame(ranges(grSummary)[subjectHits(hits)]))
colnames(tmp2) <- c('TSSminus','TSSplus','Window','Symbol','SNP_Start','SNP_End','SNP_Length','SNP')
head(tmp2)

# Calculate distance bewteen TSS and SNP
tmp2$Distance <- tmp2$TSSminus+5000000 - tmp2$SNP_Start
tmp2 <- tmp2[,c('Symbol','SNP','Distance')]
tmp2 <- merge(tmp2[,c('Symbol','SNP','Distance')], GWAS_summary[,c("SNP","Pval","BETA")], by="SNP")
tmp2 <- merge( GR,tmp2[,c('Symbol','SNP','Pval','Distance')], by="Symbol")

saveRDS(tmp2,file="SNPs_Clones.rds")
#Chisquare test for enrichment with different distances from the TSS of genes from the scRNA-seq clusters (focus on enriched genes extracted from the list "Result")
mat2 <- matrix(NA, ncol = 11, nrow = 13)
x <- 1
clusters <- paste0("Cl_", 0:10)

for (k in c(50000, 100000, 250000, seq(500000, 5000000, length = 10))) {
  for (i in seq_along(clusters)) {
    
    tmp <- Result[
      Result$column_label == clusters[i] &
        Result$Marker == "Enriched",
      "Gene"
    ]
    
    a <- length(tmp2[tmp2$Symbol %in% tmp & tmp2$Pval < 5E-8 & abs(tmp2$Distance) < k, "SNP"])
    b <- length(tmp2[tmp2$Symbol %in% tmp & tmp2$Pval > 5E-8 & abs(tmp2$Distance) < k, "SNP"])
    c <- length(tmp2[tmp2$Pval < 5E-8 & abs(tmp2$Distance) < k, "SNP"])
    d <- length(tmp2[tmp2$Pval > 5E-8 & abs(tmp2$Distance) < k, "SNP"])
    
    M <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    dimnames(M) <- list(
      group = c("in_cluster", "in_gwascat"),
      categ = c("HBMD_assoc", "not_HBMD_assoc")
    )
    
    mat2[x, i] <- chisq.test(M)$p.value
  }
  x <- x + 1
}

rownames(mat2) <- paste0("Dist_", c(50000, 100000, 250000, seq(500000, 5000000, length = 10)))
colnames(mat2) <- clusters

mat2 <- -log10(mat2[1:7,])
mat2[mat2=="Inf"]<-300
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(5,max(mat2),length=51))
heatmap.2(as.matrix(mat2),Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=rownames(mat2),labCol=colnames(mat2) )

rm(mat2,a,b,c,d,tmp,k,i,GR,mat_col, mat_col_breaks, GWAS_summary,hits,grGenes, grSummary,Result, All_genes_clones)

#Figure 4D
#Relation of cluster specific genes to osteoporotic and healthy subjects
Result <- readRDS("Clone_Cluster_markers.rds")
All_genes_clones <- readRDS("All_genes_clones.rds")
colnames(Result)[colnames(Result) == 'Gene'] <- 'Symbol'
data_Iliac <- read.delim("ReadyToUse_E-MEXP-1618.txt",h=T)

# Combine markers of single RNA-seq in one data frame and split enriched markers in a list
markers<-Result

All_genes_clones <- All_genes_clones[All_genes_clones %in% data_Iliac$SYMBOL]
data_Iliac <- data_Iliac[data_Iliac$SYMBOL %in% All_genes_clones,]

markers <- markers[markers$Symbol %in% All_genes_clones,]

Gene_groups_1 <- list()
Gene_groups_1[[1]] <- unique(data_Iliac[data_Iliac$Pval_Osteoporotic_Healthy < 0.005 & data_Iliac$LogFC_Osteoporotic_Healthy < 0,"SYMBOL"])
Gene_groups_1[[2]] <- unique(data_Iliac[data_Iliac$Pval_Osteoporotic_Healthy < 0.005 & data_Iliac$LogFC_Osteoporotic_Healthy > 0,"SYMBOL"])
names(Gene_groups_1) <- c('Iliac_down','Iliac_up')

Gene_groups_2 <- list()
Gene_groups_2[[1]] <- markers[markers$column_label== "Cl_0" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[2]] <- markers[markers$column_label== "Cl_1" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[3]] <- markers[markers$column_label== "Cl_2" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[4]] <- markers[markers$column_label== "Cl_3" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[5]] <- markers[markers$column_label== "Cl_4" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[6]] <- markers[markers$column_label== "Cl_5" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[7]] <- markers[markers$column_label== "Cl_6" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[8]] <- markers[markers$column_label== "Cl_7" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[9]] <- markers[markers$column_label== "Cl_8" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[10]] <- markers[markers$column_label== "Cl_9" & markers$Marker =="Enriched",'Symbol' ]
Gene_groups_2[[11]] <- markers[markers$column_label== "Cl_10" & markers$Marker =="Enriched",'Symbol' ]


names(Gene_groups_2) <- c('Cl_0','Cl_1','Cl_2','Cl_3','Cl_4','Cl_5','Cl_6','Cl_7','Cl_8','Cl_9','Cl_10')

# Calculate enrichment based on a hypergeometric test
mat <- matrix(NA, ncol=length(Gene_groups_1), nrow=length(Gene_groups_2))
colnames(mat) <-names(Gene_groups_1)
rownames(mat) <-names(Gene_groups_2)

for(i in 1:length(Gene_groups_1)){
  for(k in 1:length(Gene_groups_2)){
    tmp_i <- Gene_groups_1[[i]]
    tmp_k <- Gene_groups_2[[k]]
    mat[k,i] <- phyper(length(tmp_i[tmp_i %in% tmp_k]), length(tmp_k), length(All_genes_clones[!All_genes_clones %in% tmp_k]), length(tmp_i), lower.tail = F)
  }
}
mat <- -log10(mat)

#plot as heatmap
library(fields)
library(scales)
library(gplots)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(mat),length=51))
heatmap.2(mat,Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=rownames(mat),labCol=colnames(mat) )

#Figure 4E
##Autocrine pathways Analysis
Tert <- readRDS("Tert_Subtypes.rds")
Idents(Tert) <- "Subtype" ##clusters

# Convert Seurat object to CellChat input
# normalized data
data.input <- GetAssayData(Tert, assay = "RNA", slot = "data")
meta <- data.frame(labels = Idents(Tert), row.names = names(Idents(Tert)))

##Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

#Set CellChat database for humans
CellChatDB <- CellChatDB.human  
cellchat@DB <- CellChatDB

#subset signaling genes
cellchat <- subsetData(cellchat)
# identify enriched ligands/receptors                       
cellchat <- identifyOverExpressedGenes(cellchat) 
#For each overexpressed ligand and receptor obtained inthe previous step, identify over-expressed L–R interactions if either its associated ligand or receptor is over expressed      
cellchat <- identifyOverExpressedInteractions(cellchat)
#Infer cell–cell communication at a L–R pair level
cellchat <- computeCommunProb(cellchat)
#Filter the cell–cell communication, based on the number of cells in each group
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Infer cell–cell communication at a signaling pathway level (CellChat computes the communication probability at the signaling pathway level by summarizing the communication probabilities of all L–R pairs associated with each signaling pathway.)
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell–cell communication network. CellChat calculates the aggregated cell–cell communication network by counting the number of links or summarizing the communication probability across all the cell groups (option A) or a subset of cell groups (option B).
cellchat <- aggregateNet(cellchat)
comm <- subsetCommunication(cellchat)
#extract autocrine communications
autocrine <- comm[comm$source == comm$target, ]
# Extract autocrine information for all pathways
autocrine_summary <- lapply(names(cellchat@netP$centr), function(pw) {
  cdata <- cellchat@netP$centr[[pw]]
  df <- data.frame(
    Pathway = pw,
    CellType = names(cdata$outdeg),
    AutocrineScore = (unlist(cdata$outdeg) + unlist(cdata$indeg)) / 2
  )
  return(df)
})

autocrine_summary <- do.call(rbind, autocrine_summary)
head(autocrine_summary)
autocrine_matrix <- reshape2::acast(autocrine_summary, Pathway ~ CellType, value.var = "AutocrineScore", fill = 0)

#scale by pathway or by cell type
autocrine_matrix_scaled <- t(scale(t(autocrine_matrix)))
pheatmap::pheatmap(autocrine_matrix_scaled,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   fontsize_row = 8, fontsize_col = 8,
                   main = "Autocrine signaling activity across pathways and cell types in each cluster")

top_pathways <- autocrine_summary %>%
  dplyr::mutate(AutocrineScore = as.numeric(AutocrineScore)) %>%  
  dplyr::group_by(Pathway) %>%
  dplyr::summarise(TotalScore = sum(AutocrineScore, na.rm = TRUE)) %>%
  dplyr::arrange(desc(TotalScore)) %>%
  dplyr::slice(1:20) %>%
  dplyr::pull(Pathway)

pheatmap::pheatmap(autocrine_matrix_scaled[top_pathways, ],
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 20 autocrine signaling pathways",
         cluster_rows = TRUE, cluster_cols = TRUE)

#Figure 4F
Idents(Tert) <- "Dataset" ##the cell lines

##reshape the matirx
reshape2::melt(cellchat@net$count)
library(qgraph)
#matrix tmp with column A and B indcating interaction from cluster/cell type A to cluster/cell type B with an additional column showing the strength as a numeric value
tmp<-reshape2::melt(cellchat@net$count)

tmp_df <- tmp  
colnames(tmp_df) <- c("from", "to", "weight")

# Convert edge list → adjacency matrix
adj_mat <- reshape2::acast(tmp_df, from ~ to, value.var = "weight", fill = 0)

# Convert to numeric matrix
adj_mat <- as.matrix(adj_mat)
mode(adj_mat) <- "numeric"

# Create a manual layout matrix
manual_layout <- matrix(NA, nrow = nrow(adj_mat), ncol = 2)
rownames(manual_layout) <- colnames(adj_mat)

# Set x–y coordinates manually
manual_layout["AD10", ] <- c(0.3, 0.6)   # AD10 (high bone-forming)
manual_layout["DD8", ] <- c(0.5, 0.5)   # close to AD10
manual_layout["CB",  ] <- c(0.1, 0.3)   # low bone-forming
manual_layout["CD",  ] <- c(0.5, 0.2)   # close to CB
manual_layout["Msc", ] <- c(0.7, 0.6)   # central hub

# Show spring layout for all clusters - nondeterministic and can result in different visualization
qgraph(adj_mat,
       layout = manual_layout,
       directed = TRUE,
       color = "white",
       posCol = "black",
       negCol = "red",
       edge.width = 2,
       vsize = 9.7,
       label.cex = 1.2)
legend(
  "topright",
  legend = c("Strong interaction", "Weak interaction"),
  lwd = c(4, 1),
  col = "black",
  bty = "n",
  title = "Interaction strength"
)	   
   
legend(
  "bottomright",
  legend = c("Signal direction"),
  lwd = 2,
  col = "black",
  pch = 62,  # arrow-like symbol
  bty = "n"
)
	   

#Figure 4G
Implant <- read.delim("Implants_sortet.txt",h=T)
Implant$Clones <- factor(Implant$Clones, levels=c('D100','D75C25','D50C50','D25C75','C100'))
Implant$Position <- factor(Implant$Position, levels=c('LF','RF','LR','RR'))

Implant_percentage <- data.frame(matrix(NA,ncol=ncol(Implant), nrow=1))
Implant_percentage <- Implant_percentage[-1,]

for (i in unique(Implant$Position)){
  tmp <- Implant[Implant$Position==i,]
  for (k in tmp$Experiment){
    tmp[tmp$Experiment==k,'Average_Bone'] <- tmp[tmp$Experiment==k,'Average_Bone']/tmp[tmp$Experiment==k & tmp$Clones=="D100",'Average_Bone']
  }
  Implant_percentage <- rbind(Implant_percentage,tmp)
}

Implant_percentage_avg <- data.frame(matrix(NA,ncol=3, nrow=5))
rownames(Implant_percentage_avg) <- c('D100','D75C25','D50C50','D25C75','C100')
colnames(Implant_percentage_avg) <- c('mean','STD','SE')

for (i in unique(Implant_percentage$Clones)){
  tmp <- Implant_percentage[Implant_percentage$Clones==i,]
  Implant_percentage_avg[i,1] <- mean(tmp$Average_Bone)
  Implant_percentage_avg[i,2] <- sd(tmp$Average_Bone)
  Implant_percentage_avg[i,3] <- sd(tmp$Average_Bone)/sqrt(nrow(tmp))
}

plot(1:5,Implant_percentage_avg$mean, type="l")
lines(1:5,Implant_percentage_avg$mean + Implant_percentage_avg$STD, lty=2)
lines(1:5,Implant_percentage_avg$mean - Implant_percentage_avg$STD, lty=2)
abline(1.25,-0.25, lty=3)

plot(1:5,Implant_percentage_avg$mean, type="l", xlab="", ylab="", xaxt="none")
axis(1,at=c(1:5),labels = c('D100','D75C25','D50C50','D25C75','C100'),las=2)
lines(1:5,Implant_percentage_avg$mean + Implant_percentage_avg$SE, lty=2)
lines(1:5,Implant_percentage_avg$mean - Implant_percentage_avg$SE, lty=2)
abline(1.25,-0.25, lty=3)


text(2,0.9,
  t.test(Implant_percentage[Implant_percentage$Clones=="D75C25",'Average_Bone'], mu=0.75)$p.value)
text(3,0.7,
     t.test(Implant_percentage[Implant_percentage$Clones=="D50C50",'Average_Bone'], mu=0.5)$p.value)
text(4,0.5,
     t.test(Implant_percentage[Implant_percentage$Clones=="D25C75",'Average_Bone'], mu=0.25)$p.value)
  
# Final aesthetics were done in Illustrator
