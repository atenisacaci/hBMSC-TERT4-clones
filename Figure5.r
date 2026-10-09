###Figure 5
##Figure 5A
Tert <- readRDS("Tert_Subtypes.rds")
Idents(Tert) <- "Dataset" ##the cell lines
# Convert Seurat object to CellChat input
# normalized data
data.input <- GetAssayData(Tert, assay = "RNA", slot = "data")
meta <- data.frame(labels = Idents(Tert), row.names = names(Idents(Tert)))

##Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

#Set CellChat database for humans
CellChatDB <- CellChatDB.human  
cellchat@DB <- CellChatDB

##get all the ligands and receptors from the CellChat 
CellChatDB.human$interaction

db_full <- CellChatDB.human$interaction

head(db_full[, c("interaction_name", "pathway_name", "ligand", "receptor")])

##extract signaling genes from CellChat 
extract_cellchatDB_genes <- function(db) {
  
  cols <- c("ligand", "receptor", "agonist", "antagonist",
            "co_A_receptor", "co_I_receptor")
  
  genes <- unlist(db[, cols], use.names = FALSE)
  genes <- genes[!is.na(genes)]
  genes <- genes[genes != ""]
  
  # split complexes like ITGA9_ITGB1
  genes <- unlist(strsplit(genes, "_"))
  
  genes <- trimws(genes)
  genes <- unique(genes)
  genes <- sort(genes)
  
  return(genes)
}

cellchatDB_genes <- extract_cellchatDB_genes(CellChatDB.human$interaction)

length(cellchatDB_genes)
head(cellchatDB_genes, 100)

Tert<-subset(Tert, idents = c("AD10","DD8","CB","CD"))
HIPPO_genes <- c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")
Parkin_genes <-c("TUBA1A","PSMD8","PSMC5","PSMC3","HSPA8","TUBB","TUBB2B","TUBB2A","PSMC4","UBE2L6","TUBB4B","PSMD13","PSMD4","PSMC2","PSMD2","PSMD6","PSMD3","UBE2L3","HSPA5","PSMD1","TUBB6","PSMD14")
Ca_genes <-c("GNG12","GNAI2","ADCY9","YWHAB","GNG11","PRKCH","CALM1","GNAS","CAMK2D","RGS20","ADCY4","PRKAR1A","CACNA1C","ITPR1","YWHAQ","GNAQ","YWHAE","PRKD1")
TGF_genes <-c("STAT1","SMAD6","SMAD9","ZEB2","FBN1","INHBA","ENG","HRAS","SMAD3")
Gprotein_genes <-c("RRAS","GNG12","GNAI2","ADCY9","AKAP12","GNG11","PRKCH","CALM1","GNAS","PDE7B","ADCY4","PRKAR1A","PDE4A","ITPR1","HRAS","GNAQ","PRKD1","PDE1C")

##compare the overlap with our signature 
signatures <- list(
  HIPPO = HIPPO_genes,
  Parkin = Parkin_genes,
  Ca = Ca_genes,
  TGF = TGF_genes,
  Gprotein = Gprotein_genes
)


#extract all ligands and receptors from the full CellChat database
extract_cellchat_lr <- function(db) {

  ligands <- db$ligand
  receptors <- db$receptor

  ligands <- ligands[!is.na(ligands) & ligands != ""]
  receptors <- receptors[!is.na(receptors) & receptors != ""]

  # split complexes like ITGA9_ITGB1
  ligands <- unique(trimws(unlist(strsplit(ligands, "_"))))
  receptors <- unique(trimws(unlist(strsplit(receptors, "_"))))

  list(
    ligands = sort(ligands),
    receptors = sort(receptors)
  )
}
db <- CellChatDB.human$interaction
lr <- extract_cellchat_lr(db)

cellchat_ligands <- lr$ligands
cellchat_receptors <- lr$receptors
db <- CellChatDB.human$interaction
lr <- extract_cellchat_lr(db)

cellchat_ligands <- lr$ligands
cellchat_receptors <- lr$receptors
signatures <- list(
  HIPPO = HIPPO_genes,
  Parkin = Parkin_genes,
  Ca = Ca_genes,
  TGF = TGF_genes,
  Gprotein = Gprotein_genes
)
lr_overlap <- do.call(rbind, lapply(names(signatures), function(sig) {

  sig_genes <- unique(signatures[[sig]])

  ligand_overlap <- intersect(sig_genes, cellchat_ligands)
  receptor_overlap <- intersect(sig_genes, cellchat_receptors)

  data.frame(
    Signature = sig,
    Signature_Size = length(sig_genes),

    Ligand_N = length(ligand_overlap),
    Ligand_Genes = paste(ligand_overlap, collapse = ", "),

    Receptor_N = length(receptor_overlap),
    Receptor_Genes = paste(receptor_overlap, collapse = ", "),

    Total_LR_N = length(unique(c(ligand_overlap, receptor_overlap))),
    Total_LR_Genes = paste(unique(c(ligand_overlap, receptor_overlap)), collapse = ", ")
  )
}))
lr_overlap[order(-lr_overlap$Total_LR_N), ]

plot_df <- data.frame(
  Signature = lr_overlap$Signature,

  Ligands = lr_overlap$Ligand_N,

  Receptors = lr_overlap$Receptor_N,

  Other = lr_overlap$Signature_Size -
           lr_overlap$Total_LR_N
)


plot_long <- pivot_longer(
  plot_df,
  cols = c(Ligands, Receptors, Other),
  names_to = "Category",
  values_to = "Count"
)

##plot ligand/recptor composition of pathway signature
ggplot(plot_long,
       aes(x = Signature,
           y = Count,
           fill = Category)) +

  geom_bar(stat = "identity") +

  theme_classic(base_size = 14) +

  scale_fill_manual(
    values = c(
      Ligands = "tomato",
      Receptors = "steelblue",
      Other = "grey80"
    )
  ) +

  labs(
    title = "Ligand/Receptor composition of pathway signatures",
    y = "Number of genes",
    x = ""
  )

##Figure 5B
Idents(Tert) <- "Dataset"
Tert<-subset(Tert, idents = c("AD10","DD8","CB","CD"))
HIPPO_genes <- c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")
Parkin_genes <-c("TUBA1A","PSMD8","PSMC5","PSMC3","HSPA8","TUBB","TUBB2B","TUBB2A","PSMC4","UBE2L6","TUBB4B","PSMD13","PSMD4","PSMC2","PSMD2","PSMD6","PSMD3","UBE2L3","HSPA5","PSMD1","TUBB6","PSMD14")
Ca_genes <-c("GNG12","GNAI2","ADCY9","YWHAB","GNG11","PRKCH","CALM1","GNAS","CAMK2D","RGS20","ADCY4","PRKAR1A","CACNA1C","ITPR1","YWHAQ","GNAQ","YWHAE","PRKD1")
TGF_genes <-c("STAT1","SMAD6","SMAD9","ZEB2","FBN1","INHBA","ENG","HRAS","SMAD3")
Gprotein_genes <-c("RRAS","GNG12","GNAI2","ADCY9","AKAP12","GNG11","PRKCH","CALM1","GNAS","PDE7B","ADCY4","PRKAR1A","PDE4A","ITPR1","HRAS","GNAQ","PRKD1","PDE1C")

pathways <- list(
  HIPPO = HIPPO_genes,
  Parkin=Parkin_genes,
  Ca = Ca_genes,
  TGF = TGF_genes,
  Gprotein = Gprotein_genes
)
for (p in names(pathways)) {
  Tert <- AddModuleScore(
    Tert,
    features = list(pathways[[p]]),
    name = p
  )
}


##ROC curve 
roc_df <- Tert@meta.data[, c(
  "Dataset",
  "HIPPO1",
  "Parkin1",
  "Ca1",
  "TGF1",
  "Gprotein1"
)]

roc_df$Group <- ifelse(roc_df$Dataset %in% c("AD10", "DD8"), 1, 0)
roc_df$Group <- as.factor(roc_df$Group)
head(roc_df)
table(roc_df$Dataset, roc_df$Group)

library(pROC)

roc_hippo <- roc(roc_df$Group, roc_df$HIPPO1)

plot(roc_hippo, print.auc = TRUE, main = "ROC - HIPPO")


pathway_cols <- c("HIPPO1", "Parkin1", "Ca1", "TGF1", "Gprotein1")

roc_list <- lapply(pathway_cols, function(p) {
  roc_obj <- roc(roc_df$Group, roc_df[[p]])
  data.frame(
    Pathway = p,
    AUC = as.numeric(auc(roc_obj))
  )
})

roc_results <- do.call(rbind, roc_list)
roc_results <- roc_results[order(-roc_results$AUC), ]
roc_results

roc_colors <- c("red","blue","darkgreen","purple","orange","brown","black")

plot(roc(roc_df$Group, roc_df[[pathway_cols[1]]]),
     col = roc_colors[1],
     lwd = 2,
     main = "ROC Curves for Pathway Scores")

for (i in 2:length(pathway_cols)) {
  plot(roc(roc_df$Group, roc_df[[pathway_cols[i]]]),
       col = roc_colors[i],
       lwd = 2,
       add = TRUE)
}

legend("bottomright",
       legend = paste0(roc_results$Pathway, " (AUC=", round(roc_results$AUC, 2), ")"),
       col = roc_colors[match(roc_results$Pathway, pathway_cols)],
       lwd = 2,
       cex = 0.8)


roc_results$Pathway <- factor(roc_results$Pathway, levels = roc_results$Pathway)

ggplot(roc_results, aes(x = Pathway, y = AUC, fill = Pathway)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_classic() +
  ylim(0, 1) +
  labs(
    title = "Pathway classification performance",
    y = "AUC",
    x = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
##Figure 5C
#load the following object
library(R.utils)
GSE253355 <- readRDS("GSE253355_MSC_Subset_Seurat.rds")
Idents(object = GSE253355) <- "original_seurat_clusters"
Clusters  <- sort(levels(GSE253355$original_seurat_clusters))
Result <- list()

# Defining enriched genes, e.g. one cluster versus all other data points
GSE253355_markers <- FindAllMarkers(GSE253355, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
GSE253355_markers <- GSE253355_markers[GSE253355_markers$p_val_adj < 0.05,]

# Defining marker genes, e.g. one cluster versus all other clusters individually and clean up results in a list
cluster_name <- Clusters
for(clust in Clusters){
	Pairs <- data.frame(GroupA = clust, GroupB = Clusters[!(Clusters %in% clust)])
	Exclusive <- c()
	for (i in 1:nrow(Pairs)){
		tmp <- FindMarkers(GSE253355, ident.1 = clust, ident.2 = Pairs[i,2], min.pct = 0.1, only.pos = TRUE)
		tmp <- tmp[tmp$p_val_adj < 0.05,]
		Exclusive <- c(Exclusive,as.character(rownames(tmp)))
		print(paste("Cluster",clust,"out of",length(Clusters),"versus Cluster",Pairs[i,2],sep=" "))
	}
	tmp <- data.frame(table(Exclusive))
	tmp1 <- tmp[tmp$Freq==length(Clusters)-1,]
	tmp2 <- tmp[tmp$Freq > 0,]
	if(length(Exclusive) >= 0 & nrow(tmp1)>0 & length(GSE253355_markers[GSE253355_markers$cluster == clust,'gene'])>0){
			Result[[clust]] <- rbind(
			data.frame("Gene"= GSE253355_markers[GSE253355_markers$cluster == clust,'gene'],"Marker"="Enriched"),
			data.frame("Gene"= tmp1[,'Exclusive'],"Marker"="Exclusive"),
			data.frame("Gene"= tmp2[,'Exclusive'],"Marker"="Significant"))
	} else if(length(Exclusive) > 0 & length(GSE253355_markers[GSE253355_markers$cluster == clust,'gene'])>0){
			Result[[clust]] <- rbind(
			data.frame("Gene"= GSE253355_markers[GSE253355_markers$cluster == clust,'gene'],"Marker"="Enriched"),
			data.frame("Gene"= tmp2[,'Exclusive'],"Marker"="Significant"))
	} else if(length(GSE253355_markers[GSE253355_markers$cluster == clust,'gene'])>=0){
			Result[[clust]] <- rbind(
			data.frame("Gene"= GSE253355_markers[GSE253355_markers$cluster == clust,'gene'],"Marker"="Enriched"))
	} else if(length(Exclusive) > 0){
			Result[[clust]] <- rbind(
			data.frame("Gene"= tmp2[,'Exclusive'],"Marker"="Significant"))
	} else {
			cluster_name <- cluster_name[!cluster_name %in% clust]
	}

}
names(Result) <- paste("Cl_",cluster_name,sep="")
Cluster_markers <- bind_rows(Result, .id = "column_label")
# Stats on enriched and exclusive Markers
mat <- matrix(NA,ncol=3,nrow=length(Clusters))
rownames(mat) <- Clusters
colnames(mat) <- c('Exclusive','Enriched','Significant')
for (i in 1:length(Clusters)){
	tmp <- Result[[i]]
	mat[i,1] <- nrow(tmp[tmp$Marker =="Exclusive",])
	mat[i,2] <- nrow(tmp[tmp$Marker =="Enriched",])
	mat[i,3] <- nrow(tmp[tmp$Marker =="Significant",])
}
mat

###load the following object
Clones <- read.delim("Cluster_markers_GSE253355.txt",h=T)
All_genes_clones <- readRDS("All_genes_clones.rds")

# Combine markers of single RNA-seq in one data frame and split enriched markers in a list
markers<-Clones
All_genes_clones <- All_genes_clones[All_genes_clones %in% Clones$Gene]
Clones <- Clones[Clones$Gene %in% All_genes_clones,]
markers <- markers[markers$Gene %in% All_genes_clones,]

Gene_groups <- list()
Gene_groups[[1]] <- markers[markers$column_label== "Cl_Adipo-MSC" & markers$Marker =="Exclusive",'Gene' ]
Gene_groups[[2]] <- markers[markers$column_label== "Cl_Fibro-MSC" & markers$Marker =="Exclusive",'Gene' ]
Gene_groups[[3]] <- markers[markers$column_label== "Cl_Osteo-MSC" & markers$Marker =="Exclusive",'Gene' ]
Gene_groups[[4]] <- markers[markers$column_label== "Cl_Osteoblast" & markers$Marker =="Exclusive",'Gene' ]
Gene_groups[[5]] <- markers[markers$column_label== "Cl_OsteoFibro-MSC" & markers$Marker =="Exclusive",'Gene' ]
Gene_groups[[6]] <- markers[markers$column_label== "Cl_THY1+ MSC" & markers$Marker =="Exclusive",'Gene' ]


names(Gene_groups) <- c('Cl_Adipo-MSC','Cl_Fibro-MSC','Cl_Osteo-MSC','Cl_Osteoblast','Cl_OsteoFibro-MSC','Cl_THY1+ MSC')

HIPPO_genes <- c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")
Parkin_genes <-c("TUBA1A","PSMD8","PSMC5","PSMC3","HSPA8","TUBB","TUBB2B","TUBB2A","PSMC4","UBE2L6","TUBB4B","PSMD13","PSMD4","PSMC2","PSMD2","PSMD6","PSMD3","UBE2L3","HSPA5","PSMD1","TUBB6","PSMD14")
Ca_genes <-c("GNG12","GNAI2","ADCY9","YWHAB","GNG11","PRKCH","CALM1","GNAS","CAMK2D","RGS20","ADCY4","PRKAR1A","CACNA1C","ITPR1","YWHAQ","GNAQ","YWHAE","PRKD1")
TGF_genes <-c("STAT1","SMAD6","SMAD9","ZEB2","FBN1","INHBA","ENG","HRAS","SMAD3")
Gprotein_genes <-c("RRAS","GNG12","GNAI2","ADCY9","AKAP12","GNG11","PRKCH","CALM1","GNAS","PDE7B","ADCY4","PRKAR1A","PDE4A","ITPR1","HRAS","GNAQ","PRKD1","PDE1C")

Gene_signatures <- list(
  HIPPO = HIPPO_genes,
  Parkin = Parkin_genes,
  Calcium = Ca_genes,
  TGF = TGF_genes,
  Gprotein = Gprotein_genes
)

Gene_signatures <- lapply(
  Gene_signatures,
  function(x) unique(x[x %in% All_genes_clones])
)
# Test the overlap of both gene groups using a hypergeometric test
mat <- matrix(NA, ncol=length(Gene_signatures),nrow=length(Gene_groups))
colnames(mat) <- names(Gene_signatures)
rownames(mat) <- names(Gene_groups)

for (i in 1:length(Gene_groups)){
  for (k in 1:length(Gene_signatures)){
    tmp_i <- Gene_groups[[i]]
    tmp_k <- Gene_signatures[[k]]
    x <- length(tmp_k[tmp_k %in% tmp_i])
    m <- length(tmp_k)
    n <- length(All_genes_clones[!All_genes_clones %in% tmp_k])
    y <- length(tmp_i)
    mat[i,k] <- phyper(x,m,n,y,lower.tail = F)
  }
}
# transform to log scale
mat <- -log10(mat)
mat[mat=="Inf"]<-187
# Show enrichment in a heatmap
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(6),length=51))
heatmap.2(mat,Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none')

##Figure 5D
#Load the following object
genes_df <- read.csv("ebmd_con_indepen.csv")
target_genes<- unique(genes_df$C.GENE)
HIPPO_genes <- c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")
Parkin_genes <-c("TUBA1A","PSMD8","PSMC5","PSMC3","HSPA8","TUBB","TUBB2B","TUBB2A","PSMC4","UBE2L6","TUBB4B","PSMD13","PSMD4","PSMC2","PSMD2","PSMD6","PSMD3","UBE2L3","HSPA5","PSMD1","TUBB6","PSMD14")
Ca_genes <-c("GNG12","GNAI2","ADCY9","YWHAB","GNG11","PRKCH","CALM1","GNAS","CAMK2D","RGS20","ADCY4","PRKAR1A","CACNA1C","ITPR1","YWHAQ","GNAQ","YWHAE","PRKD1")
TGF_genes <-c("STAT1","SMAD6","SMAD9","ZEB2","FBN1","INHBA","ENG","HRAS","SMAD3")
Gprotein_genes <-c("RRAS","GNG12","GNAI2","ADCY9","AKAP12","GNG11","PRKCH","CALM1","GNAS","PDE7B","ADCY4","PRKAR1A","PDE4A","ITPR1","HRAS","GNAQ","PRKD1","PDE1C")
signatures <- list(
  HIPPO = HIPPO_genes,
  Parkin = Parkin_genes,
  Ca = Ca_genes,
  TGF = TGF_genes,
  Gprotein = Gprotein_genes
)
overlap_df <- do.call(rbind, lapply(names(signatures), function(sig){
  
  sig_genes <- unique(signatures[[sig]])
  
  overlap <- intersect(sig_genes, target_genes)
  
  data.frame(
    Signature = sig,
    Signature_Size = length(sig_genes),
    Overlap_N = length(overlap),
    Overlap_Fraction = length(overlap) / length(sig_genes),
    Overlap_Genes = paste(overlap, collapse = ", "),
    pvalue = phyper(length(overlap),length(sig_genes),length(Atenisa_all)-length(sig_genes), length(target_genes), lower.tail = F)
  )
}))

#plot putative casual eBMD genes overlap with rwiki signaling pathways 
ggplot(overlap_df, aes(x = Signature, y = Overlap_N, fill = Overlap_N)) + 
  geom_bar(stat = "identity", color = "black", width = 0.7, linewidth = 0.3) +
  geom_text(aes(label = Overlap_N), vjust = -0.5, size = 5) + 
  scale_fill_gradientn(colours = c("white", "white", designer.colors(n = 50, col = c("plum1","darkmagenta"))), 
                       values = scales::rescale(c(0, 2, seq(2, 8, length.out = 50)), from = c(0, 8)), limits = c(0, 8), 
                       breaks = 0:8, oob = scales::squish, name = "Overlap\ncount") + 
  expand_limits(y = max(overlap_df$Overlap_N) + 2) + theme_classic(base_size = 14) + 
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11)) + 
  labs(title = "Gene overlap with signaling signatures", y = "Number of overlapping genes", x = NULL)


##Figure 5E
#Load the following objects from https://osf.io/wxpgn/
##Dotplot expession of Hippo compentency genes 
Tert <- readRDS("Tert_Subtypes.rds")
Idents(Tert) <- "Dataset"
DotPlot(Tert, features = c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")) + RotatedAxis()

# Final aesthetics were done in Illustrator
