##Figure 3
#Figure 3A
#Stacked bar plots of fraction of cells per each cluster in each cell line
Idents(object = Tert) <- "RNA_snn_res.0.15"
ggplot(Tert@meta.data, aes(x=Dataset, fill=RNA_snn_res.0.15)) + geom_bar(position = "fill")+scale_fill_manual(values = col)

#Figure 3B
#Load the following objects from https://osf.io/wxpgn/
Cl_0 <- read.delim("Wiki_Cl_0.txt", h=T)
Cl_1 <- read.delim("Wiki_Cl_1.txt", h=T)
Cl_2 <- read.delim("Wiki_Cl_2.txt", h=T)
Cl_3 <- read.delim("Wiki_Cl_3.txt", h=T)
Cl_4 <- read.delim("Wiki_Cl_4.txt", h=T)
Cl_5 <- read.delim("Wiki_Cl_5.txt", h=T)
Cl_6 <- read.delim("Wiki_Cl_6.txt", h=T)
Cl_7 <- read.delim("Wiki_Cl_7.txt", h=T)
Cl_8 <- read.delim("Wiki_Cl_8.txt", h=T)
Cl_9 <- read.delim("Wiki_Cl_9.txt", h=T)
Cl_10 <- read.delim("Wiki_Cl_10.txt", h=T)


Cl_0 <- Cl_0[,c(2,10,11)]
names(Cl_0)[1] <- "Pathway"
names(Cl_0)[2] <- "pval_0"
names(Cl_0)[3] <- "padj_0"

Cl_0 <- subset(Cl_0, subset = Cl_0$pval_0 <0.01)
Cl_0 <- subset(Cl_0, subset = Cl_0$padj_0 <0.01)

Cl_1 <- Cl_1[,c(2,10,11)]
names(Cl_1)[1] <- "Pathway"
names(Cl_1)[2] <- "pval_1"
names(Cl_1)[3] <- "padj_1"

Cl_1 <- subset(Cl_1, subset = Cl_1$pval_1 <0.01)
Cl_1 <- subset(Cl_1, subset = Cl_1$padj_1 <0.01)

Cl_2 <- Cl_2[,c(2,10,11)]
names(Cl_2)[1] <- "Pathway"
names(Cl_2)[2] <- "pval_2"
names(Cl_2)[3] <- "padj_2"

Cl_2 <- subset(Cl_2, subset = Cl_2$pval_2 <0.01)
Cl_2 <- subset(Cl_2, subset = Cl_2$padj_2 <0.01)

Cl_3 <- Cl_3[,c(2,10,11)]
names(Cl_3)[1] <- "Pathway"
names(Cl_3)[2] <- "pval_3"
names(Cl_3)[3] <- "padj_3"

Cl_3 <- subset(Cl_3, subset = Cl_3$pval_3 <0.01)
Cl_3 <- subset(Cl_3, subset = Cl_3$padj_3 <0.01)

Cl_4 <- Cl_4[,c(2,10,11)]
names(Cl_4)[1] <- "Pathway"
names(Cl_4)[2] <- "pval_4"
names(Cl_4)[3] <- "padj_4"

Cl_4 <- subset(Cl_4, subset = Cl_4$pval_4 <0.01)
Cl_4 <- subset(Cl_4, subset = Cl_4$padj_4 <0.01)


Cl_5 <- Cl_5[,c(2,10,11)]
names(Cl_5)[1] <- "Pathway"
names(Cl_5)[2] <- "pval_5"
names(Cl_5)[3] <- "padj_5"

Cl_5 <- subset(Cl_5, subset = Cl_5$pval_5 <0.01)
Cl_5 <- subset(Cl_5, subset = Cl_5$padj_5 <0.01)

Cl_6 <- Cl_6[,c(2,10,11)]
names(Cl_6)[1] <- "Pathway"
names(Cl_6)[2] <- "pval_6"
names(Cl_6)[3] <- "padj_6"

Cl_6 <- subset(Cl_6, subset = Cl_6$pval_6 <0.01)
Cl_6 <- subset(Cl_6, subset = Cl_6$padj_6 <0.01)

Cl_7 <- Cl_7[,c(2,10,11)]
names(Cl_7)[1] <- "Pathway"
names(Cl_7)[2] <- "pval_7"
names(Cl_7)[3] <- "padj_7"

Cl_7 <- subset(Cl_7, subset = Cl_7$pval_7 <0.01)
Cl_7 <- subset(Cl_7, subset = Cl_7$padj_7 <0.01)

Cl_8 <- Cl_8[,c(2,10,11)]
names(Cl_8)[1] <- "Pathway"
names(Cl_8)[2] <- "pval_8"
names(Cl_8)[3] <- "padj_8"

Cl_8 <- subset(Cl_8, subset = Cl_8$pval_8 <0.01)
Cl_8 <- subset(Cl_8, subset = Cl_8$padj_8 <0.01)

Cl_9 <- Cl_9[,c(2,10,11)]
names(Cl_9)[1] <- "Pathway"
names(Cl_9)[2] <- "pval_9"
names(Cl_9)[3] <- "padj_9"

Cl_9 <- subset(Cl_9, subset = Cl_9$pval_9 <0.01)
Cl_9 <- subset(Cl_9, subset = Cl_9$padj_9 <0.01)

Cl_10 <- Cl_10[,c(2,10,11)]
names(Cl_10)[1] <- "Pathway"
names(Cl_10)[2] <- "pval_10"
names(Cl_10)[3] <- "padj_10"

Cl_10 <- subset(Cl_10, subset = Cl_10$pval_10 <0.01)
Cl_10 <- subset(Cl_10, subset = Cl_10$padj_10 <0.01)

ALL <- merge(Cl_0,Cl_1, all=T )
ALL <- merge(ALL,Cl_2, all=T )
ALL <- merge(ALL,Cl_3, all=T )
ALL <- merge(ALL,Cl_4, all=T )
ALL <- merge(ALL,Cl_5, all=T )
ALL <- merge(ALL,Cl_6, all=T )
ALL <- merge(ALL,Cl_7, all=T )
ALL <- merge(ALL,Cl_8, all=T )
ALL <- merge(ALL,Cl_9, all=T )
ALL <- merge(ALL,Cl_10, all=T )

ALL[is.na(ALL)] <- 0

WIKI <- ALL[,c(1,2,4,6,8,10,12,14,16,18,20,22)]
WIKI <- cbind(WIKI[,1],-log10(WIKI[,c(2:ncol(WIKI))]))
WIKI[WIKI=="Inf"]<-0

library(fields)
library(gplots)
WIKI <-unique(WIKI)
#Heatmap of selected WIKI pathways
WIKI_selected <-c("Glycolysis and Gluconeogenesis", "Cori Cycle", "VEGFA-VEGFR2 Signaling Pathway", "Focal Adhesion-PI3K-Akt-mTOR-signaling pathway", "TGF-beta Signaling Pathway", 
                     "Type I collagen synthesis in the context of Osteogenesis imperfecta", "PI3K/AKT/mTOR - VitD3 Signalling", "Pathways Regulating Hippo Signaling", "TGF-beta Receptor Signalling in Skeletal Dysplasias", 
					 "Prostaglandin Synthesis and Regulation", "Glutathione metabolism", "NRF2 pathway", "Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling", "EGF/EGFR Signaling Pathway", 
					 "Regulation of Actin Cytoskeleton", "Signaling of Hepatocyte Growth Factor Receptor", "MAPK Signaling Pathway", "Androgen receptor signaling pathway", "Integrin-mediated Cell Adhesion", 
					 "Regulation of Microtubule Cytoskeleton", "Insulin Signaling", "Ras Signaling", "Bone Morphogenic Protein (BMP) Signalling and Regulation", "Cell Cycle", "G1 to S cell cycle control", 
					 "Translation Factors", "Amino Acid metabolism", "Electron Transport Chain (OXPHOS system in mitochondria)", "Nucleotide Metabolism", "Proteasome Degradation", "Purine metabolism and related disorders", 
					 "Calcium Regulation in the Cardiac Cell", "G Protein Signaling Pathways", "Mitochondrial complex I assembly model OXPHOS system", "Mitochondrial CIV Assembly", "Mitochondrial CIII assembly", 
					 "Parkin-Ubiquitin Proteasomal System pathway", "Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds", "Senescence and Autophagy in Cancer")

p2 <-  WIKI %>%
  dplyr::filter(Pathway %in% WIKI_selected)

mat_vals <- as.matrix(p2[, 2:ncol(p2)])

rownames(mat_vals) <- p2[, 1]
mat_bin <- mat_vals
mat_bin <- ifelse(mat_vals >= 1, 1, 0)
rowSums(mat_bin)
row_dist <- dist(mat_bin, method = "binary")
row_hc   <- hclust(row_dist, method = "average")
row_dend <- as.dendrogram(row_hc)
plot(row_hc)
mat_col <- c("white",designer.colors(n = 50, col = c("plum1", "darkmagenta")))

n_col <- length(mat_col)

mat_col_breaks <- seq(
  from = min(mat_vals, na.rm = TRUE),
  to   = max(mat_vals, na.rm = TRUE),
  length.out = n_col + 1
)
mat_col_breaks <- c(0,seq(-log10(0.05),max(mat_vals[,1:ncol(mat_vals)]),length=51))

library(gplots)

heatmap.2(
  mat_vals,
  main = "WIKI",
  Rowv = as.dendrogram(row_hc),
  Colv = FALSE,
  dendrogram = "row",
  scale = "none",
  col = mat_col,
  breaks = mat_col_breaks,
  trace = "none",
  labRow = rownames(mat_vals),
  labCol = colnames(mat_vals),
  cexCol = 0.75
)

#Figure 3C
#Heatmap of selected Gene Ontology pathways
Tert.markers <- read.delim("Tert.markers.txt", h=T)

Gene_groups <- list()
Gene_groups[[1]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "0","gene"],c('gene')]
Gene_groups[[2]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "1","gene"],c('gene')]
Gene_groups[[3]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "2","gene"],c('gene')]
Gene_groups[[4]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "3","gene"],c('gene')]
Gene_groups[[5]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "4","gene"],c('gene')]
Gene_groups[[6]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "5","gene"],c('gene')]
Gene_groups[[7]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "6","gene"],c('gene')]
Gene_groups[[8]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "7","gene"],c('gene')]
Gene_groups[[9]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "8","gene"],c('gene')]
Gene_groups[[10]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "9","gene"],c('gene')]
Gene_groups[[11]] <- Tert.markers[Tert.markers$gene %in% Tert.markers[Tert.markers$cluster == "10","gene"],c('gene')]


names(Gene_groups) <- c('1','2','3','4','5','6','7','8','9','10','11')


tmp <- Gene_groups[[1]]
tmp_Cl <- make.names(tmp, unique = TRUE)

# Extracting all measured genes from the RNA-seq experiment by Symbol gene IDs 
tmp_All <- make.names(Tert.markers$gene, unique = TRUE)
# Converting DE genes to intergers, significant DE genes = 1, non-significant genes = 0
Cl <- as.integer(tmp_All %in% tmp_Cl)
names(Cl) <- tmp_All

# Fitting the probability weighting (PWF) function while correcting for gene count length
# I have used gene symbols as gene IDs, only symbol, Ensembl, or Entrez IDs are supported
pwf_Cl <- nullp(Cl, "hg19", "geneSymbol")
# GO analysis
# Here, the test is limited to biological processes (GO:BP), change to GO:CC or GO:MF for cellular components or molecular function testing
GO.BP_Cl <- goseq(pwf_Cl, "hg19", "geneSymbol", test.cats = c("GO:BP"), use_genes_without_cat = TRUE)

# Applying multiple testing by Benjamini Hochberg to overenriched p values
GO.BP_Cl$over_represented_p.adjust <- p.adjust(GO.BP_Cl$over_represented_pvalue, method = "BH")
#Collecting GO terms
GO_cluster <- GO.BP_Cl[,c(1,6,8)]
names(GO_cluster)[3] <- names(Gene_groups)[1]

# Loop for the other clusters

for(i in 2:length(Gene_groups)){
  tmp <- Gene_groups[[i]]
  tmp_Cl <- make.names(tmp, unique = TRUE)
  if(length(tmp_Cl)>0){
    Cl <- as.integer(tmp_All %in% tmp_Cl)
    names(Cl) <- tmp_All
    pwf_Cl <- nullp(Cl, "hg19", "geneSymbol")
    GO.BP_Cl <- goseq(pwf_Cl, "hg19", "geneSymbol", test.cats = c("GO:BP"), use_genes_without_cat = TRUE)
    GO.BP_Cl$over_represented_p.adjust <- p.adjust(GO.BP_Cl$over_represented_pvalue, method = "BH")
    names(GO.BP_Cl)[8] <- names(Gene_groups)[i]
    GO_cluster <- merge(GO_cluster,GO.BP_Cl[,c(1,8)],by="category")
  } else {}
}

#Plotting
p <- cbind(GO_cluster[,1:2],-log10(GO_cluster[,c(3:ncol(GO_cluster))]))
p2 <- p[order(-p[,3]),][1:20,]


for(i in 4:ncol(GO_cluster)){
  p3 <- p[order(-p[,i]),][1:20,]
  p2 <- rbind(p2,p3)
}
p2[p2=="Inf"]<-10
p2 <- unique(p2)


GOterms_selected<-c("cellular response to hypoxia", "regulation of cell migration", "hexose metabolic process", "lipid storage", "response to growth factor", "glucose metabolic process", "angiogenesis", 
                      "bone development", "extracellular matrix organization", "cell adhesion", "locomotion", "system development", "blood circulation", "cell communication", "G protein-coupled receptor signaling pathway", 
					  "lipid metabolic process", "fat cell differentiation", "mitotic cell cycle", "DNA repair", "chromatin organization", "NADH regeneration", "canonical glycolysis", "pyruvate metabolic process", 
					  "cell-cell adhesion", "proton motive force-driven mitochondrial ATP synthesis", "cytoplasmic translation")

p2 <-  GO_cluster %>%
  dplyr::filter(term %in% GOterms_selected)

mat_vals <- as.matrix(p2[, 3:ncol(p2)])
mat_vals <- -log10(mat_vals)
rownames(mat_vals) <- p2[, 2]
mat_bin <- mat_vals
mat_bin <- ifelse(mat_vals >= 1, 1, 0)
rowSums(mat_bin)
row_dist <- dist(mat_bin, method = "binary")
row_hc   <- hclust(row_dist, method = "average")
row_dend <- as.dendrogram(row_hc)
plot(row_hc)
mat_col <- c("white",designer.colors(n = 50, col = c("plum1", "darkmagenta")))

n_col <- length(mat_col)

mat_col_breaks <- seq(
  from = min(mat_vals, na.rm = TRUE),
  to   = max(mat_vals, na.rm = TRUE),
  length.out = n_col + 1
)


library(gplots)

heatmap.2(
  mat_vals,
  main = "GO cluster",
  Rowv = as.dendrogram(row_hc),
  Colv = FALSE,
  dendrogram = "row",
  scale = "none",
  col = mat_col,
  breaks = mat_col_breaks,
  trace = "none",
  labRow = rownames(mat_vals),
  labCol = colnames(mat_vals),
  cexCol = 0.75
)


#Figure 3D
#Heatmap of REACTOME analysis

# Import libraries

library(org.Hs.eg.db)
library(clusterProfiler)
library(goseq)
library(reactome.db)

HGNC_matrix <- read.delim("HGNC_matrix.txt",h=T)
colnames(HGNC_matrix) <- c("HGNC","Symbol","ENTREZID","RefSeqID")



## Pathway analysis
# Setup the reactome database
Reactome <- as.data.frame(reactomeEXTID2PATHID)
Relation <- read.delim("ReactomePathwaysRelation.txt", header=FALSE)
Pathways <- read.delim("ReactomePathways.txt", header=FALSE)
Pathways <- Pathways[ Pathways$V3 == "Homo sapiens",]
Pathways$DB_ID <- substr(Pathways$V1, 7, nchar(as.character(Pathways$V1)))

Metabolism <- Relation[ Relation[,1] %in% as.character(Pathways[ Pathways$V2 == "Metabolism",1]),]  
Metabolism <- rbind(Metabolism, Relation[ Relation[,1] %in% Metabolism[,2],])
Metabolism <- Metabolism[ duplicated(Metabolism[,2])==FALSE,]

Current <- nrow(Metabolism)
New <- Current + 1
while (New > Current) {
  Current <- nrow(Metabolism)
  Metabolism <- rbind(Metabolism, Relation[ Relation[,1] %in% Metabolism[,2],])
  Metabolism <- Metabolism[ duplicated(Metabolism[,2])==FALSE,]
  New <- nrow(Metabolism)
}

Pathways <- Pathways[ Pathways$V1 %in% Metabolism[,1] | Pathways$V1 %in% Metabolism[,2],]
Pathways <- Pathways[ duplicated(Pathways$V1)==FALSE,]
Reactome <- Reactome[ Reactome$DB_ID %in% Pathways$V1,]
Reactome <- split(Reactome$gene_id, f = Reactome$DB_ID, drop=T)


#Convert <- Convert[ duplicated(Convert$SYMBOL_Mouse)==FALSE,]

Convert <- bitr(Tert.markers$gene,fromType = "SYMBOL",toType   = "ENTREZID",OrgDb    = org.Hs.eg.db)

# Remove duplicates
Convert <- Convert[!duplicated(Convert$ENTREZID), ]

# Merge with expression data
All <- merge(Convert,Tert.markers,by.x = "SYMBOL",by.y = "gene")

All <- unique(All)

Enrichment <- Pathways[, c(1,2)]
colnames(Enrichment) <- c("category","name")

# Ensure unique ENTREZ IDs
All <- All[!duplicated(All$ENTREZID), ]
All$ENTREZID <- as.character(All$ENTREZID)

for (i in 1:length(Gene_groups)) {
  
  # 0/1 vector for the gene group
  Cl1 <- as.integer(All$ENTREZID %in% All$ENTREZID[All$SYMBOL %in% Gene_groups[[i]]])
  names(Cl1) <- All$ENTREZID
  
  # Dummy bias vector for old goseq version
  bias.vec <- seq_along(Cl1)
  
  # Build null model
  Cl1.nullp <- nullp(
    DEgenes = Cl1,
    genome = NULL,
    id = NULL,
    bias.data = bias.vec
  )
  
  # Run Reactome enrichment
  Cl1.EA <- goseq(
    pwf = Cl1.nullp,
    gene2cat = Reactome,
    use_genes_without_cat = TRUE
  )
  
  # Adjust p-values
  Cl1.EA$Cl1 <- p.adjust(Cl1.EA$over_represented_pvalue, method = "fdr")
  colnames(Cl1.EA)[6] <- names(Gene_groups)[i]
  
  # Merge with results
  Enrichment <- merge(Enrichment, Cl1.EA[, c(1,6)], by="category", all=TRUE)
}

# Set Q-value NAs to 1 
for (i in 3:ncol(Enrichment)) {Enrichment[is.na(Enrichment[, i]), i] <- 1}

# Calculate enrichment in each pathway
Result <- data.frame(matrix(ncol=length(Gene_groups)+1, nrow=299))
for (q in 1:length(Reactome)) {
  Result[q,1] <- names(Reactome[q])
  for (i in 1:length(Gene_groups)){
    A <- nrow(All[ All$SYMBOL %in% Gene_groups[[i]] & All$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
    C <- nrow(All[ All$SYMBOL %in% Gene_groups[[i]],])
    E <- nrow(All[ All$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
    G <- nrow(All)
    # Add 1/1000 of the group size as pseudocount to avoid 0's
    A <- A + (C/1000)
    C <- C + (C/1000)
    Result[q,1+i] <- (A/E)/(C/G)
  }
}

colnames(Result) <- c("category",names(Gene_groups))
Enrichment <- merge(Enrichment, Result, by="category")
rownames(Enrichment) <- Enrichment$name
p <- Enrichment[apply(Enrichment[,3:(2+length(Gene_groups))],1,min)<0.05,]
for(i in 1:length(Gene_groups)){
  p[p[,2+i] > 0.05,(2+length(Gene_groups)+i)]<- 0
}

p2 <- cbind(p[,1:2],-log10(p[,c(3:(2+length(Gene_groups)))]))
p2 <- p2[order(-p2[,3]),]
p2 <- p2[p2[,3]> -log10(0.05),]
for (i in 4:(2+length(Gene_groups))){
  p3 <- cbind(p[,1:2],-log10(p[,c(3:(2+length(Gene_groups)))]))
  p3 <- p3[order(-p3[,i]),]
  p3 <- p3[p3[,i]> -log10(0.05),]
  p2 <- rbind(p2,p3)
}
p2 <-unique(p2)

selected_terms<-c("Fatty acid metabolism","Free fatty acids regulate insulin secretion","Biological oxidations","Glycolysis","Gluconeogenesis","Glucose metabolismt",
                   "Respiratory electron transport","The citric acid (TCA) cycle and respiratory electron transport")

p2 <-  p2 %>%
  dplyr::filter(name %in% selected_terms)

mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(p2[,3:(2+length(Gene_groups))]),length=51))
heatmap.2(as.matrix(p2[,3:(2+length(Gene_groups))]),main="REACTOME", Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none' )

#Figure 3E and 3G were ploted in GraphPad

#Figure 3F
#Load the following objects 

library("xlsx")
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

##basal respiration
BR <- read.xlsx("Basalrespiration.xlsx", 
                sheetIndex = 1, header=TRUE)


df_long <- BR %>%
  pivot_longer(cols = everything(),
               names_to = "Group",
               values_to = "BasalRespiration") %>%
  na.omit()

t.test(BasalRespiration ~ Group, data = df_long)

ggplot(df_long, aes(x = Group, y = BasalRespiration)) +
  geom_jitter(width = 0.1, size = 2) +
  
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5
  ) +
  
  stat_compare_means(method = "t.test") +
  
  ylim(1, 4) +
  theme_classic()

##maximal respiration
MR <- read.xlsx("MaximalRespiration.xlsx", 
                sheetIndex = 1, header=TRUE)

df_long <- MR %>%
  pivot_longer(cols = everything(),
               names_to = "Group",
               values_to = "MaximalRespiration") %>%
  na.omit()

t.test(MaximalRespiration ~ Group, data = df_long)

ggplot(df_long, aes(x = Group, y = MaximalRespiration)) +
  geom_jitter(width = 0.1, size = 2) +
  
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5
  ) +
  
  stat_compare_means(method = "t.test") +
  theme_classic()

#Figure 3H
##basal glycolysis
BG <- read.xlsx("Basalglycolisis.xlsx", 
                sheetIndex = 1, header=TRUE)

df_long <- BG %>%
  pivot_longer(cols = everything(),
               names_to = "Group",
               values_to = "BasalGlycolysis") %>%
  na.omit()

t.test(BasalGlycolysis ~ Group, data = df_long)


ggplot(df_long, aes(x = Group, y = BasalGlycolysis)) +
  geom_jitter(width = 0.1, size = 2) +
  
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5
  ) +
  
  stat_compare_means(method = "t.test") +
  
  theme_classic()

##Compensatory glycolysis
CG <- read.xlsx("Compensatoryglycolysis.xlsx", 
                sheetIndex = 1, header=TRUE)
library(tidyr)
library(dplyr)

df_long <- CG %>%
  pivot_longer(cols = everything(),
               names_to = "Group",
               values_to = "CompensatoryGlycolysis") %>%
  na.omit()

t.test(CompesatoryGlycolysis ~ Group, data = df_long)

ggplot(df_long, aes(x = Group, y = CompesatoryGlycolysis)) +
  geom_jitter(width = 0.1, size = 2) +
  
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5
  ) +
  
  stat_compare_means(method = "t.test") +
  theme_classic()

#Figure 3I
##Eenergy map
df <- read.xlsx("Energymap.xlsx", 
                sheetIndex = 1, header=TRUE)



summary_df <- data.frame(
  Group = c("DD8","CB4"),
  
  OCR_mean = c(mean(df$DD8_OCR, na.rm=TRUE),
               mean(df$CB4_OCR, na.rm=TRUE)),
  
  OCR_sd = c(sd(df$DD8_OCR, na.rm=TRUE),
             sd(df$CB4_OCR, na.rm=TRUE)),
  
  ECAR_mean = c(mean(df$DD8_ECAR, na.rm=TRUE),
                mean(df$CB4_ECAR, na.rm=TRUE)),
  
  ECAR_sd = c(sd(df$DD8_ECAR, na.rm=TRUE),
              sd(df$CB4_ECAR, na.rm=TRUE))
)


ggplot(summary_df, aes(x = ECAR_mean, y = OCR_mean, label = Group)) +
  
  geom_point(size = 4) +
  
  geom_errorbarh(aes(xmin = ECAR_mean - ECAR_sd,
                     xmax = ECAR_mean + ECAR_sd),
                 height = 0.1) +
  
  geom_errorbar(aes(ymin = OCR_mean - OCR_sd,
                    ymax = OCR_mean + OCR_sd),
                width = 0.1) +
  
  geom_text(vjust = -1) +
  
  coord_cartesian(xlim = c(0.5, 5), ylim = c(0, 10)) +
  
  theme_classic() +
  xlab("ECAR") +
  ylab("OCR")

# Final aesthetics were done in Illustrator
