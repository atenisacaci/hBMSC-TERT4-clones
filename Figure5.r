##Figure 5
#Figure 5A
#Load the following objects from https://osf.io/wxpgn/
##Dotplot expession of Hippo compentency genes 
Tert <- readRDS("Tert_Subtypes.rds")
Idents(Tert) <- "Dataset"
DotPlot(Tert, features = c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")) + RotatedAxis()

#Figure 5B
#SPEED analysis
Tert.markers <- read.delim("Tert.markers.txt", h=T)

### Signaling Pathway enrichment
Pathways <- read_tsv("speed2_signatures.tsv")
Pathways <- Pathways[grepl("UP", Pathways$regulation),]


PathwayEnrich <- data.frame(matrix(NA, ncol=length(Gene_groups)*2, nrow=length(unique(Pathways$Pathway))))
rownames(PathwayEnrich) <- unique(Pathways$Pathway)
colnames(PathwayEnrich) <- c(paste("Pval",names(Gene_groups),sep="_"), paste("Enrich",names(Gene_groups),sep="_"))

for( k in 1:length(Gene_groups)){
  tmp <- Gene_groups[[k]]
  tmp_length <- length(tmp)
  tmp_length_not <- length(unique(Tert.markers[!Tert.markers$gene %in% tmp, 'Symbol']))
  for (i in 1:length(unique(Pathways$Pathway))){
    tmp_Symbol <- Pathways[Pathways$Pathway == unique(Pathways$Pathway)[i] & Pathways$qval < 0.01,]
    tmp_Symbol <- tmp_Symbol[tmp_Symbol$SYMBOL %in% Tert.markers$gene,]
    tmp_Symbol_length <- length(tmp_Symbol$SYMBOL)
    tmp_tmp <- nrow(tmp_Symbol[tmp_Symbol$SYMBOL %in% tmp,])
    if(tmp_tmp>0){
      PathwayEnrich[i,k] <- phyper(tmp_tmp,tmp_length, tmp_length_not,tmp_Symbol_length,lower.tail=FALSE)
      PathwayEnrich[i,k+length(Gene_groups)] <- log2((tmp_tmp/tmp_Symbol_length)/(tmp_length/length(Tert.markers$gene)))
    } else{
      PathwayEnrich[i,k] <- 1
      PathwayEnrich[i,k+length(Gene_groups)] <- 0
    }
  }
}
for (k in 1:length(Gene_groups)){
  PathwayEnrich[PathwayEnrich[,k] > 0.05,k+length(Gene_groups)] <- 0
}

p <- -log10(PathwayEnrich[,1:length(Gene_groups)])
p <- p[!is.infinite(rowSums(p)),] 
# Plot the pathways of interest
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.01),max(p),length=51))

heatmap.2(as.matrix(p),main="RNASPEED", Rowv = T, Colv=F, dendrogram='none',cexCol = 0.5, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none' )

rm(i,PathwayEnrich, Pathways, col1, tmp, tmp2, tmp_down, tmp_up, down, down_not, up, up_not)

clusters <- sort(unique(Tert.markers$cluster))

Gene_groups <- lapply(clusters, function(cl){
  unique(Tert.markers$gene[Tert.markers$cluster == cl])
})

names(Gene_groups) <- clusters

Pathway_list <- unique(Pathways$Pathway)

PathwayEnrich <- data.frame(
  matrix(NA,
         ncol = length(Gene_groups)*2,
         nrow = length(Pathway_list))
)

rownames(PathwayEnrich) <- Pathway_list
colnames(PathwayEnrich) <- c(
  paste("Pval", names(Gene_groups), sep="_"),
  paste("Enrich", names(Gene_groups), sep="_")
)

all_genes <- unique(Tert.markers$gene)

for(k in seq_along(Gene_groups)){
  
  tmp <- Gene_groups[[k]]
  tmp <- tmp[!is.na(tmp)]
  
  tmp_length <- length(tmp)
  tmp_length_not <- length(setdiff(all_genes, tmp))
  
  for(i in seq_along(Pathway_list)){
    
    tmp_Symbol <- Pathways[
      Pathways$Pathway == Pathway_list[i] &
        Pathways$qval < 0.01, ]
    
    tmp_Symbol <- tmp_Symbol[tmp_Symbol$SYMBOL %in% all_genes, ]
    
    tmp_Symbol_length <- length(unique(tmp_Symbol$SYMBOL))
    
    tmp_tmp <- length(intersect(tmp, tmp_Symbol$SYMBOL))
    
    if(tmp_tmp > 0 && tmp_length > 0 && tmp_Symbol_length > 0){
      
      PathwayEnrich[i,k] <- phyper(
        tmp_tmp,
        tmp_length,
        tmp_length_not,
        tmp_Symbol_length,
        lower.tail = FALSE
      )
      
      PathwayEnrich[i,k + length(Gene_groups)] <-
        log2((tmp_tmp/tmp_Symbol_length) /
               (tmp_length/length(all_genes)))
      
    } else {
      
      PathwayEnrich[i,k] <- 1
      PathwayEnrich[i,k + length(Gene_groups)] <- 0
      
    }
    
  }
}


for(k in seq_along(Gene_groups)){
  PathwayEnrich[
    !is.na(PathwayEnrich[,k]) &
      PathwayEnrich[,k] > 0.05,
    k + length(Gene_groups)
  ] <- 0
}



max_val <- max(PathwayEnrich[is.finite(as.matrix(PathwayEnrich))], na.rm = TRUE)

# Replace Inf and -Inf
PathwayEnrich[is.infinite(as.matrix(PathwayEnrich))] <- max_val

p <- -log10(PathwayEnrich[,1:length(Gene_groups)])

p <- as.matrix(p)

max_val <- max(p[is.finite(p)], na.rm = TRUE)
p[!is.finite(p)] <- max_val


p <- p[rowSums(p, na.rm = TRUE) > 0, ]

mat_col <- c('white', designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0, seq(-log10(0.01), max(p), length = 51))

heatmap.2(
  as.matrix(p),
  main = "SPEED",
  Rowv = TRUE,
  Colv = FALSE,
  dendrogram = 'none',
  cexCol = 0.6,
  scale = 'none',
  col = mat_col,
  breaks = mat_col_breaks,
  trace = 'none'
)

#Figure 5C
###Figure 5C
# Cell: Transcriptional targets of Hippo signaling in mammalian cells
# Goal here Hippo pathway pertubation in 3 cell lines, identify NF2-LATS2 sensitive genes that are reversed by YAP-TAZ knockdown
library(GEOquery)
library(limma)

HIPPO_signature <- c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")

HIPPO_target <- c("CCN1", "PEA15", "NPPB","EPHA2", "NUAK2", "FAM107B","ANKRD1", "MYOF", "TSPAN4", "PARVA", "RBM14","CENATAC","GPRC5A","KRT7", "KRT18","HSPB8","CRY1","SAMD4A","SNAPC1","TPM1","DUSP14","LDLR",
"SYDE1", "NFKBID","CRIM1","KMT5A","RND3","AMOTL2","PIM1","CPA4","CYRIB","MIR622","UGCG","TCEAL9","PIM2","FLNA")

GSE49384_results <- read.delim("GSE49384_results.txt", h=T)


# helper function: median + shaded interval
add_median_band <- function(mat, x, col_line="black", col_fill=rgb(0,0,1,0.2)) {
  med <- apply(mat, 2, median)
  sdev <- apply(mat, 2, sd)
  n <- nrow(mat)
  se <- sdev / sqrt(n)
  
  # approximate 95% CI
  lower <- med - 1.96 * se
  upper <- med + 1.96 * se
  
  polygon(c(x, rev(x)),
          c(lower, rev(upper)),
          col = col_fill, border = NA)
  
  lines(x, med, col=col_line, lwd=2)
}

# Making line plots for the three cell lines for HIPPO_signature

y <- GSE49384_results[GSE49384_results$Symbol_canonical %in% HIPPO_signature,10:27]
colnames(y) <- paste(GSE49384_pheno$cellline, GSE49384_pheno$group, sep="_")

y <- (y[, seq(1, ncol(y), by = 2)] +
        y[, seq(2, ncol(y), by = 2)]) / 2

y <- t(scale(t(y)))

plot(0,0,pch="", xlim=c(1,9), ylim=c(-1,1))
add_median_band(y[,1:3], x=1:3, col_line="black",   col_fill="lightgrey")
add_median_band(y[,4:6], x=4:6, col_line="black",  col_fill="lightgrey")
add_median_band(y[,7:9], x=7:9, col_line="black", col_fill="lightgrey")

# Making line plots for the three cell lines for HIPPO_target signature

y <- GSE49384_results[GSE49384_results$Symbol_canonical %in% HIPPO_target,10:27]
colnames(y) <- paste(GSE49384_pheno$cellline, GSE49384_pheno$group, sep="_")

y <- (y[, seq(1, ncol(y), by = 2)] +
        y[, seq(2, ncol(y), by = 2)]) / 2

y <- t(scale(t(y)))

plot(0,0,pch="", xlim=c(1,9), ylim=c(-1,1))
add_median_band(y[,1:3], x=1:3, col_line="black",   col_fill="lightgrey")
add_median_band(y[,4:6], x=4:6, col_line="black",  col_fill="lightgrey")
add_median_band(y[,7:9], x=7:9, col_line="black", col_fill="lightgrey")


#Figure 5D & 5E
Tert <- readRDS("Tert_Subtypes.rds")
Idents(Tert) <- "Dataset"
#Vlnplot expression of Hippo compentency genes per each cell lines
genes_of_interest<- c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")

##Calculate the average expression levels of each program (cluster) on single cell level
Tert <- AddModuleScore(
  Tert,
  features = list(genes_of_interest),
  name = "ProgramScore"
)

VlnPlot(Tert, features="ProgramScore1", group.by="Dataset", pt.size=0)

#Vlnplot expression of Hippo target genes per each cell lines
hippo_targets <- c("CCN1", "PEA15", "NPPB","EPHA2", "NUAK2", "FAM107B","ANKRD1", "MYOF", "TSPAN4", "PARVA", "RBM14","CENATAC","GPRC5A","KRT7", "KRT18","HSPB8","CRY1","SAMD4A","SNAPC1","TPM1","DUSP14","LDLR",
"SYDE1", "NFKBID","CRIM1","KMT5A","RND3","AMOTL2","PIM1","CPA4","CYRIB","MIR622","UGCG","TCEAL9","PIM2","FLNA")

##Calculate the average expression levels of each program (cluster) on single cell level
Tert <- AddModuleScore(
  Tert,
  features = list(hippo_targets),
  name = "ProgramScore"
)
VlnPlot(Tert, features="ProgramScore1", group.by="Dataset", pt.size=0)

#Featureplot of Hippo compentency and target genes
Tert <- AddModuleScore(
  Tert,
  features = list(hippo_targets, genes_of_interest),
  name = c("HippoScore","SignatureScore")
)

umap <- Embeddings(Tert, reduction = "umap")

df <- data.frame(
  UMAP_1 = umap[,1],
  UMAP_2 = umap[,2],
  HippoScore1 = Tert$HippoScore1,
  SignatureScore2 = Tert$SignatureScore2
)
col_man <- rev(RColorBrewer::brewer.pal(8, "Spectral"))


p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = HippoScore1)) +
  geom_point(size = 0.2) +
  scale_color_gradientn(
    colours = colorRampPalette(col_man)(100),
    limits = quantile(df$HippoScore1, c(0.05, 0.95), na.rm = TRUE),
    name = "HippoScore"
  ) +
  theme_classic()

p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = HippoScore1)) +
  geom_point(size = 0.2) +
  scale_color_gradientn(
    colours = colorRampPalette(col_man)(100),
    limits = quantile(df$HippoScore1, c(0.05, 0.95), na.rm = TRUE),
    name = "HippoScore"
  ) +
  theme_classic()
p2 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = SignatureScore2)) +
  geom_point(size = 0.2) +
  scale_color_gradientn(
    colours = colorRampPalette(col_man)(100),
    limits = quantile(df$SignatureScore2, c(0.05, 0.95), na.rm = TRUE),
    name = "SignatureScore"
  ) +
  theme_classic()

library(patchwork)

p1 + p2

#Figure 5F
##Truli treatment
Truli_ALP <- read.xlsx("Truli_Rep1.xlsx", 
                       sheetIndex = 1, header=TRUE)

library(dplyr)
library(tidyverse)
df.summary <- Truli_ALP %>%
  group_by(Conditions) %>%
  summarise(
    sd = sd(Data, na.rm = TRUE),
    Data = mean(Data)
  )
df.summary

library(ggplot2)
# Default bar plot
conditions <- c("UT", "AD10_DMSO", "AD10_0.125uM", "DD8_DMSO", "DD8_0.125uM","CB4_DMSO", "CB4_0.125uM","CD8_DMSO", "CD8_0.125uM")

ggplot(Truli_ALP, aes(Conditions, Data)) + scale_x_discrete(limits = conditions)+ 
  geom_bar(stat = "identity", data = df.summary,
           fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.2),
               color = "black") + 
  geom_errorbar(
    aes(ymin = Data-sd, ymax = Data+sd),
    data = df.summary, width = 0.2) 
#check significance
tmp <- read.xlsx("Truli_Rep1_tmp.xlsx", 
                 sheetIndex = 1, header=TRUE)

t.test(tmp$AD10_0.125uM,tmp$AD10_DMSO)
t.test(tmp$DD8_0.125uM,tmp$DD8_DMSO)
t.test(tmp$CB4_0.125uM,tmp$CB4_DMSO)
t.test(tmp$CD8_0.125uM,tmp$CD8_DMSO)

#Celastrol treatment
Celastrol_ALP <- read.xlsx("Celastrol_Rep1.xlsx", 
                           sheetIndex = 1, header=TRUE)
library(dplyr)
df.summary <- Celastrol_ALP %>%
  group_by(Conditions) %>%
  summarise(
    sd = sd(Data, na.rm = TRUE),
    Data = mean(Data)
  )
df.summary

library(ggplot2)
# Default bar plot
conditions <- c("UT", "AD10_DMSO", "AD10_0.03uM", "DD8_DMSO", "DD8_0.03uM","CB4_DMSO", "CB4_0.03uM","CD8_DMSO", "CD8_0.03uM")

ggplot(Celastrol_ALP, aes(Conditions, Data)) + scale_x_discrete(limits = conditions)+ 
  geom_bar(stat = "identity", data = df.summary,
           fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.2),
               color = "black") + 
  geom_errorbar(
    aes(ymin = Data-sd, ymax = Data+sd),
    data = df.summary, width = 0.2) 
	
tmp <- read.xlsx("Celastrol_Rep1_tmp.xlsx", 
                 sheetIndex = 1, header=TRUE)

t.test(tmp$AD10_0.03uM,tmp$AD10_DMSO)
t.test(tmp$DD8_0.03uM,tmp$DD8_DMSO)
t.test(tmp$CB4_0.03uM,tmp$CB4_DMSO)
t.test(tmp$CD8_0.03uM,tmp$CD8_DMSO)


#Figure 5G
library("xlsx")
library(ggplot2)
library(dplyr)
library(tidyverse)

Truli_ALP <- read.xlsx("Truli.xlsx", 
                       sheetIndex = 1, header=TRUE)


df.summary <- Truli_ALP %>%
  group_by(Conditions) %>%
  summarise(
    sd = sd(Data, na.rm = TRUE),
    Data = mean(Data)
  )
df.summary


df.summary <- df.summary %>%
  mutate(
    cell_line = sub("_.*", "", Conditions),
    treatment = sub(".*_", "", Conditions)
  )


gain_df <- df.summary %>%
  filter(treatment %in% c("DMSO", "0.125uM")) %>%
  select(cell_line, treatment, Data) %>%
  pivot_wider(
    names_from = treatment,
    values_from = Data
  ) %>%
  mutate(gain = `0.125uM` / DMSO)

gain_df

gain_fixed <- gain_df %>%
  group_by(cell_line) %>%
  summarise(
    DMSO = max(DMSO, na.rm = TRUE),
    `0.125uM` = max(`0.125uM`, na.rm = TRUE),
    gain = `0.125uM` / DMSO,
    .groups = "drop"
  )

gain_fixed <- gain_fixed %>%
  mutate(
    gain_pct = (gain - 1) * 100
  )
gain_fixed$cell_line <- factor(
  gain_fixed$cell_line,
  levels = c("AD10", "DD8", "CB4", "CD8")
)
ggplot(gain_fixed, aes(x = cell_line, y = gain_pct)) +
  geom_col(fill = "grey70", color = "black") +
  theme_minimal() +
  labs(
    title = "Gain relative to DMSO",
    y = "Gain (%) vs DMSO",
    x = "Cell line"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed")

gain_fixed_truli <- gain_fixed


##Celastrol
Celastrol_ALP <- read.xlsx("Celastrol.xlsx", 
                       sheetIndex = 1, header=TRUE)
df.summary <- Celastrol_ALP %>%
  group_by(Conditions) %>%
  summarise(
    sd = sd(Data, na.rm = TRUE),
    Data = mean(Data)
  )
df.summary

df.summary <- df.summary %>%
  mutate(
    cell_line = sub("_.*", "", Conditions),
    treatment = sub(".*_", "", Conditions)
  )
gain_df <- df.summary %>%
  filter(treatment %in% c("DMSO", "0.03uM")) %>%
  select(cell_line, treatment, Data) %>%
  pivot_wider(
    names_from = treatment,
    values_from = Data
  ) %>%
  mutate(gain = `0.03uM` / DMSO)

gain_df

gain_fixed <- gain_df %>%
  group_by(cell_line) %>%
  summarise(
    DMSO = max(DMSO, na.rm = TRUE),
    `0.03uM` = max(`0.03uM`, na.rm = TRUE),
    gain = `0.03uM` / DMSO,
    .groups = "drop"
  )

gain_fixed <- gain_fixed %>%
  mutate(
    gain_pct = (gain - 1) * 100
  )
gain_fixed$cell_line <- factor(
  gain_fixed$cell_line,
  levels = c("AD10", "DD8", "CB4", "CD8")
)
ggplot(gain_fixed, aes(x = cell_line, y = gain_pct)) +
  geom_col(fill = "grey70", color = "black") +
  theme_minimal() +
  labs(
    title = "Loss relative to DMSO",
    y = "Loss (%) vs DMSO",
    x = "Cell line"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed")

gain_fixed_cel   <- gain_fixed

##plot gain vs loss

gain_merged <- gain_fixed_truli %>%
  select(cell_line, truli_gain_pct = gain_pct) %>%
  left_join(
    gain_fixed_cel %>%
      select(cell_line, cel_loss_pct = gain_pct),
    by = "cell_line"
  )

gain_merged


ggplot(gain_merged,
       aes(x = truli_gain_pct, y = cel_loss_pct, label = cell_line)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 4) +
  
  # zero reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # expected opposing-response line
  geom_abline(intercept = 0, slope = -1,
              linetype = "dotted", color = "black") +
  
  theme_minimal() +
  labs(
    x = "Gain with Truli (% vs DMSO)",
    y = "Loss with Celastrol (% vs DMSO)",
    title = "Opposing effects of Hippo pathway activation and inhibition on ALP activity"
  )



# Figure 5H

HIPPO_signature <- c("NTRK2","TEAD1","PLCB1","PDGFRB","GNAS","TCF7L2","TCF7L1","PRKAR1A","PRKCE","IGF1R","GNAI2","PRKCH","SMAD3","CDH6","EGFR","SMAD2","GNAQ","SAV1","MAPK10","MAP4K3","CTNNA1","YWHAQ")

HIPPO_target <- c("CCN1", "PEA15", "NPPB","EPHA2", "NUAK2", "FAM107B","ANKRD1", "MYOF", "TSPAN4", "PARVA", "RBM14","CENATAC","GPRC5A","KRT7", "KRT18","HSPB8","CRY1","SAMD4A","SNAPC1","TPM1","DUSP14","LDLR",
"SYDE1", "NFKBID","CRIM1","KMT5A","RND3","AMOTL2","PIM1","CPA4","CYRIB","MIR622","UGCG","TCEAL9","PIM2","FLNA")


# TERT cells
final.TERT <- read.delim("final.TERT.txt",h=T)
# Update annotation of TERT object
library(org.Hs.eg.db)
library(AnnotationDbi)

alias_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(final.TERT$SYMBOL),
  columns = c("SYMBOL"),
  keytype = "ALIAS"
)
  
alias_map <- alias_map[!duplicated(alias_map$ALIAS), ]
  
final.TERT$Symbol_canonical <- alias_map$SYMBOL[match(final.TERT$SYMBOL, alias_map$ALIAS)]
final.TERT$Symbol_canonical[is.na(final.TERT$Symbol_canonical)] <- final.TERT$SYMBOL[is.na(final.TERT$Symbol_canonical)]

length(unique(final.TERT$SYMBOL))
length(unique(final.TERT$Symbol_canonical))

boxplot(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob14d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob14d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob7d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob7d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob3d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob3d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob1d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob1d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob4h"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob4h"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad4h"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad4h"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad1d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad1d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad3d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad3d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad7d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad7d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad14d"],
        final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad14d"],
        col=c('blue','lightblue','blue','lightblue','blue','lightblue','blue','lightblue','blue','lightblue','red','salmon','red','salmon','red','salmon','red','salmon','red','salmon')
        )
abline(h=0, lty=2)
legend(
  "topright",
  legend = c("Target (Osteogenic)", "Signature (Osteogenic)",
             "Target (Adipogenic)", "Signature (Adipogenic)"),
  fill = c("blue", "lightblue", "red", "salmon"),
  border = "black",
  cex = 0.8
)


wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob14d"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob14d"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob7d"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob7d"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob3d"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob3d"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob1d"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob1d"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ob4h"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ob4h"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad4h"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad4h"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad1d"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad1d"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad3d"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad3d"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad7d"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad7d"])
wilcox.test(final.TERT[final.TERT$Symbol_canonical %in% HIPPO_target,"log2FC_Ad14d"],
            final.TERT[final.TERT$Symbol_canonical %in% HIPPO_signature,"log2FC_Ad14d"])


#Figure 5I
#load the following object
library(R.utils)
GSE253355 <- readRDS("GSE253355_MSC_Subset_Seurat.rds")
DimPlot(GSE253355)

ggplot(df, aes(mscumapdim50_1, mscumapdim50_2, color = CytoTRACE_score)) +
  geom_point(size = 0.15) +
  scale_color_gradientn(
    colours = colorRampPalette(col_man)(100),
    limits = range(df$CytoTRACE_score, na.rm = TRUE),
    name = "CytoTRACE_score"
  ) +
  theme_classic()



#Figure 5J
# Adding HIPPO_target and HIPPO signature to the object
GSE253355 <- AddModuleScore(
  object = GSE253355,
  features = list(HIPPO_target),
  name = "HIPPO_target"
)

GSE253355 <- AddModuleScore(
  object = GSE253355,
  features = list(HIPPO_signature),
  name = "HIPPO_signature"
)

VlnPlot(GSE253355, "HIPPO_signature1", pt.size = 0)
VlnPlot(GSE253355, "HIPPO_target1", pt.size = 0)

#Figure 5K
col_man <- c(rev(RColorBrewer::brewer.pal(8,"Spectral"))[8])

library(ggplot2)

df <- FetchData(
  GSE253355,
  vars = c("mscumapdim50_1", "mscumapdim50_2", "HIPPO_signature1","HIPPO_target1","CytoTRACE_score")
)

ggplot(df, aes(mscumapdim50_1, mscumapdim50_2, color = HIPPO_target1)) +
  geom_point(size = 0.15) +
  scale_color_gradientn(
    colours = colorRampPalette(col_man)(100),
    limits = range(df$HIPPO_target1, na.rm = TRUE),
    name = "HIPPO_target1"
  ) +
  theme_classic()

ggplot(df, aes(mscumapdim50_1, mscumapdim50_2, color = Atenisa1)) +
  geom_point(size = 0.15) +
  scale_color_gradientn(
    colours = colorRampPalette(col_man)(100),
    limits = range(df$Atenisa1, na.rm = TRUE),
    name = "HIPPO_signature1"
  ) +
  theme_classic()
  
rm(df,p,col_man)
# Final aesthetics were done in Illustrator

