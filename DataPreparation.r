##Data preparation
######## Get packages for R

conda activate r-4.0
R
# Get all the libraries

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
BiocManager::install("simpleSingleCell")
BiocManager::install("BiocSingular")
BiocManager::install("umap")
BiocManager::install("Seurat")
BiocManager::install("Matrix")
BiocManager::install("scran")
BiocManager::install("scater")
BiocManager::install("DropletUtils")
BiocManager::install("batchelor")
BiocManager::install("circlize")
BiocManager::install("biomaRt")
BiocManager::install("MAST")
BiocManager::install("limma")
BiocManager::install("RANN")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("dplyr")
BiocManager::install("clusterProfiler")
BiocManager::install("rWikiPathways")
BiocManager::install("GEOquery")
BiocManager::install("edgeR")
BiocManager::install("glmnet")
BiocManager::install("phateR")
BiocManager::install("pcaMethods")
BiocManager::install("TSCAN")
BiocManager::install("tradeSeq")
BiocManager::install("scDblFinder")
install.packages("seriation")

BiocManager::install("devtools")
library(devtools)
install_github("immunogenomics/harmony")
install_github('theislab/kBET')
install_github("Albluca/distutils")
install_github("Albluca/ElPiGraph.R")
install_github("jokergoo/ComplexHeatmap")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install_github('cole-trapnell-lab/leidenbase')
install_github("immunogenomics/lisi")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnsDb.Hsapiens.v79")

q()

######## Analysis in R
# Start R from Terminal

conda activate r-4.0
R
# Loading packages for analysis
library(CellChat) #load Cell Chat first and then the other packages
library(SingleCellExperiment) 
library(simpleSingleCell) 
library(BiocSingular) 
library(umap) 
library(Seurat) 
library(Matrix) 
library(scran) 
library(scater) 
library(DropletUtils) 
library(batchelor)
library(harmony)
library(ComplexHeatmap)
library(circlize)
library(MAST)
library(limma)
library(RANN)
library(biomaRt)
library(kBET)
library(lisi)
library(org.Mm.eg.db)
library(dplyr)
library(clusterProfiler)
library(rWikiPathways)
library(GEOquery)
library(edgeR)
library(glmnet)
library(velociraptor)
library(phateR)
library(ElPiGraph.R)
library(TSCAN)
library(tradeSeq)
library(seriation)
library(scuttle)
library(mvoutlier)
library(dplyr)
library(scDblFinder)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(DEFormats)
library(harmony)

AD10_data <- readRDS("AD10.dgecounts.rds")
CB4_data <- readRDS("CB4.dgecounts.rds")
CD_data <- readRDS("CD8.dgecounts.rds")
CD8_old_data <- readRDS("CD8_old.dgecounts.rds")
DD8_data <- readRDS("DD8.dgecounts.rds")
Msc_BM_data <- readRDS("TERT_Msc.dgecounts.rds")

# Extracting intron+exons counts 
AD10_data <- AD10_data$umicount$inex$all
CB4_data <- CB4_data$umicount$inex$all
CD_data <- CD_data$umicount$inex$all
CD8_old_data <- CD8_old_data$umicount$inex$all
DD8_data <- DD8_data$umicount$inex$all
Msc_BM_data <- Msc_BM_data$umicount$inex$all


#Delete rows with 0 values
AD10_data <- AD10_data[which(rowSums(AD10_data)>0),]
CB4_data <- CB4_data[which(rowSums(CB4_data)>0),]
CD_data <- CD_data[which(rowSums(CD_data)>0),]
CD8_old_data <- CD8_old_data[which(rowSums(CD8_old_data)>0),]
DD8_data <- DD8_data[which(rowSums(DD8_data)>0),]
Msc_BM_data <- Msc_BM_data[which(rowSums(Msc_BM_data)>0),]

#Delete columns with less than 2
AD10_data <- AD10_data[,which(colSums(AD10_data)>2)] 
CB4_data <- CB4_data[,which(colSums(CB4_data)>2)]
CD_data <- CD_data[,which(colSums(CD_data)>2)]
CD8_old_data <- CD8_old_data[,which(colSums(CD8_old_data)>2)]
DD8_data <- DD8_data[,which(colSums(DD8_data)>2)]
Msc_BM_data <- Msc_BM_data[,which(colSums(Msc_BM_data)>2)]





## Fill in the matrices to give them all the same dimensions
## RATIONALE: Genes with 0 counts across all barcodes in a particular experiment are left out from zUMIs.
# Find all non-zero genes in all conditions
Genes <- unique(c(rownames(AD10_data),rownames(CB4_data),rownames(CD_data),rownames(CD8_old_data),rownames(DD8_data),rownames(Msc_BM_data)))

# Insert empty lines with missing genes to get same dimensions
Tmp <- as.sparse(matrix(ncol=ncol(AD10_data), nrow=length(Genes[!(Genes %in% rownames(AD10_data))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(AD10_data))]
AD10_data <- Matrix::rbind2(AD10_data, Tmp)
AD10_data <- AD10_data[ order(rownames(AD10_data)),]
dim(AD10_data) 
# 43631  6654

Tmp <- as.sparse(matrix(ncol=ncol(CB4_data), nrow=length(Genes[!(Genes %in% rownames(CB4_data))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(CB4_data))]
CB4_data <- Matrix::rbind2(CB4_data, Tmp)
CB4_data <- CB4_data[ order(rownames(CB4_data)),]
dim(CB4_data) 
# 43631 10908

Tmp <- as.sparse(matrix(ncol=ncol(CD_data), nrow=length(Genes[!(Genes %in% rownames(CD_data))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(CD_data))]
CD_data <- Matrix::rbind2(CD_data, Tmp)
CD_data <- CD_data[ order(rownames(CD_data)),]
dim(CD_data) 
# 43631  3962

Tmp <- as.sparse(matrix(ncol=ncol(CD8_old_data), nrow=length(Genes[!(Genes %in% rownames(CD8_old_data))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(CD8_old_data))]
CD8_old_data <- Matrix::rbind2(CD8_old_data, Tmp)
CD8_old_data <- CD8_old_data[ order(rownames(CD8_old_data)),]
dim(CD8_old_data) 
# 43631 38736

Tmp <- as.sparse(matrix(ncol=ncol(DD8_data), nrow=length(Genes[!(Genes %in% rownames(DD8_data))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(DD8_data))]
DD8_data <- Matrix::rbind2(DD8_data, Tmp)
DD8_data <- DD8_data[ order(rownames(DD8_data)),]
dim(DD8_data) 
#43631  7256


Tmp <- as.sparse(matrix(ncol=ncol(Msc_BM_data), nrow=length(Genes[!(Genes %in% rownames(Msc_BM_data))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(Msc_BM_data))]
Msc_BM_data <- Matrix::rbind2(Msc_BM_data, Tmp)
Msc_BM_data <- Msc_BM_data[ order(rownames(Msc_BM_data)),]
dim(Msc_BM_data)
# 43631  2597


## Paste in the experiment name into the column names to make barcodes/cell IDs unique.
colnames(AD10_data) <- paste("AD10_",colnames(AD10_data), sep="")
colnames(CB4_data) <- paste("CB_",colnames(CB4_data), sep="")
colnames(CD_data) <- paste("CD_",colnames(CD_data), sep="")
colnames(CD8_old_data) <- paste("CD8_",colnames(CD8_old_data), sep="")
colnames(DD8_data) <- paste("DD8_",colnames(DD8_data), sep="")
colnames(Msc_BM_data) <- paste("Msc_",colnames(Msc_BM_data), sep="")

##merge the data from te CD8 (new and previous sequencing)
CD_data <- cbind(CD_data, CD8_old_data)
colnames(CD_data) <- paste("CD_",colnames(CD_data), sep="")
# Process gene list (downloaded from BioMart), keep only non-empty gene symbols and deduplicate.
Genes <- read.delim("mart_export_human_Ens_GRCh38.txt")
Genes <- Genes[Genes$Gene.stable.ID %in% rownames(AD10_data),]
Genes <- Genes[Genes$Gene.name !=  "",]
Genes <- Genes[ duplicated(Genes$Gene.stable.ID) == F,]
Genes <- Genes[ duplicated(Genes$Gene.name) == F,]
Genes <- Genes[ order(Genes$Gene.stable.ID),]
dim(Genes)
# 41475     5

# Subset count matrices and replace Ensemble IDs with gene symbols
AD10_data <- AD10_data[ rownames(AD10_data) %in% Genes$Gene.stable.ID,]
AD10_data <- AD10_data[ order(rownames(AD10_data)),]
rownames(AD10_data) <- as.character(Genes$Gene.name)
dim(AD10_data)
#41475   6654


CB4_data <- CB4_data[ rownames(CB4_data) %in% Genes$Gene.stable.ID,]
CB4_data <- CB4_data[ order(rownames(CB4_data)),]
rownames(CB4_data) <- as.character(Genes$Gene.name)
dim(CB4_data)
#41475  10908

CD_data <- CD_data[ rownames(CD_data) %in% Genes$Gene.stable.ID,]
CD_data <- CD_data[ order(rownames(CD_data)),]
rownames(CD_data) <- as.character(Genes$Gene.name)
dim(CD_data)
#41475   42698

DD8_data <- DD8_data[ rownames(DD8_data) %in% Genes$Gene.stable.ID,]
DD8_data <- DD8_data[ order(rownames(DD8_data)),]
rownames(DD8_data) <- as.character(Genes$Gene.name)
dim(DD8_data)
# 41475   7256

Msc_BM_data <- Msc_BM_data[ rownames(Msc_BM_data) %in% Genes$Gene.stable.ID,]
Msc_BM_data <- Msc_BM_data[ order(rownames(Msc_BM_data)),]
rownames(Msc_BM_data) <- as.character(Genes$Gene.name)
dim(Msc_BM_data)
# 41475   2597

# Check mito genes
length(grep("MT-",rownames(AD10_data)))
#37
length(grep("MT-",rownames(CB4_data)))
#37
length(grep("MT-",rownames(CD_data)))
#37
length(grep("MT-",rownames(DD8_data)))
#37
length(grep("MT-",rownames(Msc_BM_data)))
#37

# Check proliferation genes
par(mfrow=c(5,3))
tmp <- AD10_data
hist(tmp[rownames(tmp)=="MKI67",], main="AD10 MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="TOP2A",], main="AD10 TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="PCNA",], main="AD10 PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))
tmp <- CB4_data
hist(tmp[rownames(tmp)=="MKI67",], main="CB MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="TOP2A",], main="CB TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="PCNA",], main="CB PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))
tmp <- CD_data
hist(as.numeric(tmp[rownames(tmp)=="MKI67",]), main="CD MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(as.numeric(tmp[rownames(tmp)=="TOP2A",]), main="CD TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(as.numeric(tmp[rownames(tmp)=="PCNA",]), main="CD PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))
tmp <- DD8_data
hist(tmp[rownames(tmp)=="MKI67",], main="DD8 MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="TOP2A",], main="DD8 TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="PCNA",], main="DD8 PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))
tmp <- Msc_BM_data
hist(tmp[rownames(tmp)=="MKI67",], main="Msc MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="TOP2A",], main="Msc TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="PCNA",], main="Msc PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))



# Create SingleCellExperiment objects to use DropUtils, scran, scater, etc.
AD10_data_sce <- SingleCellExperiment(list(counts=AD10_data))
CB_data_sce <- SingleCellExperiment(list(counts=CB4_data))
CD_data_sce <- SingleCellExperiment(list(counts=CD_data))
DD8_data_sce <- SingleCellExperiment(list(counts=DD8_data))
Msc_BM_data_sce <- SingleCellExperiment(list(counts=Msc_BM_data))

# Testing proliferative genes
par(mfrow=c(5,4))
tmp <- counts(AD10_data_sce)
hist(tmp[rownames(tmp)=="MKI67",], main="AD10 MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="TOP2A",], main="AD10 TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="PCNA",], main="AD10 PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))
tmp <- counts(CB_data_sce)
hist(tmp[rownames(tmp)=="MKI67",], main="CB MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="TOP2A",], main="CB TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="PCNA",], main="CB PCNA", ylim=c(0,50), breaks=seq(00,100,length=25))
tmp <- counts(CD_data_sce)
hist(as.numeric(tmp[rownames(tmp)=="MKI67",]), main="CD MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(as.numeric(tmp[rownames(tmp)=="TOP2A",]), main="CD TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(as.numeric(tmp[rownames(tmp)=="PCNA",]), main="CD PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))
tmp <- counts(DD8_data_sce)
hist(tmp[rownames(tmp)=="MKI67",], main="DD8 MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="TOP2A",], main="DD8 TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="PCNA",], main="DD8 PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))
tmp <- counts(Msc_BM_data_sce)
hist(tmp[rownames(tmp)=="MKI67",], main="Msc MKI67", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="TOP2A",], main="Msc TOP2A", ylim=c(0,50), breaks=seq(0,100,length=25))
hist(tmp[rownames(tmp)=="PCNA",], main="Msc PCNA", ylim=c(0,50), breaks=seq(0,100,length=25))


## Find empty droplets. NOTE: THIS STEP IS NON-DETERMINSTIC - RESULTS VARY FROM RUN TO RUN

# Exclude mitochondial and ribosomal genes for this part of the analysis, but do not remove them from the objects
mito_ribo <- Genes[c(grep("RPL",Genes$Gene.type),grep("MT-",Genes$Gene.type)),'Gene.name']

##  Observational plots
par(mfrow=c(5,2))

#AD10
br.out <- barcodeRanks(AD10_data)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="AD10")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(AD10_data_sce))
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

br.out <- barcodeRanks(AD10_data[!rownames(AD10_data) %in% mito_ribo,])
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="AD10 filtered")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(AD10_data_sce[!rownames(AD10_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

#CB
br.out <- barcodeRanks(CB4_data)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="CB")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(CB_data_sce))
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

br.out <- barcodeRanks(CB4_data[!rownames(CB4_data) %in% mito_ribo,])
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="CB filtered")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(CB_data_sce[!rownames(CB_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

# CD	
br.out <- barcodeRanks(CD_data)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="CD")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(CD_data_sce))
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

br.out <- barcodeRanks(CD_data[!rownames(CD_data) %in% mito_ribo,])
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="CD filtered")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(CD_data_sce[!rownames(CD_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

#DD8	
br.out <- barcodeRanks(DD8_data)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="DD8")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(DD8_data_sce))
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

br.out <- barcodeRanks(DD8_data[!rownames(DD8_data) %in% mito_ribo,])
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="DD8 filtered")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(DD8_data_sce[!rownames(DD8_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

# Msc	
br.out <- barcodeRanks(Msc_BM_data)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="Msc")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(Msc_BM_data_sce))
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)

br.out <- barcodeRanks(Msc_BM_data[!rownames(Msc_BM_data) %in% mito_ribo,])
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total",main="Msc filtered")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))

e.out <- emptyDrops(counts(Msc_BM_data_sce[!rownames(Msc_BM_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
abline(v=table(is.cell)[2])
title(paste("Cells:",table(is.cell)[2]),line=0)
is.cell <- e.out$FDR <= 0.001
abline(v=table(is.cell)[2])
title(paste("Cells FDR < 0.001:",table(is.cell)[2]),line=1)


	
## Threshold was taken after correction for wrong knee points using reatin = Inf and removal of ribosomal and mitochondial genes

par(mfrow=c(2,3), pty="s")

e.out <- emptyDrops(counts(AD10_data_sce[!rownames(AD10_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability",main="AD10")
AD10_data_sce <- AD10_data_sce[,which(e.out$FDR <= 0.01)]

e.out <- emptyDrops(counts(CB_data_sce[!rownames(CB_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability",main="CB4")
CB_data_sce <- CB_data_sce[,which(e.out$FDR <= 0.01)]

e.out <- emptyDrops(counts(CD_data_sce[!rownames(CD_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability",main="CD8")
CD_data_sce <- CD_data_sce[,which(e.out$FDR <= 0.01)]

e.out <- emptyDrops(counts(DD8_data_sce[!rownames(DD8_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability",main="DD8")
DD8_data_sce <- DD8_data_sce[,which(e.out$FDR <= 0.01)]

e.out <- emptyDrops(counts(Msc_BM_data_sce[!rownames(Msc_BM_data_sce) %in% mito_ribo,]), retain = Inf)
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability",main="Msc_BM")
Msc_BM_data_sce <- Msc_BM_data_sce[,which(e.out$FDR <= 0.01)]

dim(AD10_data_sce)
#41475  6648
dim(CB_data_sce)
#41475 10904
dim(CD_data_sce)
#41475 42693
dim(DD8_data_sce)
#41475  7238
dim(Msc_BM_data_sce)
#41475  2594

rm(tmp1,tmp2,br.out,Tmp,Genes,o,e.out,is.cell)

## Calculate QC parameters (throws a warning, that can be ignored)

AD10_data_sce <- addPerCellQC(AD10_data_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(AD10_data_sce), value = GPSLSE)))
CB_data_sce <- addPerCellQC(CB_data_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(CB_data_sce), value = GPSLSE)))
CD_data_sce <- addPerCellQC(CD_data_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(CD_data_sce), value = GPSLSE)))
DD8_data_sce <- addPerCellQC(DD8_data_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(DD8_data_sce), value = GPSLSE)))
Msc_BM_data_sce <- addPerCellQC(Msc_BM_data_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(Msc_BM_data_sce), value = GPSLSE)))

### Histograms of quality measures for each dataset

par(mfcol=c(5,4),pty="s")
hist(AD10_data_sce$subsets_Mito_percent,main="AD10", breaks=100)
hist(CB_data_sce$subsets_Mito_percent,main="CB4", breaks=100)
hist(CD_data_sce$subsets_Mito_percent,main="CD8", breaks=100)
hist(DD8_data_sce$subsets_Mito_percent,main="DD8", breaks=100)
hist(Msc_BM_data_sce$subsets_Mito_percent,main="Msc_BM", breaks=100)

hist(AD10_data_sce$sum,main="AD10", breaks=100)
hist(CB_data_sce$sum,main="CB4", breaks=100)
hist(CD_data_sce$sum,main="CD8", breaks=100)
hist(DD8_data_sce$sum,main="DD8", breaks=100)
hist(Msc_BM_data_sce$sum,main="Msc_BM", breaks=100)

hist(AD10_data_sce$detected,main="AD10", breaks=100)
hist(CB_data_sce$detected,main="CB4", breaks=100)
hist(CD_data_sce$detected,main="CD8", breaks=100)
hist(DD8_data_sce$detected,main="DD8", breaks=100)
hist(Msc_BM_data_sce$detected,main="Msc_BM", breaks=100)

hist(AD10_data_sce$sum/AD10_data_sce$detected,main="AD10", breaks=100)
hist(CB_data_sce$sum/CB_data_sce$detected,main="CB4", breaks=100)
hist(CD_data_sce$sum/CD_data_sce$detected,main="CD8", breaks=100)
hist(DD8_data_sce$sum/DD8_data_sce$detected,main="DD8", breaks=100)
hist(Msc_BM_data_sce$sum/Msc_BM_data_sce$detected,main="Msc_BM", breaks=100)

par(mfcol=c(5,4),pty="s")
hist(AD10_data_sce$subsets_Mito_percent,main="AD10")
hist(CB_data_sce$subsets_Mito_percent,main="CB4")
hist(CD_data_sce$subsets_Mito_percent,main="CD8")
hist(DD8_data_sce$subsets_Mito_percent,main="DD8")
hist(Msc_BM_data_sce$subsets_Mito_percent,main="Msc_BM")

hist(AD10_data_sce$sum,main="AD10")
hist(CB_data_sce$sum,main="CB4")
hist(CD_data_sce$sum,main="CD8")
hist(DD8_data_sce$sum,main="DD8")
hist(Msc_BM_data_sce$sum,main="Msc_BM")

hist(AD10_data_sce$detected,main="AD10")
hist(CB_data_sce$detected,main="CB4")
hist(CD_data_sce$detected,main="CD8")
hist(DD8_data_sce$detected,main="DD8")
hist(Msc_BM_data_sce$detected,main="Msc_BM")

hist(AD10_data_sce$sum/AD10_data_sce$detected,main="AD10")
hist(CB_data_sce$sum/CB_data_sce$detected,main="CB4")
hist(CD_data_sce$sum/CD_data_sce$detected,main="CD8")
hist(DD8_data_sce$sum/DD8_data_sce$detected,main="DD8")
hist(Msc_BM_data_sce$sum/Msc_BM_data_sce$detected,main="Msc_BM")


saveRDS(AD10_data_sce, "AD10_data_sce.rds")
saveRDS(CB_data_sce, "CB_data_sce.rds")
saveRDS(CD_data_sce, "CD_data_sce.rds")
saveRDS(DD8_data_sce, "DD8_data_sce.rds")
saveRDS(Msc_BM_data_sce, "Msc_BM_data_sce.rds")


## Threshold filtering of droplets in each dataset
## REMOVE: Droplets with more than 10% mitochondrial reads, less than 1000 UMIs, less than 500 genes or extremely high ratio between counts and genes (low complexity))
AD10_data_sce <- AD10_data_sce[,!(AD10_data_sce$subsets_Mito_percent > 15)] 
AD10_data_sce <- AD10_data_sce[,!(AD10_data_sce$sum/AD10_data_sce$detected > 6)] 
AD10_data_sce <- AD10_data_sce[,(AD10_data_sce$sum >= 1000 & AD10_data_sce$detected >= 1500)] 

CB_data_sce <- CB_data_sce[,!(CB_data_sce$subsets_Mito_percent > 15)] 
CB_data_sce <- CB_data_sce[,!(CB_data_sce$sum/CB_data_sce$detected > 6)] 
CB_data_sce <- CB_data_sce[,(CB_data_sce$sum >= 1000 & CB_data_sce$detected >= 1500)] 

CD_data_sce <- CD_data_sce[,!(CD_data_sce$subsets_Mito_percent > 15)] 
CD_data_sce <- CD_data_sce[,!(CD_data_sce$sum/CD_data_sce$detected > 6)] 
CD_data_sce <- CD_data_sce[,(CD_data_sce$sum >= 1000 & CD_data_sce$detected >= 1500)] 

DD8_data_sce <- DD8_data_sce[,!(DD8_data_sce$subsets_Mito_percent > 15)] 
DD8_data_sce <- DD8_data_sce[,!(DD8_data_sce$sum/DD8_data_sce$detected > 6)] 
DD8_data_sce <- DD8_data_sce[,(DD8_data_sce$sum >= 1000 & DD8_data_sce$detected >= 1500)] 

Msc_BM_data_sce <- Msc_BM_data_sce[,!(Msc_BM_data_sce$subsets_Mito_percent > 15)] 
Msc_BM_data_sce <- Msc_BM_data_sce[,!(Msc_BM_data_sce$sum/Msc_BM_data_sce$detected > 6)] 
Msc_BM_data_sce <- Msc_BM_data_sce[,(Msc_BM_data_sce$sum >= 1000 & Msc_BM_data_sce$detected >= 1500)] 

dim(AD10_data_sce)
#41475  5769
dim(CB_data_sce)
#41475 10859
dim(CD_data_sce)
#41475 40314
dim(DD8_data_sce)
#41475  6472
dim(Msc_BM_data_sce)
#41475  2404

#1500 threshold for detected genes
dim(AD10_data_sce)
#41475  4709
dim(CB_data_sce)
#41475  7879
dim(CD_data_sce)
#41475 3532
dim(DD8_data_sce)
#41475  5357
dim(Msc_BM_data_sce)
#41475  2246

## Automatic filtering droplets in each dataset using PCA across all QC metrics
# Calculate outliers
AD10_data_sce <- logNormCounts(AD10_data_sce)
CB_data_sce <- logNormCounts(CB_data_sce)
CD_data_sce <- logNormCounts(CD_data_sce)
DD8_data_sce <- logNormCounts(DD8_data_sce)
Msc_BM_data_sce <- logNormCounts(Msc_BM_data_sce)

AD10_data_sce <- runColDataPCA(AD10_data_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)
CB_data_sce <- runColDataPCA(CB_data_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)
CD_data_sce <- runColDataPCA(CD_data_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)
DD8_data_sce <- runColDataPCA(DD8_data_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)
Msc_BM_data_sce <- runColDataPCA(Msc_BM_data_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)

table(AD10_data_sce$outlier) 
# GPSLSE  TRUE
# 5757    12
# 4656    53 <1500detected genes theshhold
table(CB_data_sce$outlier)
#GPSLSE  TRUE
#10769    90
#7787    92

table(CD_data_sce$outlier)
#GPSLSE  TRUE
#39709   605
#3117   415

table(DD8_data_sce$outlier)
#GPSLSE  TRUE
#6447    25
#5272    85

table(Msc_BM_data_sce$outlier)
#GPSLSE  TRUE
# 2376    28
#2080   166


# Remove outliers with PCA
AD10_data_sce <- AD10_data_sce[ ,!AD10_data_sce$outlier] 
CB_data_sce <- CB_data_sce[ ,!CB_data_sce$outlier] 
CD_data_sce <- CD_data_sce[ ,!CD_data_sce$outlier] 
DD8_data_sce <- DD8_data_sce[ ,!DD8_data_sce$outlier] 
Msc_BM_data_sce <- Msc_BM_data_sce[ ,!Msc_BM_data_sce$outlier] 

## Threshold filtering of genes in each dataset
## REMOVE: Genes expressed in less than 10 nuclei in all datasets
# Find lowly expressed genes and get the intersection
AD10_low <- names(which(nexprs(AD10_data_sce, byrow=T) <= 1))
CB_low <- names(which(nexprs(CB_data_sce, byrow=T) <= 1))
CD_low <- names(which(nexprs(CD_data_sce, byrow=T) <= 1))
DD8_low <- names(which(nexprs(DD8_data_sce, byrow=T) <= 1))
Msc_BM_low <- names(which(nexprs(Msc_BM_data_sce, byrow=T) <= 1))
Low <- Reduce(intersect, list(AD10_low,CB_low,CD_low,DD8_low,Msc_BM_low))

# Remove lowly expressed genes
AD10_data_sce <- AD10_data_sce[ which(!(rownames(AD10_data_sce) %in% Low)),] 
CB_data_sce <- CB_data_sce[ which(!(rownames(CB_data_sce) %in% Low)),] 
CD_data_sce <- CD_data_sce[ which(!(rownames(CD_data_sce) %in% Low)),] 
DD8_data_sce <- DD8_data_sce[ which(!(rownames(DD8_data_sce) %in% Low)),] 
Msc_BM_data_sce <- Msc_BM_data_sce[ which(!(rownames(Msc_BM_data_sce) %in% Low)),] 


# Save sce and seurat objects and remove them
saveRDS(AD10_data_sce, "AD10_data_sce.rds")
saveRDS(CB_data_sce, "CB_data_sce.rds")
saveRDS(CD_data_sce, "CD_data_sce.rds")
saveRDS(DD8_data_sce, "DD8_data_sce.rds")
saveRDS(Msc_BM_data_sce, "Msc_BM_data_sce.rds")


## Filtering genes based on biotype and transcript level support (not done yet)
## REMOVE: Non-protein coding genes, keep the one with high transcript level support (tsl)
Genes <- read.delim("mart_export_human_Ens_GRCh38.txt")
Genes <- Genes[ Genes$Gene.name %in% rownames(AD10_data_sce),]
Genes <- Genes[ c(grep("tsl1", Genes$Transcript.support.level..TSL.), grep("tsl2", Genes$Transcript.support.level..TSL.)),]
Genes <- Genes[ Genes$Gene.name !="" ,]
Genes <- Genes[ Genes$Transcript.type == "protein_coding",]
Genes <- Genes[ Genes$Gene.type == "protein_coding",]
Genes <- Genes[ !is.na(Genes$Gene.name),]
Genes <- Genes[ duplicated(Genes$Gene.stable.ID) == F,]

# Filter the sce objects 
AD10_data_sce <- AD10_data_sce[ which(rownames(AD10_data_sce) %in% Genes$Gene.name),] #15665 genes left
CB_data_sce <- CB_data_sce[ which(rownames(CB_data_sce) %in% Genes$Gene.name),] # 15665 genes left
CD_data_sce <- CD_data_sce[ which(rownames(CD_data_sce) %in% Genes$Gene.name),] # 15665 genes left
DD8_data_sce <- DD8_data_sce[ which(rownames(DD8_data_sce) %in% Genes$Gene.name),] # 15665 genes left
Msc_BM_data_sce <- Msc_BM_data_sce[ which(rownames(Msc_BM_data_sce) %in% Genes$Gene.name),] # 15665 genes left

dim(AD10_data_sce)
# 15665  5757
dim(CB_data_sce)
# 15665 10769
dim(CD_data_sce)
# 15665  39709
dim(DD8_data_sce)
#15665  6447
dim(Msc_BM_data_sce)
#15665  2376

#1500 theshold
dim(AD10_data_sce)
# 15601  4656
dim(CB_data_sce)
# 15601  7787
dim(CD_data_sce)
# 15601  3117
dim(DD8_data_sce)
# 15601  5272
dim(Msc_BM_data_sce)
#15601  2080


## Normalize the count matrices # This is how GPSr I am
# Cluster each data. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
AD10_data_clusters <- quickCluster(AD10_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())
CB_data_clusters <- quickCluster(CB_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())
CD_data_clusters <- quickCluster(CD_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())
DD8_data_clusters <- quickCluster(DD8_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())
Msc_BM_data_clusters <- quickCluster(Msc_BM_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())

# Compute scaling GPSctors
AD10_data_sce <- computeSumGPSctors(AD10_data_sce, min.mean=0.1, cluster=AD10_data_clusters)
CB_data_sce <- computeSumGPSctors(CB_data_sce, min.mean=0.1, cluster=CB_data_clusters)
CD_data_sce <- computeSumGPSctors(CD_data_sce, min.mean=0.1, cluster=CD_data_clusters)
DD8_data_sce <- computeSumGPSctors(DD8_data_sce, min.mean=0.1, cluster=DD8_data_clusters)
Msc_BM_data_sce <- computeSumGPSctors(Msc_BM_data_sce, min.mean=0.1, cluster=Msc_BM_data_clusters)

# Normalize the counts
AD10_data_sce <- logNormCounts(AD10_data_sce)
CB_data_sce <- logNormCounts(CB_data_sce)
CD_data_sce <- logNormCounts(CD_data_sce)
DD8_data_sce <- logNormCounts(DD8_data_sce)
Msc_BM_data_sce <- logNormCounts(Msc_BM_data_sce)

## Calculate doublet scores. NOTE: THIS STEP IS NON-DETERMINISTI - RESULTS VARY FROM RUN TO RUN
AD10_data_sce$DoubletScore <- scDblFinder::computeDoubletDensity(AD10_data_sce, BSPARAM=IrlbaParam())
CB_data_sce$DoubletScore <- scDblFinder::computeDoubletDensity(CB_data_sce, BSPARAM=IrlbaParam())
CD_data_sce$DoubletScore <- scDblFinder::computeDoubletDensity(CD_data_sce, BSPARAM=IrlbaParam())
DD8_data_sce$DoubletScore <- scDblFinder::computeDoubletDensity(DD8_data_sce, BSPARAM=IrlbaParam())
Msc_BM_data_sce$DoubletScore <- scDblFinder::computeDoubletDensity(Msc_BM_data_sce, BSPARAM=IrlbaParam())

### QC by clustering - Individual datasets. NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
### RATIONaLE: Low quality nuclei may be included due to threshold effects. Deep clustering can help to reveal if there are groups of low quality nuclei that does not mix with the remaining nuclei, and thus can be removed.
## Create Seurat objects
AD10_data_seurat <- as.Seurat(AD10_data_sce, counts = "counts", data = "logcounts") #15665  5757
CB_data_seurat <- as.Seurat(CB_data_sce, counts = "counts", data = "logcounts") # 15665 10769
CD_data_seurat <- as.Seurat(CD_data_sce, counts = "counts", data = "logcounts") #15665 39709
DD8_data_seurat <- as.Seurat(DD8_data_sce, counts = "counts", data = "logcounts") # 15665  6447
Msc_BM_data_seurat <- as.Seurat(Msc_BM_data_sce, counts = "counts", data = "logcounts") # 15665  2376


AD10_data_seurat <- as.Seurat(AD10_data_sce, counts = "counts", data = "logcounts") #15601  4656
CB_data_seurat <- as.Seurat(CB_data_sce, counts = "counts", data = "logcounts") #15601  7787
CD_data_seurat <- as.Seurat(CD_data_sce, counts = "counts", data = "logcounts") #15601  3117
DD8_data_seurat <- as.Seurat(DD8_data_sce, counts = "counts", data = "logcounts") #15601  5272
Msc_BM_data_seurat <- as.Seurat(Msc_BM_data_sce, counts = "counts", data = "logcounts") #15601  2080

saveRDS(AD10_data_seurat, "AD10_data_seurat.rds")
saveRDS(CB_data_seurat, "CB_data_seurat.rds")
saveRDS(CD_data_seurat, "CD_data_seurat.rds")
saveRDS(DD8_data_seurat, "DD8_data_seurat.rds")
saveRDS(Msc_BM_data_seurat, "Msc_BM_data_seurat.rds")

saveRDS(AD10_data_sce, "AD10_data_sce.rds")
saveRDS(CB_data_sce, "CB_data_sce.rds")
saveRDS(CD_data_sce, "CD_data_sce.rds")
saveRDS(DD8_data_sce, "DD8_data_sce.rds")
saveRDS(Msc_BM_data_sce, "Msc_BM_data_sce.rds")

rm(AD10_data_seurat,CB_data_seurat,CD_data_seurat,DD8_data_seurat,Msc_BM_data_seurat,AD10_data_sce,CB_data_sce,CD_data_sce,DD8_data_sce,Msc_BM_data_sce,AD10_data,CB_data,CD_data,DD8_data,Msc_BM_data)
rm(ExprsFun,CD_low, CB_low, DD8_low, Low, Msc_BM_low, Msc_BM_data_clusters, DD8_data_clusters, CD_data_clusters, AD10_low, AD10_data_clusters)


AD10_data_sce  <- readRDS("AD10_data_sce.rds")
CB_data_sce  <- readRDS("CB_data_sce.rds")
CD_data_sce  <- readRDS("CD_data_sce.rds")
DD8_data_sce  <- readRDS("DD8_data_sce.rds")
Msc_BM_data_sce  <- readRDS("Msc_BM_data_sce.rds")

# AD10
# Iteration 1 - Use all genes
##nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell. Low nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet. High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in GPSct be a doublet (or multiplet). 
#In combination with %mitochondrial reads, removing outliers from these groups removes most doublets/dead cells/empty droplets, hence why filtering is a common pre-processing step.
AD10_data_seurat  <- readRDS("AD10_data_seurat.rds")
VariableFeatures(AD10_data_seurat) <- rownames(AD10_data_seurat)
AD10_data_seurat <- ScaleData(AD10_data_seurat)
AD10_data_seurat <- RunPCA(AD10_data_seurat)
AD10_data_seurat <- RunUMAP(AD10_data_seurat, dims=1:20, reduction="pca")
AD10_data_seurat <- FindNeighbors(object = AD10_data_seurat, dims = 1:20, reduction = "pca")
AD10_data_seurat <- FindClusters(AD10_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(AD10_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(AD10_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(AD10_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(AD10_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(AD10_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(AD10_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# Clusters 0,28,35,39,44,45 have high doublet scores or counts/genes
AD10_data_seurat <- subset(AD10_data_seurat, cells=rownames(AD10_data_seurat@meta.data[ !(AD10_data_seurat@meta.data$seurat_clusters %in% c(1,3,5,8,9,10,11,12,13,14,15,16,19,20,21,22,29,30,32,35,37,38,41,42,43,44,45,48)),]))



# Iteration 2 - Exclude ambient genes
AD10_data_seurat <- ScaleData(AD10_data_seurat)
AD10_data_seurat <- RunPCA(AD10_data_seurat)
AD10_data_seurat <- RunUMAP(AD10_data_seurat, dims=1:20, reduction="pca")
AD10_data_seurat <- FindNeighbors(object = AD10_data_seurat, dims = 1:20, reduction = "pca")
AD10_data_seurat <- FindClusters(AD10_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(AD10_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(AD10_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(AD10_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(AD10_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(AD10_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(AD10_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# Clusters 35 have high doublet scores or counts/genes
AD10_data_seurat <- subset(AD10_data_seurat, cells=rownames(AD10_data_seurat@meta.data[ !(AD10_data_seurat@meta.data$seurat_clusters %in% c(0,15,17,18,28,29,34)),]))

# Iteration 2 - Exclude ambient genes
AD10_data_seurat <- ScaleData(AD10_data_seurat)
AD10_data_seurat <- RunPCA(AD10_data_seurat)
AD10_data_seurat <- RunUMAP(AD10_data_seurat, dims=1:20, reduction="pca")
AD10_data_seurat <- FindNeighbors(object = AD10_data_seurat, dims = 1:20, reduction = "pca")
AD10_data_seurat <- FindClusters(AD10_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(AD10_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(AD10_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(AD10_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(AD10_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(AD10_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(AD10_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# No more clusters to remove

#CB4
CB_data_seurat  <- readRDS("CB_data_seurat.rds")
VariableFeatures(CB_data_seurat) <- rownames(CB_data_seurat)
CB_data_seurat <- ScaleData(CB_data_seurat)
CB_data_seurat <- RunPCA(CB_data_seurat)
CB_data_seurat <- RunUMAP(CB_data_seurat, dims=1:20, reduction="pca")
CB_data_seurat <- FindNeighbors(object = CB_data_seurat, dims = 1:20, reduction = "pca")
CB_data_seurat <- FindClusters(CB_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(CB_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(CB_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(CB_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(CB_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(CB_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(CB_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# Clusters 38 have high doublet scores or counts/genes
CB_data_seurat <- subset(CB_data_seurat, cells=rownames(CB_data_seurat@meta.data[ !(CB_data_seurat@meta.data$seurat_clusters %in% c(0,1,3,4,5,6,7,8,10,15,16,17,18,20,21,22,23,24,26,27,28,29,30,32,33,37,38,40,42,43,44,45,47,49)),]))

# Iteration 2 - Exclude ambient genes
VariableFeatures(CB_data_seurat) <- rownames(CB_data_seurat)
CB_data_seurat <- ScaleData(CB_data_seurat)
CB_data_seurat <- RunPCA(CB_data_seurat)
CB_data_seurat <- RunUMAP(CB_data_seurat, dims=1:20, reduction="pca")
CB_data_seurat <- FindNeighbors(object = CB_data_seurat, dims = 1:20, reduction = "pca")
CB_data_seurat <- FindClusters(CB_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(CB_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(CB_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(CB_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(CB_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(CB_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(CB_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# Clusters 40 have high doublet scores or counts/genes
CB_data_seurat <- subset(CB_data_seurat, cells=rownames(CB_data_seurat@meta.data[ !(CB_data_seurat@meta.data$seurat_clusters %in% c(40)),]))

# Iteration 3 - Exclude ambient genes
VariableFeatures(CB_data_seurat) <- rownames(CB_data_seurat)
CB_data_seurat <- ScaleData(CB_data_seurat)
CB_data_seurat <- RunPCA(CB_data_seurat)
CB_data_seurat <- RunUMAP(CB_data_seurat, dims=1:20, reduction="pca")
CB_data_seurat <- FindNeighbors(object = CB_data_seurat, dims = 1:20, reduction = "pca")
CB_data_seurat <- FindClusters(CB_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(CB_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(CB_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(CB_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(CB_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(CB_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(CB_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# No more clusters to remove

#CD8
CD_data_seurat  <- readRDS("CD_data_seurat.rds")
VariableFeatures(CD_data_seurat) <- rownames(CD_data_seurat)
CD_data_seurat <- ScaleData(CD_data_seurat)
CD_data_seurat <- RunPCA(CD_data_seurat)
CD_data_seurat <- RunUMAP(CD_data_seurat, dims=1:20, reduction="pca")
CD_data_seurat <- FindNeighbors(object = CD_data_seurat, dims = 1:20, reduction = "pca")
CD_data_seurat <- FindClusters(CD_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(CD_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(CD_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(CD_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(CD_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(CD_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(CD_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# Clusters 1,2,7,12,32,109 have high doublet scores or counts/genes
CD_data_seurat <- subset(CD_data_seurat, cells=rownames(CD_data_seurat@meta.data[ !(CD_data_seurat@meta.data$seurat_clusters %in% c(0,3,4,6,13,14,16,18,19,20,22,23,24,25,26,27,28,29,30,38,40,41,43)),]))


VariableFeatures(CD_data_seurat) <- rownames(CD_data_seurat)
CD_data_seurat <- ScaleData(CD_data_seurat)
CD_data_seurat <- RunPCA(CD_data_seurat)
CD_data_seurat <- RunUMAP(CD_data_seurat, dims=1:20, reduction="pca")
CD_data_seurat <- FindNeighbors(object = CD_data_seurat, dims = 1:20, reduction = "pca")
CD_data_seurat <- FindClusters(CD_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(CD_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(CD_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(CD_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(CD_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(CD_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(CD_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))

CD_data_seurat <- subset(CD_data_seurat, cells=rownames(CD_data_seurat@meta.data[ !(CD_data_seurat@meta.data$seurat_clusters %in% c(5)),]))


VariableFeatures(CD_data_seurat) <- rownames(CD_data_seurat)
CD_data_seurat <- ScaleData(CD_data_seurat)
CD_data_seurat <- RunPCA(CD_data_seurat)
CD_data_seurat <- RunUMAP(CD_data_seurat, dims=1:20, reduction="pca")
CD_data_seurat <- FindNeighbors(object = CD_data_seurat, dims = 1:20, reduction = "pca")
CD_data_seurat <- FindClusters(CD_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(CD_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(CD_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(CD_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(CD_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(CD_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(CD_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))

# DD8
DD8_data_seurat  <- readRDS("DD8_data_seurat.rds")
VariableFeatures(DD8_data_seurat) <- rownames(DD8_data_seurat)
DD8_data_seurat <- ScaleData(DD8_data_seurat)
DD8_data_seurat <- RunPCA(DD8_data_seurat)
DD8_data_seurat <- RunUMAP(DD8_data_seurat, dims=1:20, reduction="pca")
DD8_data_seurat <- FindNeighbors(object = DD8_data_seurat, dims = 1:20, reduction = "pca")
DD8_data_seurat <- FindClusters(DD8_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(DD8_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(DD8_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(DD8_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(DD8_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(DD8_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(DD8_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# Clusters 30,37,39) have high doublet scores or counts/genes
DD8_data_seurat <- subset(DD8_data_seurat, cells=rownames(DD8_data_seurat@meta.data[ !(DD8_data_seurat@meta.data$seurat_clusters %in% c(0,1,2,3,4,5,6,8,10,11,13,14,17,18,27,29,37,38,43,44,45)),]))


# Iteration 2 - Exclude clusters
VariableFeatures(DD8_data_seurat) <- rownames(DD8_data_seurat)
DD8_data_seurat <- ScaleData(DD8_data_seurat)
DD8_data_seurat <- RunPCA(DD8_data_seurat)
DD8_data_seurat <- RunUMAP(DD8_data_seurat, dims=1:20, reduction="pca")
DD8_data_seurat <- FindNeighbors(object = DD8_data_seurat, dims = 1:20, reduction = "pca")
DD8_data_seurat <- FindClusters(DD8_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(DD8_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(DD8_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(DD8_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(DD8_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(DD8_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(DD8_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))

DD8_data_seurat <- subset(DD8_data_seurat, cells=rownames(DD8_data_seurat@meta.data[ !(DD8_data_seurat@meta.data$seurat_clusters %in% c(16)),]))


# Iteration 2 - Exclude clusters
VariableFeatures(DD8_data_seurat) <- rownames(DD8_data_seurat)
DD8_data_seurat <- ScaleData(DD8_data_seurat)
DD8_data_seurat <- RunPCA(DD8_data_seurat)
DD8_data_seurat <- RunUMAP(DD8_data_seurat, dims=1:20, reduction="pca")
DD8_data_seurat <- FindNeighbors(object = DD8_data_seurat, dims = 1:20, reduction = "pca")
DD8_data_seurat <- FindClusters(DD8_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(DD8_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(DD8_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(DD8_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(DD8_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(DD8_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(DD8_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# No more clusters to remove


#Msc TERT
Msc_BM_data_seurat  <- readRDS("Msc_BM_data_seurat.rds")
VariableFeatures(Msc_BM_data_seurat) <- rownames(Msc_BM_data_seurat)
Msc_BM_data_seurat <- ScaleData(Msc_BM_data_seurat)
Msc_BM_data_seurat <- RunPCA(Msc_BM_data_seurat)
Msc_BM_data_seurat <- RunUMAP(Msc_BM_data_seurat, dims=1:20, reduction="pca")
Msc_BM_data_seurat <- FindNeighbors(object = Msc_BM_data_seurat, dims = 1:20, reduction = "pca")
Msc_BM_data_seurat <- FindClusters(Msc_BM_data_seurat, resolution = 6, algorithm = 1)
plot1 <- VlnPlot(Msc_BM_data_seurat, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(Msc_BM_data_seurat, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(Msc_BM_data_seurat, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(Msc_BM_data_seurat, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(Msc_BM_data_seurat, label=T) + NoLegend()
plot6 <- FeaturePlot(Msc_BM_data_seurat, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))
# No more clusters to remove


saveRDS(AD10_data_seurat, "AD10_data_seurat.rds")
saveRDS(CB_data_seurat, "CB_data_seurat.rds")
saveRDS(CD_data_seurat, "CD_data_seurat.rds")
saveRDS(DD8_data_seurat, "DD8_data_seurat.rds")
saveRDS(Msc_BM_data_seurat, "Msc_BM_data_seurat.rds")


### QC by clustering - Replicates NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
# Merge the datasets
Tert <- merge(AD10_data_seurat, y= c(CB_data_seurat,CD_data_seurat, DD8_data_seurat, Msc_BM_data_seurat))
Tert$Dataset <-gsub("_.*","", colnames(Tert))
dim(Tert) 
15665 60496

35523 10036 #1500 threshold detected


AD10   CB   CD  DD8  Msc
1600 2414 1482 2460 2080




# Iteration 1. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
VariableFeatures(Tert) <- rownames(Tert)
Tert <- ScaleData(Tert)
Tert <- RunPCA(Tert)
Tert <- RunHarmony(Tert, group.by.vars="Dataset", dims.use=1:20)
Tert <- RunUMAP(Tert, dims=1:20, reduction="harmony")
Tert <- FindNeighbors(object = Tert, dims = 1:20, reduction = "harmony")
Tert <- FindClusters(Tert, resolution = 5, algorithm = 1)
plot1 <- VlnPlot(Tert, c("detected"), pt.size=0) +NoLegend()
plot2 <- VlnPlot(Tert, c("sum"), pt.size=0) +NoLegend()
plot3 <- VlnPlot(Tert, c("DoubletScore"), pt.size=0) +NoLegend()
plot4 <- VlnPlot(Tert, c("subsets_Mito_percent"), pt.size=0) +NoLegend()
plot5 <- DimPlot(Tert, label=T) + NoLegend()
plot6 <- FeaturePlot(Tert, features = "PCNA")
CombinePlots(list(plot1,plot2,plot3,plot4,plot5,plot6))



DimPlot(Tert, group.by = "Dataset", pt.size = .1)


Idents(object = Tert) <- "Dataset"
DimPlot(Tert, group.by = "Dataset", pt.size = .1)
DimPlot(Tert, split.by = "Dataset", pt.size = .1)



saveRDS(AD10_data_seurat, "AD10_data_seurat.rds")
saveRDS(CB_data_seurat, "CB_data_seurat.rds")
saveRDS(DD8_data_seurat, "DD8_data_seurat.rds")
saveRDS(Msc_BM_data_seurat, "Msc_BM_data_seurat.rds")

saveRDS(Tert, "Tert_zUMI_mix.rds")


Tert  <- readRDS("Tert_zUMI_mix.rds")

### Subset each SCE object to the identified high-quality droplets. 
AD10_data_sce <- AD10_data_sce[,which(colnames(AD10_data_sce) %in% colnames(Tert))] #15665  5081
CB_data_sce <- CB_data_sce[,which(colnames(CB_data_sce) %in% colnames(Tert))] #15665 10495
CD_data_sce <- CD_data_sce[,which(colnames(CD_data_sce) %in% colnames(Tert))] #15665 36393
DD8_data_sce <- DD8_data_sce[,which(colnames(DD8_data_sce) %in% colnames(Tert))] #15665  6151
Msc_BM_data_sce <- Msc_BM_data_sce[,which(colnames(Msc_BM_data_sce) %in% colnames(Tert))] #15665  2376

AD10_data_sce <- AD10_data_sce[,which(colnames(AD10_data_sce) %in% colnames(Tert))] #15601  1600
CB_data_sce <- CB_data_sce[,which(colnames(CB_data_sce) %in% colnames(Tert))] #15601  2414
CD_data_sce <- CD_data_sce[,which(colnames(CD_data_sce) %in% colnames(Tert))] #15601  1482
DD8_data_sce <- DD8_data_sce[,which(colnames(DD8_data_sce) %in% colnames(Tert))] #15601  2460
Msc_BM_data_sce <- Msc_BM_data_sce[,which(colnames(Msc_BM_data_sce) %in% colnames(Tert))] #15601  2080


### Normalize the datasets. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
# Cluster the nuclei in each dataset
AD10_data_clusters <- quickCluster(AD10_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())
CB_data_clusters <- quickCluster(CB_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())
CD_data_clusters <- quickCluster(CD_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())
DD8_data_clusters <- quickCluster(DD8_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam())
Msc_BM_data_clusters <- quickCluster(Msc_BM_data_sce, use.ranks=GPSLSE, BSPARAM=IrlbaParam()) 


# Compute scaling GPSctors
AD10_data_sce <- computeSumGPSctors(AD10_data_sce, min.mean=0.1, cluster=AD10_data_clusters)
CB_data_sce <- computeSumGPSctors(CB_data_sce, min.mean=0.1, cluster=CB_data_clusters)
CD_data_sce <- computeSumGPSctors(CD_data_sce, min.mean=0.1, cluster=CD_data_clusters)
DD8_data_sce <- computeSumGPSctors(DD8_data_sce, min.mean=0.1, cluster=DD8_data_clusters)
Msc_BM_data_sce <- computeSumGPSctors(Msc_BM_data_sce, min.mean=0.1, cluster=Msc_BM_data_clusters) 
# Normalize the counts
AD10_data_sce <- logNormCounts(AD10_data_sce)
CB_data_sce <- logNormCounts(CB_data_sce)
CD_data_sce <- logNormCounts(CD_data_sce)
DD8_data_sce <- logNormCounts(DD8_data_sce)
Msc_BM_data_sce <- logNormCounts(Msc_BM_data_sce)

### Reduce batch effects by rescaling across datasets
# Normalize across samples
rescaled <- batchelor::multiBatchNorm(
  AD10_data_sce, 
  CB_data_sce,
  CD_data_sce,   
  DD8_data_sce,
  Msc_BM_data_sce
)

### Merge all the data and embed
# Create seurat objects
AD10_data_seurat <- as.Seurat(rescaled[[1]], counts = "counts", data = "logcounts")
CB_data_seurat <- as.Seurat(rescaled[[2]], counts = "counts", data = "logcounts")
CD_data_seurat <- as.Seurat(rescaled[[3]], counts = "counts", data = "logcounts")
DD8_data_seurat <- as.Seurat(rescaled[[4]], counts = "counts", data = "logcounts")
Msc_BM_data_seurat <- as.Seurat(rescaled[[5]], counts = "counts", data = "logcounts")

# Merge the datasets
Tert <- merge(AD10_data_seurat, y=c(CB_data_seurat,CD_data_seurat,DD8_data_seurat, Msc_BM_data_seurat))
Tert$Dataset <-gsub("_.*","", colnames(Tert))


AD10   CB   CD  DD8  Msc
1600 2414 1482 2460 2080




saveRDS(Tert, "Tert_zUMI_mix.rds")


Tert  <- readRDS("Tert_zUMI_mix.rds")
Tert <- ScaleData(Tert,features=rownames(Tert))
Tert <- FindVariableFeatures(Tert)
VariableFeatures(Tert) <- rownames(Tert)
Tert <- RunPCA(Tert, verbose = GPSLSE)
Tert <- Tert %>% RunHarmony("Dataset", plot_convergence = TRUE, do.return = TRUE)
ElbowPlot(Tert, reduction="pca")
ElbowPlot(Tert, reduction="harmony")
Tert <- Tert %>% 
    RunUMAP(reduction = "harmony", dim=1:10) %>% 
    FindNeighbors(reduction = "umap", dim= 1:2) %>% 
    FindClusters(resolution = 0.15) %>% 
    identity()

DimPlot(Tert, reduction = "umap", group.by = "RNA_snn_res.0.15",label = TRUE, pt.size = .1)
FeaturePlot(Tert, features=c("MKI67","TOP2A", "PCNA", "MCM2"))

Idents(object = Tert) <- "Dataset"	
DimPlot(Tert, group.by = "Dataset", pt.size = .1)
DimPlot(Tert, split.by = "Dataset", pt.size = .1)

Idents(object = Tert) <- "Dataset"

plotAD <- DimPlot(Tert, cells=WhichCells(Tert, idents="AD10"), group.by= "RNA_snn_res.0.15", pt.size = 0.2)+ xlim(c(-15,10)) +ylim(c(-8,16))+ NoLegend()
plotCB <- DimPlot(Tert, cells=WhichCells(Tert, idents="CB"), group.by= "RNA_snn_res.0.15", pt.size = 0.2)+ xlim(c(-15,10)) +ylim(c(-8,16))+ NoLegend()
plotCD <- DimPlot(Tert, cells=WhichCells(Tert, idents="CD"), group.by= "RNA_snn_res.0.15", pt.size = 0.2)+ xlim(c(-15,10)) +ylim(c(-8,16))+ NoLegend()
plotDD8 <- DimPlot(Tert, cells=WhichCells(Tert, idents="DD8"), group.by= "RNA_snn_res.0.15", pt.size = 0.2)+ xlim(c(-15,10)) +ylim(c(-8,16))+ NoLegend()
plotMsc <- DimPlot(Tert, cells=WhichCells(Tert, idents="Msc"), group.by= "RNA_snn_res.0.15", pt.size = 0.2)+ xlim(c(-15,10)) +ylim(c(-8,16))+ NoLegend()
DimPlot(Tert, cells=WhichCells(Tert, idents="Dataset"), group.by= "RNA_snn_res.0.15", pt.size = 0.2)+ xlim(c(-15,10)) +ylim(c(-8,16))+ NoLegend()

CombinePlots(plots = list ( plotAD,plotCB, plotCD, plotDD8, plotMsc))


plotAD <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="AD10"), "detected", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotCB <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="CB"),"detected", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotCD <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="CD"),"detected", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotDD8 <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="DD8"), "detected", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotMsc <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="Msc"), "detected", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotTert <- FeaturePlot(Tert, cells=WhichCells(Tert), "detected", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))

CombinePlots(plots = list ( plotAD,plotCB, plotCD, plotDD8, plotMsc))

plotAD <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="AD10"), "sum", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotCB <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="CB"),"sum", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotCD <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="CD"),"sum", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotDD8 <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="DD8"), "sum", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
plotMsc <- FeaturePlot(Tert, cells=WhichCells(Tert, idents="Msc"), "sum", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
FeaturePlot(Tert, cells=WhichCells(Tert), "sum", pt.size = 0.2, cols = c("pink", "black"))+ xlim(c(-15,10)) +ylim(c(-8,16))
CombinePlots(plots = list ( plotAD,plotCB, plotCD, plotDD8, plotMsc))


DoHeatmap(Tert, features=as.character(Cluster_markers[Cluster_markers$Marker=="Exclusive",'Gene']))
Idents(object = Tert) <- "RNA_snn_res.0.15"
Clusters  <- sort(levels(Tert$RNA_snn_res.0.15))
Result <- list()
# Defining enriched genes, e.g. one cluster versus all other data points
Tert_markers <- FindAllMarkers(Tert, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
Tert_markers <- Tert_markers[Tert_markers$p_val_adj < 0.05,]

# Defining marker genes, e.g. one cluster versus all other clusters individually and clean up results in a list
cluster_name <- Clusters
for(clust in Clusters){
	Pairs <- data.frame(GroupA = clust, GroupB = Clusters[!(Clusters %in% clust)])
	Exclusive <- c()
	for (i in 1:nrow(Pairs)){
		tmp <- FindMarkers(Tert, ident.1 = clust, ident.2 = Pairs[i,2], min.pct = 0.1, only.pos = TRUE)
		tmp <- tmp[tmp$p_val_adj < 0.05,]
		Exclusive <- c(Exclusive,as.character(rownames(tmp)))
		print(paste("Cluster",clust,"out of",length(Clusters),"versus Cluster",Pairs[i,2],sep=" "))
	}
	tmp <- data.frame(table(Exclusive))
	tmp1 <- tmp[tmp$Freq==length(Clusters)-1,]
	tmp2 <- tmp[tmp$Freq > 0,]
	if(length(Exclusive) >= 0 & nrow(tmp1)>0 & length(Tert_markers[Tert_markers$cluster == clust,'gene'])>0){
			Result[[clust]] <- rbind(
			data.frame("Gene"= Tert_markers[Tert_markers$cluster == clust,'gene'],"Marker"="Enriched"),
			data.frame("Gene"= tmp1[,'Exclusive'],"Marker"="Exclusive"),
			data.frame("Gene"= tmp2[,'Exclusive'],"Marker"="Significant"))
	} else if(length(Exclusive) > 0 & length(Tert_markers[Tert_markers$cluster == clust,'gene'])>0){
			Result[[clust]] <- rbind(
			data.frame("Gene"= Tert_markers[Tert_markers$cluster == clust,'gene'],"Marker"="Enriched"),
			data.frame("Gene"= tmp2[,'Exclusive'],"Marker"="Significant"))
	} else if(length(Tert_markers[Tert_markers$cluster == clust,'gene'])>=0){
			Result[[clust]] <- rbind(
			data.frame("Gene"= Tert_markers[Tert_markers$cluster == clust,'gene'],"Marker"="Enriched"))
	} else if(length(Exclusive) > 0){
			Result[[clust]] <- rbind(
			data.frame("Gene"= tmp2[,'Exclusive'],"Marker"="Significant"))
	} else {
			cluster_name <- cluster_name[!cluster_name %in% clust]
	}

}
names(Result) <- paste("Cl_",cluster_name,sep="")
Cluster_markers <- bind_rows(Result, .id = "column_label")

saveRDS(Result, "Result_nonprolif.rds")
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
    Exclusive Enriched Significant
0          0      235        1158
1          7      780        1457
10        56      815        1719
2         14      674        1544
3        174     1065        1722
4         43      889        1387
5          1      708        1475
6         14     1539        1811
7         47     1974        1847
8          6      537        1428
9         20     1620        1917

DoHeatmap(Tert, features=as.character(Cluster_markers[Cluster_markers$Marker=="Exclusive",'Gene']))

## Pathway analysis. NOTE: As WikiPathways is updated results may change!
# Prepare for WikiPathways analysis
gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", format = "gmt")
wp2gene <- clusterProfiler::read.gmt("wikipathways-20210210-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name", "version", "wpid", "org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) # TERM2NAME
Pathways <- list()

library(rWikiPathways)
library(org.Hs.eg.db)
# Loop through all clusters

for (i in 1:length(Result)) {
  # Extract list
  Tmp <- Result[[i]]
  
  # Convert gene names
  Enriched_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker == "Enriched",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  Exclusive_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker == "Exclusive",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  Significant_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker == "Significant",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Do the pathway analysis
  if (sum(wp2gene$gene %in% Enriched_Entrez[[2]]) > 0) { wiki_enriched <- clusterProfiler::enricher(Enriched_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_enriched <- data.frame() }
  if (sum(wp2gene$gene %in% Exclusive_Entrez[[2]]) > 0) { wiki_exclusive <- clusterProfiler::enricher(Exclusive_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_exclusive <- data.frame() }
  if (sum(wp2gene$gene %in% Significant_Entrez[[2]]) > 0) { wiki_significant <- clusterProfiler::enricher(Significant_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_significant <- data.frame() }

  # Adding gene symbols to the resulting pathway file
  if (nrow(wiki_enriched) > 0) { wiki_enriched <- as.data.frame(DOSE::setReadable(wiki_enriched, org.Hs.eg.db, keyType = "ENTREZID")) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive <- as.data.frame(DOSE::setReadable(wiki_exclusive, org.Hs.eg.db, keyType = "ENTREZID")) }
  if (nrow(wiki_significant) > 0) { wiki_significant <- as.data.frame(DOSE::setReadable(wiki_significant, org.Hs.eg.db, keyType = "ENTREZID")) }

  # Calculate gene ratios (i.e. the fraction of genes in the pathway)
  if (nrow(wiki_enriched) > 0) { wiki_enriched$Enriched_GeneRatio <- as.numeric(substr(wiki_enriched$GeneRatio, 1, regexpr("/", wiki_enriched$GeneRatio)-1)) / as.numeric(substr(wiki_enriched$GeneRatio, regexpr("/", wiki_enriched$GeneRatio)+1, nchar(wiki_enriched$GeneRatio))) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive$Exclusive_GeneRatio <- as.numeric(substr(wiki_exclusive$GeneRatio, 1, regexpr("/", wiki_exclusive$GeneRatio)-1)) / as.numeric(substr(wiki_exclusive$GeneRatio, regexpr("/", wiki_exclusive$GeneRatio)+1, nchar(wiki_exclusive$GeneRatio))) }
  if (nrow(wiki_significant) > 0) { wiki_significant$Exclusive_GeneRatio <- as.numeric(substr(wiki_significant$GeneRatio, 1, regexpr("/", wiki_significant$GeneRatio)-1)) / as.numeric(substr(wiki_significant$GeneRatio, regexpr("/", wiki_significant$GeneRatio)+1, nchar(wiki_significant$GeneRatio))) }

  # Calculate pathway ratios (i.e. the fraction of the pathway in the genes)
  if (nrow(wiki_enriched) > 0) { wiki_enriched$Enriched_PathRatio <- as.numeric(substr(wiki_enriched$GeneRatio, 1, regexpr("/", wiki_enriched$GeneRatio)-1)) / as.numeric(substr(wiki_enriched$BgRatio, 1, regexpr("/", wiki_enriched$BgRatio)-1)) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive$Exclusive_PathRatio <- as.numeric(substr(wiki_exclusive$GeneRatio, 1, regexpr("/", wiki_exclusive$GeneRatio)-1)) / as.numeric(substr(wiki_exclusive$BgRatio, 1, regexpr("/", wiki_exclusive$BgRatio)-1)) }
  if (nrow(wiki_significant) > 0) { wiki_significant$Exclusive_PathRatio <- as.numeric(substr(wiki_significant$GeneRatio, 1, regexpr("/", wiki_significant$GeneRatio)-1)) / as.numeric(substr(wiki_significant$BgRatio, 1, regexpr("/", wiki_significant$BgRatio)-1)) }

  # Set column names
  if (nrow(wiki_enriched) > 0) { colnames(wiki_enriched)[c(5,7,8)] <- c("Enriched_Pvalue","Enriched_FDR","Enriched_Genes") }
  if (nrow(wiki_exclusive) > 0) { colnames(wiki_exclusive)[c(5,7,8)] <- c("Exclusive_Pvalue","Exclusive_FDR","Exclusive_Genes") }
  if (nrow(wiki_significant) > 0) { colnames(wiki_significant)[c(5,7,8)] <- c("Significant_Pvalue","Significant_FDR","Significant_Genes") }

  # Combine the results
  significant <- c(wiki_significant[ wiki_significant$Significant_FDR <= 0.05,1],wiki_enriched[ wiki_enriched$Enriched_FDR <= 0.05,1],  wiki_exclusive[ wiki_exclusive$Exclusive_FDR <= 0.05,1])
  if (length(significant) > 0) {
    # Subset from all pathways to only significant ones
    all <- wpid2name[ wpid2name[,1] %in% significant,]
    all <- all[ duplicated(all[,1])==F,]
    colnames(all) <- c("ID", "Description")
    
    # Merge statistics
    if (nrow(wiki_significant) > 0) { all <- merge(all, wiki_significant[,c(1,10,11,5,7,8)], all.x=T, by="ID") }
	if (nrow(wiki_enriched) > 0) { all <- merge(all, wiki_enriched[,c(1,10,11,5,7,8)], all.x=T, by="ID") }
    if (nrow(wiki_exclusive) > 0) { all <- merge(all, wiki_exclusive[,c(1,10,11,5,7,8)], all.x=T, by="ID") }

    # Handle NAs appropriately (ratios to 0, pvalues/FDRs to 1 and gene lists to "")
    for (m in grep("Ratio", colnames(all))) { all[ is.na(all[,m]),m] <- 0 }
    for (m in grep("Pvalue", colnames(all))) { all[ is.na(all[,m]),m] <- 1 }
    for (m in grep("FDR", colnames(all))) { all[ is.na(all[,m]),m] <- 1 }
    for (m in grep("Genes", colnames(all))) { all[ is.na(all[,m]),m] <- "" }
    
    # Store the results
    Pathways[[length(Pathways)+1]] <- all # The list is the basis for the Supplementary Tables
	names(Pathways)[length(Pathways)] <- paste('Cl_',Clusters[i],sep="")
  }
}

for (i in 1:length(Pathways)){
	write.table(Pathways[[i]], file=paste("Wiki_",names(Pathways)[i],".txt",sep=""),quote=F, col.names=T, row.names=F,sep="\t")
}




Tert$Subtype <- "NA"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 0, "Subtype"] <- "PI3K/AKT/mTOR - VitD3 signaling"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 1, "Subtype"] <- "Gluthathione metabolism"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 2, "Subtype"] <- "Ras signaling"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 3, "Subtype"] <- "Proliferative"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 4, "Subtype"] <- "Purine metabolism"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 5, "Subtype"] <- "Calcium reguation"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 6, "Subtype"] <- "Parkin−Ubiquitin Proteasomal System"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 7, "Subtype"] <- "Oxidative phosphorylation"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 8, "Subtype"] <- "TGF-b / BMP signalling"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 9, "Subtype"] <- "Translation factors"
Tert@meta.data[Tert@meta.data$RNA_snn_res.0.15 %in% 10, "Subtype"] <- "Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds"

saveRDS(Tert, "Tert_Subtypes.rds")



##find specific markers for each cell line and cluster 
AD10 <- subset(Tert, subset = Dataset == c('AD10'))
AD10.markers <- FindAllMarkers(AD10,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DD8 <- subset(Tert, subset = Dataset == c('DD8'))
DD8.markers <- FindAllMarkers(DD8,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CB4 <- subset(Tert, subset = Dataset == c('CB'))
CB4.markers <- FindAllMarkers(CB4,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CD8 <- subset(Tert, subset = Dataset == c('CD'))
CD8.markers <- FindAllMarkers(CD8,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Msc <- subset(Tert, subset = Dataset == c('Msc'))
Msc.markers <- FindAllMarkers(Msc,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(AD10.markers, "AD10.markers.rds")
saveRDS(DD8.markers, "DD8.markers.rds")
saveRDS(CB4.markers, "CB4.markers.rds")
saveRDS(CD8.markers, "CD8.markers.rds")
saveRDS(Msc.markers, "Msc.markers.rds")

AD10.markers  <- readRDS("AD10.markers.rds")
DD8.markers  <- readRDS("DD8.markers.rds")
CB4.markers  <- readRDS("CB4.markers.rds")
CD8.markers  <- readRDS("CD8.markers.rds")
Msc.markers  <- readRDS("Msc.markers.rds")

Msc.markers <- Msc.markers[Msc.markers$p_val_adj < 0.05,]
AD10.markers <- AD10.markers[AD10.markers$p_val_adj < 0.05,]
DD8.markers <- DD8.markers[DD8.markers$p_val_adj < 0.05,]
CB4.markers <- CB4.markers[CB4.markers$p_val_adj < 0.05,]
CD8.markers <- CD8.markers[CD8.markers$p_val_adj < 0.05,]



AD10_PI3K <- rownames(AD10.markers[AD10.markers$cluster=="PI3K/AKT/mTOR - VitD3 signaling",])
AD10_Glutathione <- rownames(AD10.markers[AD10.markers$cluster=="Glutathione metabolism",])
AD10_Ras <- rownames(AD10.markers[AD10.markers$cluster=="Ras signaling",])
AD10_Proliferative <- rownames(AD10.markers[AD10.markers$cluster=="Proliferative",])
AD10_Purine <- rownames(AD10.markers[AD10.markers$cluster=="Purine metabolism",])
AD10_Calcium <- rownames(AD10.markers[AD10.markers$cluster=="Calcium reguation",])
AD10_PUPS <- rownames(AD10.markers[AD10.markers$cluster=="Parkin−Ubiquitin Proteasomal System",])
AD10_OP <- rownames(AD10.markers[AD10.markers$cluster=="Oxidative phosphorylation",])
AD10_BMP <- rownames(AD10.markers[AD10.markers$cluster=="TGF-b / BMP signalling",])
AD10_TF <- rownames(AD10.markers[AD10.markers$cluster=="Translation GPSctors",])
AD10_WNT <- rownames(AD10.markers[AD10.markers$cluster=="Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds",])

DD8_PI3K <- rownames(DD8.markers[DD8.markers$cluster=="PI3K/AKT/mTOR - VitD3 signaling",])
DD8_Glutathione <- rownames(DD8.markers[DD8.markers$cluster=="Glutathione metabolism",])
DD8_Ras <- rownames(DD8.markers[DD8.markers$cluster=="Ras signaling",])
DD8_Proliferative <- rownames(DD8.markers[DD8.markers$cluster=="Proliferative",])
DD8_Purine <- rownames(DD8.markers[DD8.markers$cluster=="Purine metabolism",])
DD8_Calcium <- rownames(DD8.markers[DD8.markers$cluster=="Calcium reguation",])
DD8_PUPS <- rownames(DD8.markers[DD8.markers$cluster=="Parkin−Ubiquitin Proteasomal System",])
DD8_OP <- rownames(DD8.markers[DD8.markers$cluster=="Oxidative phosphorylation",])
DD8_BMP <- rownames(DD8.markers[DD8.markers$cluster=="TGF-b / BMP signalling",])
DD8_TF <- rownames(DD8.markers[DD8.markers$cluster=="Translation GPSctors",])
DD8_WNT <- rownames(DD8.markers[DD8.markers$cluster=="Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds",])

CB4_PI3K <- rownames(CB4.markers[CB4.markers$cluster=="PI3K/AKT/mTOR - VitD3 signaling",])
CB4_Glutathione <- rownames(CB4.markers[CB4.markers$cluster=="Glutathione metabolism",])
CB4_Ras <- rownames(CB4.markers[CB4.markers$cluster=="Ras signaling",])
CB4_Proliferative <- rownames(CB4.markers[CB4.markers$cluster=="Proliferative",])
CB4_Purine <- rownames(CB4.markers[CB4.markers$cluster=="Purine metabolism",])
CB4_Calcium <- rownames(CB4.markers[CB4.markers$cluster=="Calcium reguation",])
CB4_PUPS <- rownames(CB4.markers[CB4.markers$cluster=="Parkin−Ubiquitin Proteasomal System",])
CB4_OP <- rownames(CB4.markers[CB4.markers$cluster=="Oxidative phosphorylation",])
CB4_BMP <- rownames(CB4.markers[CB4.markers$cluster=="TGF-b / BMP signalling",])
CB4_TF <- rownames(CB4.markers[CB4.markers$cluster=="Translation GPSctors",])
CB4_WNT <- rownames(CB4.markers[CB4.markers$cluster=="Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds",])

CD8_PI3K <- rownames(CD8.markers[CD8.markers$cluster=="PI3K/AKT/mTOR - VitD3 signaling",])
CD8_Glutathione <- rownames(CD8.markers[CD8.markers$cluster=="Glutathione metabolism",])
CD8_Ras <- rownames(CD8.markers[CD8.markers$cluster=="Ras signaling",])
CD8_Proliferative <- rownames(CD8.markers[CD8.markers$cluster=="Proliferative",])
CD8_Purine <- rownames(CD8.markers[CD8.markers$cluster=="Purine metabolism",])
CD8_Calcium <- rownames(CD8.markers[CD8.markers$cluster=="Calcium reguation",])
CD8_PUPS <- rownames(CD8.markers[CD8.markers$cluster=="Parkin−Ubiquitin Proteasomal System",])
CD8_OP <- rownames(CD8.markers[CD8.markers$cluster=="Oxidative phosphorylation",])
CD8_BMP <- rownames(CD8.markers[CD8.markers$cluster=="TGF-b / BMP signalling",])
CD8_TF <- rownames(CD8.markers[CD8.markers$cluster=="Translation GPSctors",])
CD8_WNT <- rownames(CD8.markers[CD8.markers$cluster=="Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds",])

Msc_PI3K <- rownames(Msc.markers[Msc.markers$cluster=="PI3K/AKT/mTOR - VitD3 signaling",])
Msc_Glutathione <- rownames(Msc.markers[Msc.markers$cluster=="Gluthathione metabolism",])
Msc_Ras <- rownames(Msc.markers[Msc.markers$cluster=="Ras signaling",])
Msc_Proliferative <- rownames(Msc.markers[Msc.markers$cluster=="Proliferative",])
Msc_Purine <- rownames(Msc.markers[Msc.markers$cluster=="Purine metabolism",])
Msc_Calcium <- rownames(Msc.markers[Msc.markers$cluster=="Calcium reguation",])
Msc_PUPS <- rownames(Msc.markers[Msc.markers$cluster=="Parkin−Ubiquitin Proteasomal System",])
Msc_OP <- rownames(Msc.markers[Msc.markers$cluster=="Oxidative phosphorylation",])
Msc_BMP <- rownames(Msc.markers[Msc.markers$cluster=="TGF-b / BMP signalling",])
Msc_TF <- rownames(Msc.markers[Msc.markers$cluster=="Translation GPSctors",])
Msc_WNT <- rownames(Msc.markers[Msc.markers$cluster=="Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds",])



AD10.markers  <- readRDS("AD10.markers.rds")
DD8.markers  <- readRDS("DD8.markers.rds")
CB4.markers  <- readRDS("CB4.markers.rds")
CD8.markers  <- readRDS("CD8.markers.rds")
Msc.markers  <- readRDS("Msc.markers.rds")

x = list(AD10_WNT,
         DD8_WNT,
		 CB4_WNT,
		 CD8_WNT,
		 Msc_WNT)

# Flatten the list into a single vector
all_genes <- unlist(x)

# Count occurrences of each gene
gene_counts <- table(all_genes)

# Convert to a data frame for better visualization
gene_count_df <- as.data.frame(gene_counts)

# Rename columns
colnames(gene_count_df) <- c("Gene", "Count")

# Print result
print(gene_count_df)

gene_count_df <- gene_count_df[order(-gene_count_df$Count),]

write.table(gene_count_df, "Counts_PI3K_AKT_mTOR.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_Glutathione metabolism.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_Ras signaling.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_Proliferative.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_Purine metabolism.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_Calcium reguation.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_Parkin-Ubiquitin.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_OP.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_TGF-b_BMP signalling.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_TF.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(gene_count_df, "Counts_WNT.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(All_markers_final, "All_markers_final.txt", quote=F, sep="\t", col.names=T, row.names=F)

write.table(AD10.markers, "AD10.markers.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(DD8.markers, "DD8.markers.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(CB4.markers, "CB4.markers.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(CD8.markers, "CD8.markers.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(Msc.markers, "Msc.markers.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(Tert.markers, "Tert.markers.txt", quote=F, sep="\t", col.names=T, row.names=F)

##Raw data are generated from the in vitro assay by measuring cell viability and ALP activity 
###Load the following objects from https://osf.io/wxpgn/
##Celastrol treatment
AD10_viability <- read.xlsx("AD10_Celastrol_viability.xlsx", 
                       sheetIndex = 1, header=TRUE)
AD10_celastrol <- read.xlsx("AD10_Celastrol.xlsx", 
                       sheetIndex = 1, header=TRUE)
DD8_viability <- read.xlsx("DD8_Celastrol_viability.xlsx", 
                        sheetIndex = 1, header=TRUE)
DD8_celastrol <- read.xlsx("DD8_Celastrol.xlsx", 
                           sheetIndex = 1, header=TRUE)
CB4_viability <- read.xlsx("CB4_Celastrol_viability.xlsx", 
                           sheetIndex = 1, header=TRUE)
CB4_celastrol <- read.xlsx("CB4_Celastrol.xlsx", 
                           sheetIndex = 1, header=TRUE)
CD8_viability <- read.xlsx("CD8_Celastrol_viability.xlsx", 
                           sheetIndex = 1, header=TRUE)
CD8_celastrol <- read.xlsx("CD8_Celastrol.xlsx", 
                           sheetIndex = 1, header=TRUE)

AD10_Celastrol_ALP <- AD10_celastrol / AD10_viability
AD10_Celastrol_ALP_norm<- AD10_Celastrol_ALP[, c("DMSO", "X0.03uM")] / AD10_Celastrol_ALP$UT

DD8_Celastrol_ALP <- DD8_celastrol / DD8_viability
DD8_Celastrol_ALP_norm<- DD8_Celastrol_ALP[, c("DMSO", "X0.03uM")] / DD8_Celastrol_ALP$UT

CB4_Celastrol_ALP <- CB4_celastrol / CB4_viability
CB4_Celastrol_ALP_norm<- CB4_Celastrol_ALP[, c("DMSO", "X0.03uM")] / CB4_Celastrol_ALP$UT

CD8_Celastrol_ALP <- CD8_celastrol / CD8_viability
CD8_Celastrol_ALP_norm<- CD8_Celastrol_ALP[, c("DMSO", "X0.03uM")] / CD8_Celastrol_ALP$UT

library(dplyr)
library(tidyr)

# Add an ID to each dataframe
AD10_Celastrol_ALP_norm$CellLine <- "AD10"
DD8_Celastrol_ALP_norm$CellLine <- "DD8"
CB4_Celastrol_ALP_norm$CellLine <- "CB4"
CD8_Celastrol_ALP_norm$CellLine <- "CD8"


Celastrol_all <- bind_rows(AD10_Celastrol_ALP_norm, DD8_Celastrol_ALP_norm, CB4_Celastrol_ALP_norm, CD8_Celastrol_ALP_norm)


Celastrol <- Celastrol_all %>%
  pivot_longer(
    cols = -CellLine,
    names_to = "Treatment",
    values_to = "Data"
  )

Celastrol <- Celastrol %>%
  mutate(Conditions = paste(CellLine, Treatment, sep = "_"))

Celastrol <- Celastrol %>%
  select(Conditions, Data)


Celastrol$Conditions <- gsub("X0\\.03uM", "0.03uM", Celastrol$Conditions)

UT_df <- data.frame(
  Conditions = "UT",
  Data = rep(1, 4)  # number of replicates
)

Celastrol <- bind_rows(UT_df, Celastrol)
Celastrol$Data <- as.numeric(gsub(",", ".", Celastrol$Data))
library(writexl)
write_xlsx(Celastrol, "C:\\Users\\acaci\\Desktop\\Clones\\Celastrol.xlsx")

##Truli treatment 
AD10_viability <- read.xlsx("AD10_Truli_viability.xlsx", 
                       sheetIndex = 1, header=TRUE)
AD10_Truli <- read.xlsx("AD10_Truli.xlsx", 
                       sheetIndex = 1, header=TRUE)
DD8_viability <- read.xlsx("DD8_Truli_viability.xlsx", 
                        sheetIndex = 1, header=TRUE)
DD8_Truli <- read.xlsx("DD8_Truli.xlsx", 
                           sheetIndex = 1, header=TRUE)
CB4_viability <- read.xlsx("CB4_Truli_viability.xlsx", 
                           sheetIndex = 1, header=TRUE)
CB4_Truli <- read.xlsx("CB4_Truli.xlsx", 
                           sheetIndex = 1, header=TRUE)
CD8_viability <- read.xlsx("CD8_Truli_viability.xlsx", 
                           sheetIndex = 1, header=TRUE)
CD8_Truli <- read.xlsx("CD8_Truli.xlsx", 
                           sheetIndex = 1, header=TRUE)

AD10_Truli_ALP <- AD10_Truli / AD10_viability
AD10_Truli_ALP_norm<- AD10_Truli_ALP[, c("DMSO", "X0.125uM")] / AD10_Truli_ALP$UT

DD8_Truli_ALP <- DD8_Truli / DD8_viability
DD8_Truli_ALP_norm<- DD8_Truli_ALP[, c("DMSO", "X0.125uM")] / DD8_Truli_ALP$UT

CB4_Truli_ALP <- CB4_Truli / CB4_viability
CB4_Truli_ALP_norm<- CB4_Truli_ALP[, c("DMSO", "X0.125uM")] / CB4_Truli_ALP$UT

CD8_Truli_ALP <- CD8_Truli / CD8_viability
CD8_Truli_ALP_norm<- CD8_Truli_ALP[, c("DMSO", "X0.125uM")] / CD8_Truli_ALP$UT

library(dplyr)
library(tidyr)

# Add an ID to each dataframe
AD10_Truli_ALP_norm$CellLine <- "AD10"
DD8_Truli_ALP_norm$CellLine <- "DD8"
CB4_Truli_ALP_norm$CellLine <- "CB4"
CD8_Truli_ALP_norm$CellLine <- "CD8"


Truli_all <- bind_rows(AD10_Truli_ALP_norm, DD8_Truli_ALP_norm, CB4_Truli_ALP_norm, CD8_Truli_ALP_norm)


Truli <- Truli_all %>%
  pivot_longer(
    cols = -CellLine,
    names_to = "Treatment",
    values_to = "Data"
  )

Truli <- Truli %>%
  mutate(Conditions = paste(CellLine, Treatment, sep = "_"))

Truli <- Truli %>%
  select(Conditions, Data)


Truli$Conditions <- gsub("X0\\.125uM", "0.125uM", Truli$Conditions)

UT_df <- data.frame(
  Conditions = "UT",
  Data = rep(1, 4)  # number of replicates
)

Truli <- bind_rows(UT_df, Truli)
Truli$Data <- as.numeric(gsub(",", ".", Truli$Data))
library(writexl)
write_xlsx(Truli, "C:\\Users\\acaci\\Desktop\\Clones\\Truli.xlsx")
