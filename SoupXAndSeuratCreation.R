install.packages('RTools')

install.packages('SoupX')
install.packages('Matrix')

library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(DoubletFinder)
library(knitr)

setwd('C:\\Users\\megda\\Thesis_4460Z\\UnSplicedRNARemoval')

source('R\\OpenData.R')

filtered.mat <- OpenSplicedData('C:\\Users\\megda\\Thesis_4460Z\\Remapped_for_ParaAmb\\P1TLH\\cellranger_exons_only\\filtered_feature_bc_matrix')

raw.mat <- OpenUnsplicedData('C:\\Users\\megda\\Thesis_4460Z\\Remapped_for_ParaAmb\\P1TLH\\cellranger_exons_only\\raw_feature_bc_matrix')

#Creates a Seurat
seuratObj <- CreateSeuratObject(counts = filtered.mat)

seuratObj <- seuratObj[!duplicated(rownames(seuratObj)),] #Removal of all duplicate rows found within the seurat Object. These duplicate rows do not allow the SoupX Channel to be formed.

soupx.chan <- SoupChannel(tod=raw.mat, toc=filtered.mat, calcSoupProfile=FALSE)

seuratObj <- SCTransform(seuratObj, verbose = F)
seuratObj <- RunPCA(seuratObj, verbose = F)
seuratObj <- RunUMAP(seuratObj, dims = 1:30, verbose = F)
seuratObj <- FindNeighbors(seuratObj, dims = 1:30, verbose = F)
seuratObj <- FindClusters(seuratObj, verbose = T)

meta <- seuratObj@meta.data
umap <- seuratObj@reductions$umap@cell.embeddings
soupx.chan <- setClusters(soupx.chan, setNames(meta$seurat_clusters, rownames(meta)))
soupx.chan <- setDR(soupx.chan, umap)
head(meta)

soupx.chan <- autoEstCont(soupx.chan)
head(soupx.chan$soupProfile[order(soupx.chan$soupProfile$est, decreasing = T), ], n = 20)

require(Seurat)
seur_obj <- CreateSeuratObject(counts = all_expression, min.cells = 3, min.features = 300)
seur_obj <- NormalizeData(seur_obj)
seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seur_obj)
seur_obj <- ScaleData(seur_obj, features = all.genes)
seur_obj  <- RunPCA(seur_obj , features = VariableFeatures(object = seur_obj))
seur_obj <- FindNeighbors(seur_obj, dims = 1:10)
seur_obj <- FindClusters(seur_obj, resolution = 0.5)
markers <- FindAllMarkers(seur_obj, only.pos = TRUE)
