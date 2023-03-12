##Megan Da Costa
##Supervisor: Tallulah Andrews

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

###Function to create a seurat object that compares all the RNA data from this dataset
##@param: inputFile - - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@return: seur_obj - A seurat object containing all metadata for a set of RNA
MakeSeuratObject <- function(inputFile){
  all_expression <- inputFile$spliced_cells_matrix+inputFile$unspliced_cells_matrix
  all_expr <- all_expression[all_expression > 0]
  seur_obj <- CreateSeuratObject(counts = all_expression, min.cells = 3, min.features = 300)
  seur_obj <- NormalizeData(seur_obj)
  seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seur_obj)
  seur_obj <- ScaleData(seur_obj, features = all.genes)
  seur_obj  <- RunPCA(seur_obj , features = VariableFeatures(object = seur_obj))
  seur_obj <- FindNeighbors(seur_obj, dims = 1:10)
  seur_obj <- FindClusters(seur_obj, resolution = 0.5)
  return(seur_obj)
}


###Function to perform seurat comparison on all cell matrices to find marker genes, noting if the genes returned in SplicedProcessing occur within them
##@param: seuratObj - A seurat object containing all metadata for a set of RNA
##@return: markers - A dataframe containing all information about genes recognized as markers from all spliced and unspliced genes
FindMarkerGenes <- function(seuratObj){
  markers <- FindAllMarkers(seuratObj, only.pos = TRUE)
  return(markers)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
