##Megan Da Costa
##Supervisor: Tallulah Andrews
#R package for ambient RNA removal

#-------------------------------------------------------------------------------------------------------------------------------------------------------------

###Function for opening filtered (spliced) data that is broken into matrix, features, and barcodes; creates a matrix
##@param: pathToDataSet - should be to level 'Remapped_for_ParaAmb\<Data name>\cellranger_exons_only\filtered_feature_bc_matrix'
##@return: filtered.mat - The UMI count matrix for a given set of data
OpenSplicedData <- function(pathToDataSet){
  setwd(pathToDataSet)
  filtered.mat <- Matrix::readMM(file = paste(pathToDataSet, '\\matrix.mtx.gz', sep = '')) #makes a sparse matrix readMM
  row.names(filtered.mat) <- read.table(paste(pathToDataSet, '\\features.tsv.gz', sep = '')) [,2] #Only uses the second column which holds the features of interest
  colnames(filtered.mat) <- read.table(paste(pathToDataSet, '\\barcodes.tsv.gz', sep = '')) [,1] #Only uses the first column which holds the barcodes of interest
  return(filtered.mat)
}

###Function for opening raw (unspliced) data that is broken into matrix, features, and barcodes; creates a matrix
##@param: pathToDataSet - should be to level 'Remapped_for_ParaAmb\<Data name>\cellranger_exons_only\filtered_feature_bc_matrix'
##@return: raw.mat - The UMI count matrix for a given set of data
OpenUnsplicedData <- function(pathToDataSet){
  setwd(pathToDataSet)
  raw.mat <- Matrix::readMM(file = paste(pathToDataSet, '\\matrix.mtx.gz', sep = '')) #makes a sparse matrix readMM
  row.names(raw.mat) <- read.table(paste(pathToDataSet, '\\features.tsv.gz', sep = '')) [,2] #Only uses the second column which holds the features of interest
  colnames(raw.mat) <- read.table(paste(pathToDataSet, '\\barcodes.tsv.gz', sep = '')) [,1] #Only uses the first column which holds the barcodes of interest

  return(raw.mat)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
