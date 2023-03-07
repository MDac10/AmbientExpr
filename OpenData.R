##Megan Da Costa
#R package for ambient RNA removal

##pathToDataSet should be to level 'Remapped_for_ParaAmb\<Data name>\cellranger_exons_only\filtered_feature_bc_matrix'
#makes a sparse matrix readMM
OpenSplicedData <- function(pathToDataSet){
  setwd(pathToDataSet)
  filtered.mat <- Matrix::readMM(file = paste(pathToDataSet, '\\matrix.mtx.gz', sep = ''))
  row.names(filtered.mat) <- read.table(paste(pathToDataSet, '\\features.tsv.gz', sep = '')) [,2] #Only uses the second column which holds the features of interest
  colnames(filtered.mat) <- read.table(paste(pathToDataSet, '\\barcodes.tsv.gz', sep = '')) [,1] #Only uses the first column which holds the barcodes of interest

  return(filtered.mat)
}

##pathToDataSet should be to level 'Remapped_for_ParaAmb\<Data name>\cellranger_exons_only\filtered_feature_bc_matrix'
#makes a sparse matrix readMM
OpenUnsplicedData <- function(pathToDataSet){
  setwd(pathToDataSet)
  raw.mat <- Matrix::readMM(file = paste(pathToDataSet, '\\matrix.mtx.gz', sep = ''))
  row.names(raw.mat) <- read.table(paste(pathToDataSet, '\\features.tsv.gz', sep = '')) [,2] #Only uses the second column which holds the features of interest
  colnames(raw.mat) <- read.table(paste(pathToDataSet, '\\barcodes.tsv.gz', sep = '')) [,1] #Only uses the first column which holds the barcodes of interest

  return(raw.mat)
}
