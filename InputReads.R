###Megan Da Costa###
##Open and read in input cell/gene data for comparison of spliced to unspliced ambient RNA

install.packages("Matrix")
install.packages("proxy")
install.packages("mclust")
install.packages("ggplot2")

library(Matrix)
library(proxy)
library(mclust)
library(ggplot2)

require(Seurat)

setwd('C:\\Users\\megda\\Thesis_4460Z\\ExampleInput')

input <- readRDS('KorEar_Cells1_ForMegan_inputtrimmed.rds')



###Function to organize ambient RNA samples for a given spliced ambient RNA
##@param: inputFile - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile - the attachedment after inputFile$, numeric 1 or 2 representing data ambient1_spliced or ambient2_spliced, respectively.
##@return: Ordered_Data_d_as - an ordered vector of the ambient gene expressions; This vector is multiplied by 10000 so that the data can be compared easily
Vector_ambient_spliced <- function( inputFile, ambientFile ) {
  if ( ambientFile == 1 ) {
    Ordered_Data_d_as <- sort(inputFile$ambient1_spliced, decreasing = FALSE )
  } else if ( ambientFile == 2 ) {
    Ordered_Data_d_as <- sort( inputFile$ambient2_spliced, decreasing = FALSE )
  } else {
    stop( "ambientFile parameter must be numeric 1 or 2." )
  }
  return( Ordered_Data_d_as * 10000 )
}

###Function to organize ambient RNA samples for a given unspliced ambient RNA
##@param: inputFile - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile - the attachedment after inputFile$, numeric 1 or 2 representing data ambient1_unspliced or ambient2_unspliced, respectively.
##@return: Ordered_Data_d_au
Vector_ambient_unspliced <- function( inputFile, ambientFile ) {
  if ( ambientFile == 1 ) {
    Ordered_Data_d_au <- sort( inputFile$ambient1_unspliced, decreasing = FALSE )
  } else if ( ambientFile == 2 ) {
    Ordered_Data_d_au <- sort( inputFile$ambient2_unspliced, decreasing = FALSE )
  } else {
    stop( "ambientFile parameter must be numeric 1 or 2." )
  }
  return( Ordered_Data_d_au * 10000)
}

#--------------------------------------------------------------------------------------------------------------

###Function to find all genes in a given droplet that indicate expression at their ambient RNA level
##@param: input - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile - Indication of ambient file 1 or 2
##@param: droplet - Indication of droplet being tested
Spliced_Processing <- function( input, ambientFile, droplet ){

  spliced_matrix <- input$spliced_cells_matrix


  ##Ordering the data found within trimmed file of spliced and unspliced ambient genes
  ambient <- Vector_ambient_spliced( input, ambientFile)


  ##Create a new vector that contains every 50th value to create a more evenly distributed range of expression percentages
  tempamb <- ambient[ ambient > 0 ] #Finds all the genes that only contain percentages greater than 0
  landmarks <- tempamb[seq(1, length(tempamb), length.out = 200)]


  ##Percentages of expression per droplet
  expr_mat = spliced_matrix
  expr_mat = t( t( 10000 * spliced_matrix ) / as.numeric( colSums( spliced_matrix, na.rm = FALSE, dims = 1L ) ) )
  expr_mat = expr_mat[expr_mat[ , droplet] > 0,]
  expr_mat[ is.nan( expr_mat ) ] = 0


  #Readjust landmarks to only show those genes that are present in the expression matrix of all spliced genes
  landmarks <- landmarks[names(landmarks) %in% rownames(expr_mat)]

  #If 0 or 1 landmarks are left after filtering, then dataset will not be accurately compared, this operation will return NA
  if(length(landmarks) <= 3){
    return(NA)
  }


  ###Does comparison between expression matrix and the expression matrix at specific landmarks
  spliced_dist_expression <- apply( expr_mat, 1, function(x)drop_distance_vector( gene = x, matr = expr_mat, lnd = landmarks, drop = droplet ) ) #Takes a couple of seconds to run


  ##Makes the ambient array contain only those genes that are also present in the expression matrix of all spliced genes
  ambient <- ambient[names(ambient) %in% rownames(expr_mat)]


  ###Ambient comparison
  ambient_dist_expression <- sapply(ambient, function(x)amb_distance_vector( gene = x, lnd = landmarks))

  ambient_dist_expression <- ambient_dist_expression[ , match( colnames(spliced_dist_expression), colnames(ambient_dist_expression) ) ]


  ###Perform correlation between the expression values
  cor_spliced_ambient <- c()

  #for(x in 1:ncol(spliced_dist_expression)){
   # corx <- setNames(cor( spliced_dist_expression[ , x ], ambient_dist_expression[ , colnames( spliced_dist_expression )[ x ] ] , method = "spearman" ), colnames(spliced_dist_expression)[x])
    #cor_spliced_ambient <- c( cor_spliced_ambient, corx )
  #}

  for(x in 1:ncol(spliced_dist_expression)){
    if(length(unique(spliced_dist_expression[,x])) == 1) {
      return(NA)
    } else {
      corx <- setNames(cor( spliced_dist_expression[ , x ], ambient_dist_expression[ , colnames( spliced_dist_expression )[ x ] ] , method = "spearman" ), colnames(spliced_dist_expression)[x])
    }
    cor_spliced_ambient <- c( cor_spliced_ambient, corx )
  }

  #Checks to see that all correlations are different; else Mclust will not differentiate
  if(length(unique(cor_spliced_ambient)) == 1){

    return(NA)

  }

  threshold <- clustering(as.numeric(cor_spliced_ambient))

  genes_above_thresh <- list(setNames(threshold, droplet), cor_spliced_ambient[cor_spliced_ambient >= threshold])

  print(droplet)

  return(genes_above_thresh)
}


###Function to
##@param: input - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile -
##@param: droplet -
Unspliced_Processing <- function( input, ambientFile, droplet ){

  unspliced_matrix <- input$unspliced_cells_matrix #Most of these values are 0


  ##Ordering the data found within trimmed file of spliced and unspliced ambient genes
  ambient <- Vector_ambient_unspliced( input, ambientFile)


  ##Create a new vector that contains 200 values over 0 to create a more evenly distributed range of expression percentages
  tempamb <- ambient[ ambient > 0 ] #Finds all the genes that only contain percentages greater than 0
  landmarks <- tempamb[seq(1, length(tempamb), length.out = 200)]


  ##Percentages of expression per droplet
  un_expr_mat = unspliced_matrix
  un_expr_mat = t( t( 10000 * unspliced_matrix ) / as.numeric( colSums( unspliced_matrix, na.rm = FALSE, dims = 1L ) ) )
  un_expr_mat = un_expr_mat[un_expr_mat[ , droplet] > 0,]
  un_expr_mat[ is.nan( un_expr_mat ) ] = 0


  #Readjust landmarks to only show those genes that are present in the expression matrix of all spliced genes
  landmarks <- landmarks[names(landmarks) %in% rownames(un_expr_mat)]


  ###Does comparison between expression matrix and the expression matrix at specific landmarks
  unspliced_dist_expression <- apply( un_expr_mat, 1, function(x)drop_distance_vector( gene = x, matr = un_expr_mat, lnd = landmarks, drop = droplet ) ) #Takes a couple of seconds to run


  ##Makes the ambient array contain only those genes that are also present in the expression matrix of all spliced genes
  ambient <- ambient[names(ambient) %in% rownames(un_expr_mat)]


  ###Ambient comparison
  ambient_dist_expression <- sapply(ambient, function(x)amb_distance_vector(gene = x, lnd = landmarks))

  ambient_dist_expression <- ambient_dist_expression[ , match( colnames(unspliced_dist_expression), colnames(ambient_dist_expression) ) ]


  ###Perform correlation between the expression values
  cor_unspliced_ambient <- c()

  #for(x in 1:ncol(unspliced_dist_expression)){
   # corx <- setNames(cor( unspliced_dist_expression[ , x ], ambient_dist_expression[ , colnames( unspliced_dist_expression )[ x ] ] , method = "spearman" ), colnames(unspliced_dist_expression)[x])
    #cor_unspliced_ambient <- c( cor_unspliced_ambient, corx )
  #}

  for(x in 1:ncol(spliced_dist_expression)){
    if(length(unique(spliced_dist_expression[,x])) == 1) {
      #corx <- setNames(0, colnames(spliced_dist_expression)[x])
      x <- x + 1
    } else {
      corx <- setNames(cor( spliced_dist_expression[ , x ], ambient_dist_expression[ , colnames( spliced_dist_expression )[ x ] ] , method = "spearman" ), colnames(spliced_dist_expression)[x])
    }
    cor_spliced_ambient <- c( cor_spliced_ambient, corx )
  }

  if(length(unique(cor_spliced_ambient)) == 1){

    return(NA)

  }

  threshold <- clustering(as.numeric(cor_unspliced_ambient))

  genes_above_thresh <- cor_unspliced_ambient[cor_unspliced_ambient >= threshold]

  return(genes_above_thresh)
}

Sys.time()
unspliced_goi <- Unspliced_Processing(input, 1, 1) #significantly less are returned as most data falls below 0, and are excluded from dataset
Sys.time()

#----------------------------------------------------------------------------------------------------------

###Function to calculate the euclidean distances between the expression percentage with those of the landmark genes of a given droplet
##@param: gene - a row of the expression matrix
##@param: matr -
##@param: lnd -
##@param: drop -
drop_distance_vector <- function(gene, matr, lnd, drop ) {

  droplet_vec <- c()

  droplet_vec <- dist( gene[drop], matr[names(lnd), drop])

  return(droplet_vec)
}


###Function to calculate the euclidean distances between the expression percentage with those of the landmark genes of a given droplet
##@param: gene -
##@param: lnd -
amb_distance_vector <- function(gene, lnd){

  ambient_vec <- c()

  ambient_vec <- dist( gene, lnd )

  return(ambient_vec)
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------

###Function to
##@param: correlation_vector -
clustering <- function(correlation_vector){

  fit <- Mclust(correlation_vector, G=1:2)

  #h <- hist(fit$data, prob=TRUE, col="grey85", xlab="Score", main="Correlation of RNA Clusters To Their Ambient Expression Levels", breaks=20)
  #cols <- c("forestgreen", "darkorange", "blue")

  #for (i in 1:fit$G) {
    #if (is.na(fit$parameters$variance$sigmasq[i])) {
      #fit$parameters$variance$sigmasq[i] <- fit$parameters$variance$sigmasq[1]
    #}
    #distr_fxn <- function(x) { dnorm(x, mean=fit$parameters$mean[i], sd=sqrt(fit$parameters$variance$sigmasq[i]))*fit$parameters$pro[i] }
    #curve(distr_fxn(x), col=cols[i], add=TRUE, lwd=1.5)
  #}

  return(min(correlation_vector[fit$classification == 2]))
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

Sys.time()
droplet_thresholds <- c()
for(y in 1:ncol(input$spliced_cells_matrix)){
  gois <- Spliced_Processing(input, 1, y)
  droplet_thresholds <- c(droplet_thresholds, gois[1])
}
Sys.time()

dataDrop <- droplet_thresholds

dataDrop[is.infinite(as.numeric(dataDrop))] <- NA
cleaned_data <- dataDrop[!is.na(as.numeric(dataDrop))]

b <- boxplot(as.numeric(cleaned_data), col="orange", xlab="Droplets", ylab = "Correlation Values", main="", border="brown", notch = TRUE)
points(mean(as.numeric(cleaned_data)), col = "blue", pch = 20)

###Function to perform seurat comparison on all cell matrices to find marker genes, noting if the genes returned in SplicedProcessing occur within them
##@param: inputFile - - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: genesOfInterest - List of all genes that exceed the given threshold value, indicating that they are "ambiently" expressed
markerMatch <- function(inputFile, genesOfInterest){
  all_expression <- inputFile$spliced_cells_matrix+inputFile$unspliced_cells_matrix

  seur_obj <- CreateSeuratObject(counts = all_expression, min.cells = 3, min.features = 300)
  seur_obj <- NormalizeData(seur_obj)
  seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seur_obj)
  seur_obj <- ScaleData(seur_obj, features = all.genes)
  seur_obj  <- RunPCA(seur_obj , features = VariableFeatures(object = seur_obj))
  seur_obj <- FindNeighbors(seur_obj, dims = 1:10)
  seur_obj <- FindClusters(seur_obj, resolution = 0.5)
  markers <- FindAllMarkers(seur_obj, only.pos = TRUE)

  matchingGenes <- sum(names(genesOfInterest)%in%markers$gene)

  return(matchingGenes)
}
