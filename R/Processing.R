##Megan Da Costa
##Supervisor: Tallulah Andrews

#--------------------------------------------------------------------------------------------------------------

###Function to find all spliced genes in a given droplet that indicate expression at their ambient RNA level
##@param: input - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile - Indication of ambient file 1 or 2
##@param: droplet - Indication of droplet being tested in numerical or character format
##@return: genes_above_thresh - a list containing the threshold value calculated, and a vector of all spliced genes found to be above the threshold value
Spliced_Processing <- function( input, ambientFile, droplet ){
  spliced_matrix <- input$spliced_cells_matrix
  if(is.numeric(droplet)){
    droplet <- colnames(spliced_matrix)[droplet]
  }
  ##Ordering the data found within trimmed file of spliced and unspliced ambient genes
  ambient <- Vector_ambient_spliced( input, ambientFile)
  ##Create a new vector that contains every 50th value to create a more evenly distributed range of expression percentages
  tempamb <- ambient[ ambient > 0 ] #Finds all the genes that only contain percentages greater than 0
  landmarks <- tempamb[seq(1, length(tempamb), length.out = 200)]
  ##Percentages of expression per droplet
  expr_mat = spliced_matrix
  expr_mat = t( t( 10000 * spliced_matrix ) / as.numeric( colSums( spliced_matrix, na.rm = FALSE, dims = 1L ) ) )
  expr_mat[ is.nan( expr_mat ) ] = 0
  expr_mat = expr_mat[expr_mat[ , droplet] > 0,] #For efficiency of algorithm, all variables that are set at 0 are removed
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
  fitValues <- clustering(as.numeric(cor_spliced_ambient))
  threshold <- suppressWarnings(min(cor_spliced_ambient[fitValues$classification == 2]))
  if(is.infinite(threshold)){
    return(NA)
  }
  genes_above_thresh <- list(setNames(threshold, droplet), cor_spliced_ambient[cor_spliced_ambient >= threshold])
  print(droplet)
  return(genes_above_thresh)
}



###Function to find all unspliced genes in a given droplet that indicate expression at their ambient RNA level
##@param: input - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile - Indication of ambient file 1 or 2
##@param: droplet -Indication of droplet being tested in numerical or character format
##@return: genes_above_thresh - a list containing the threshold value calculated, and a vector of all spliced genes found to be above the threshold value
Unspliced_Processing <- function( input, ambientFile, droplet ){
  unspliced_matrix <- input$unspliced_cells_matrix #Most of these values are likely to be 0
  if(is.numeric(droplet)){
    droplet <- colnames(unspliced_matrix)[droplet]
  }
  ##Ordering the data found within trimmed file of spliced and unspliced ambient genes
  ambient <- Vector_ambient_unspliced( input, ambientFile )
  ##Create a new vector that contains 200 values over 0 to create a more evenly distributed range of expression percentages
  tempamb <- ambient[ ambient > 0 ] #Finds all the genes that only contain percentages greater than 0
  landmarks <- tempamb[seq(1, length(tempamb), length.out = 200)]
  ##Percentages of expression per droplet
  un_expr_mat = unspliced_matrix
  un_expr_mat = t( t( 10000 * unspliced_matrix ) / as.numeric( colSums( unspliced_matrix, na.rm = FALSE, dims = 1L ) ) )
  un_expr_mat[ is.nan( un_expr_mat ) ] = 0
  un_expr_mat = un_expr_mat[un_expr_mat[ , droplet] > 0,]
  #Readjust landmarks to only show those genes that are present in the expression matrix of all spliced genes
  landmarks <- landmarks[names(landmarks) %in% rownames(un_expr_mat)]
  #If 0 or 1 landmarks are left after filtering, then dataset will not be accurately compared, this operation will return NA
  if(length(landmarks) <= 3){
    return(NA)
  }
  ###Does comparison between expression matrix and the expression matrix at specific landmarks
  unspliced_dist_expression <- apply( un_expr_mat, 1, function(x)drop_distance_vector( gene = x, matr = un_expr_mat, lnd = landmarks, drop = droplet ) ) #Takes a couple of seconds to run
  ##Makes the ambient array contain only those genes that are also present in the expression matrix of all spliced genes
  ambient <- ambient[names(ambient) %in% rownames(un_expr_mat)]
  ###Ambient comparison
  ambient_dist_expression <- sapply(ambient, function(x)amb_distance_vector(gene = x, lnd = landmarks))
  ambient_dist_expression <- ambient_dist_expression[ , match( colnames(unspliced_dist_expression), colnames(ambient_dist_expression) ) ]
  ###Perform correlation between the expression values
  cor_unspliced_ambient <- c()
  for(x in 1:ncol(unspliced_dist_expression)){
    if(length(unique(unspliced_dist_expression[,x])) == 1) {
      return(NA)
    } else {
      corx <- setNames(cor( unspliced_dist_expression[ , x ], ambient_dist_expression[ , colnames( unspliced_dist_expression )[ x ] ] , method = "spearman" ), colnames(unspliced_dist_expression)[x])
    }
    cor_unspliced_ambient <- c( cor_unspliced_ambient, corx )
  }
  if(length(unique(cor_unspliced_ambient)) == 1){
    return(NA)
  }
  fitValues <- clustering(as.numeric(cor_unspliced_ambient))
  threshold <- suppressWarnings(min(cor_unspliced_ambient[fitValues$classification == 2]))
  if(is.infinite(threshold)){
    return(NA)
  }
  genes_above_thresh <- list(setNames(threshold, droplet), cor_unspliced_ambient[cor_unspliced_ambient >= threshold])
  print(droplet)
  return(genes_above_thresh)
}

#----------------------------------------------------------------------------------------------------------
