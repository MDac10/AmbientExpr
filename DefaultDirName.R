install.packages("Matrix")
install.packages("proxy")
install.packages("mclust")
install.packages("ggplot2")

library(Matrix)
library(proxy)
library(mclust)
library(ggplot2)

setwd('C:\\Users\\megda\\Thesis_4460Z\\ExampleInput')

input <- readRDS('KorEar_Cells1_ForMegan_inputtrimmed.rds')

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

testFunc <- function( input, ambientFile, droplet ){
  spliced_matrix <- input$spliced_cells_matrix


##Ordering the data found within trimmed file of spliced and unspliced ambient genes
  ambient <- Vector_ambient_spliced( input, ambientFile)


##Create a new vector that contains every 50th value to create a more evenly distributed range of expression percentages
  tempamb <- ambient[ ambient > 0 ]
#landmarks <- ambient[ ambient > 0 ]#Finds all the genes that only contain percentages greater than 0
  landmarks <- tempamb[seq(1, length(tempamb), length.out = 200)]


##Percentages of expression per droplet
  expr_mat = spliced_matrix
  expr_mat = t( t( 10000 * spliced_matrix ) / as.numeric( colSums( spliced_matrix, na.rm = FALSE, dims = 1L ) ) )
  expr_mat = expr_mat[expr_mat[ , droplet] > 0,]
  expr_mat[ is.nan( expr_mat ) ] = 0



#Readjust landmarks to only show those genes that are present in the expression matrix of all spliced genes
  landmarks <- landmarks[names(landmarks) %in% rownames(expr_mat)]

  if(length(landmarks) <= 1){

    return(NA)
  }

  cor_spliced_ambient <- c()

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

  threshold <- clustering(as.numeric(cor_spliced_ambient))

  genes_above_thresh <- list(setNames(threshold, droplet), cor_spliced_ambient[cor_spliced_ambient >= threshold])

  return(expr_mat)
}

expr_mat <- testFunc(input, 1, 1)
