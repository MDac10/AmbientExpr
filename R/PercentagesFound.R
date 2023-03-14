#'##Megan Da Costa
#'##Supervisor: Tallulah Andrews

#'#--------------------------------------------------------------------------------------------------------------------------------------------------------------

#'
#'@export
#'###Function to determine the percentage of identified ambient genes in a droplet that are found within its corresponding cluster of marker genes
#'##@param: seur_obj - A Seurat object containing all genes, spliced,and unspliced expressions, and their corresponding informtaion (clusters, id, neighbours etc.)
#'##@param: markers - A dataframe contaning all marker genes (those suspected to be ambient) in the total expression dataset
#'##@param: gois - The vector containing the droplet names with its threshold value for determining ambient expression, and the list of ambient genes identified in the droplet
#'##@param: verbose - A boolean variable where, when set to true will print the results in the console; default setting is false
#'##@return: ambPerCluster - A list that keeps track of all ambient genes found within each cluster
PerClusterAmbient <- function( seur_obj, markers, gois, verbose=FALSE ) {
  clusts <- list()
  #Using the number of clusters
  for( x in levels( seur_obj@meta.data$seurat_clusters ) ) {
    droplets_in_this_cluster <- rownames( seur_obj@meta.data )[ seur_obj@meta.data$seurat_clusters == x ]
    this_cluster_marker <- markers[ markers$cluster == x,"gene" ]
    cl <- list( "droplets" = droplets_in_this_cluster, "markers" = this_cluster_marker )
    clusts[[length( clusts ) + 1]] <- cl
  }
  ambPerCluster <- list()
  allGeneMatches <- c()
  for( x in 1:length( clusts ) ) {
    perClust = c()
    for( y in 1:length( gois ) ) {
      if( names( gois[y] )%in%clusts[[x]]$droplets ) {
        perClust <- c( perClust, gois[[y]] )
        perClust <- perClust[!duplicated( perClust )]
      }
    }
    ambPerCluster <- append( ambPerCluster, list( perClust ) )
    matchingGenes <- round( sum( perClust%in%clusts[[x]]$markers )/length( perClust ),4 )
    allGeneMatches <- c( allGeneMatches, matchingGenes )
    if( is.na( matchingGenes ) ) {
      next
    }
    if( verbose ) {
      a <- "For cluster "
      b <- " the percent ambient genes found in Seurat markers is: "
      c <- " x100%"
      cat( a,x,b,matchingGenes,c, '\n' )
    }
  }
  ambPerCluster <- setNames( ambPerCluster, levels( seur_obj@meta.data$seurat_clusters ) )
  return( ambPerCluster )
}


#'--------------------------------------------OBSOLETE-----------------------------------------------
#'###Function to find percentages of marker genes compared to ambient genes
#'##@param: inputFile -  of form readRDS('xyz_Cells__inputtrimmed.rds')
#'##@param: genesOfInterest - An array of the genes identified as ambient
#'##@return: matchingGenes - A percentage of the ambient genes identified by the algorithm that are found as marker genes in the dataset
AllPercentAmbient <- function( foundMarkers, genesOfInterest ) {
  matchingGenes <- sum( genesOfInterest%in%foundMarkers$gene )/length( genesOfInterest ) #percent of the ambient genes identified that are marker genes
  matching <- sum( genesOfInterest%in%foundMarkers$gene )/nrow( foundMarkers ) #percent of the marker genes that are also identified as ambient by our algorithm
  markersPercent <- nrow( foundMarkers )/nrow( input$spliced_cells_matrix+input$unspliced_cells_matrix ) #percent of rows total number of genes that are marker genes
  return( matchingGenes )
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
