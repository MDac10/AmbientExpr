##Megan Da Costa
##Supervisor: Tallulah Andrews

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

###Function to calculate the euclidean distances between the expression percentage with those of the landmark genes of a given droplet
##@param: gene - A row of the expression matrix
##@param: matr - The expression matrix being used (either spliced or unspliced)
##@param: lnd - A vector containing all landmark genes for a given droplet
##@param: drop - The droplet being considered in numerical or character format
##@return: droplet_vec - A vector containing the distances between the expression percentage of a gene at a given droplet and all of the landmark genes at the same droplet
drop_distance_vector <- function(gene, matr, lnd, drop ) {

  droplet_vec <- c()

  droplet_vec <- dist( gene[drop], matr[names(lnd), drop])

  return(droplet_vec)
}


###Function to calculate the euclidean distances between the expression percentage with those of the landmark genes of a given droplet
##@param: gene - An entry of the ambient genes vector
##@param: lnd - A vector containing all landmark genes for a given droplet
##@return: ambient_vec - A vector containing the distances between the expression percentages of an ambient identified gene and all of the landmark genes within the ambiently defined genes
amb_distance_vector <- function(gene, lnd){

  ambient_vec <- c()

  ambient_vec <- dist( gene, lnd )

  return(ambient_vec)
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------

###Function to cluster the correlation values for a droplet
##@param: correlation_vector - A vector containing all the genes
##@return: min(correlation_vector[fit$classification == 2])
##The comments can be removed from the histogram section below to see the distribution of correlation values for  droplet
clustering <- function(correlation_vector){

  fit <- Mclust(correlation_vector, G=1:2)

  ##For plotting all correlation values to show distribution in the form of a scatter plot
  #plot(fit$data, col="blue", xlab="Genes in Expression Matrix", ylab="Correlation To Ambient Expression", main="Clustering for Droplet 1 Correlations")

  ##For plotting all correlation values over a histogram to show bimodal distribution of data, where the intersection of the distributions acts as a threshold value for what is to be considered ambient
  #h <- hist(fit$data, prob=TRUE, col="grey85", xlab="Score", main="Correlation of RNA Clusters To Their Ambient Expression Levels", breaks=20)
  #cols <- c("forestgreen", "darkorange", "blue")
  #for (i in 1:fit$G) {
  # if (is.na(fit$parameters$variance$sigmasq[i])) {
  #  fit$parameters$variance$sigmasq[i] <- fit$parameters$variance$sigmasq[1]
  #}
  #distr_fxn <- function(x) { dnorm(x, mean=fit$parameters$mean[i], sd=sqrt(fit$parameters$variance$sigmasq[i]))*fit$parameters$pro[i] }
  #curve(distr_fxn(x), col=cols[i], add=TRUE, lwd=1.5)
  #}

  return(fit)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
