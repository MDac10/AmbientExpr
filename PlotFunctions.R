#'##Megan Da Costa
#'##Supervisor: Tallulah Andrews

#'#--------------------------------------------------------------------------------------------------------------

#'
#'@export
#'###Plotting function to show distribution of threshold values seen across both Spliced and Unspliced data in Box Plot form
#'##@param: Spliced_Data - The output from AllSplicedDroplets function
#'##@param: Unspliced_Data - The output from AllUnsplicedDroplets function
ThresholdBoxPlot <- function( Spliced_Data, Unspliced_Data ) {
  Spliced_Thresholds <- Spliced_Data[[1]]
  Unspliced_Threshold <- Unspliced_Data[[1]]
  df <- data.frame( Spliced_Thresholds, Unspliced_Thresholds )
  #Box plot distribution - Showing mean values
  comparativebox <- boxplot( df, na.rm = TRUE, xlab = "Gene Types", ylab = "Correlation Threshold" )
  points( colMeans( df, na.rm=TRUE ), col = "red", pch = "+", cex = 2 )
}


#'
#'@export
#'###Plotting function to show distribution of threshold values seen across both Spliced and Unspliced data in Violin Plot form
#'##@param: Spliced_Data - The output from AllSplicedDroplets function
#'##@param: Unspliced_Data - The output from AllUnsplicedDroplets function
ThresholdViolinPlot <- function( Spliced_Data, Unspliced_Data ) {
  Spliced_Thresholds <- Spliced_Data[[1]]
  Unspliced_Threshold <- Unspliced_Data[[1]]
  df <- data.frame( Spliced_Thresholds, Unspliced_Thresholds )
  #Violin plot showing mean and median values
  comparativev <- vioplot( df, col=c( "cadetblue3", "mediumpurple2" ), border = c( 1, 1 ), xlab = "Gene Types", ylab = "Threshold Value" ) + title( main = "Distribution of Threshold Correlation Values Across All Droplets" )
  points( colMeans( df, na.rm = TRUE ), pch = 15, col = "red", cex = 1.0 )
  legend( "bottomright", pch = c( 21, 19 ), col = c( "black", "red" ), bg = "white", legend = c( "Median", "Mean" ), cex = 1.0 )
  text( c( 0.63, 0.7 ), labels = round( colMeans( df, na.rm = TRUE ), 4 ), cex = 0.9, font = 2, col = "black" )
}


#'
#'@export
#'###Plotting function to show distribution of number of ambient genes returned per droplet across both Spliced and Unspliced data in Scatter Plot form
#'##@param: Spliced_Data - The output from AllSplicedDroplets function
#'##@param: Unspliced_Data - The output from AllUnsplicedDroplets function
NumAmbPlot <- function( Spliced_Data, Unspliced_Data ) {
  #Collection of total number of ambient genes returned per droplet of Spliced data sample
  sg <- c()
  for( x in 1:length( Spliced_Data[[2]] ) ) {
    sg <- c( sg, length( Spliced_Data[[2]][[x]] ) )
  }
  #Collection of total number of ambient genes returned per droplet of Unspliced data sample
  ug <- c()
  for( x in 1:length( Unspliced_Data[[2]] ) ) {
    ug <- c( ug, length( Unspliced_Data[[2]][[x]] ) )
  }
  #Plot all Spliced and Unspliced data against each other
  par( mfrow = c(1:2) )
  plot( sg, ylim = c(0,1500), xlab = "Droplet", ylab="Number of Ambient genes Returned For Spliced Dataset" )
  abline( h=mean( sg, na.rm=TRUE ), col="red" )
  plot( ug, ylim = c(0,1500), xlab = "Droplet", ylab="Number of Ambient genes Returned For Unspliced Dataset" )
  abline( h=mean( ug, na.rm=TRUE ), col="red" )
}

#'#--------------------------------------------------------------------------------------------------------------
