###Megan Da Costa
###Supervisor: Tallulah Andrews

#-------------------------------------------------------------------------------------------------------------------------------------------------------------

###Function to organize ambient RNA samples for a given spliced ambient RNA
##@param: inputFile - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile - the attachedment after inputFile$, numeric 1 or 2 representing data ambient1_spliced or ambient2_spliced, respectively.
##@return: Ordered_Data_d_as - an ordered vector of the ambient gene expressions of the spliced set of genes; This vector is multiplied by 10000 so that the data can be compared easily
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
##@return: Ordered_Data_d_au - an ordered vector pf the ambient gene expressions of the unspliced set of genes; This vector is multiplied by 10000 so that the data can be compared easily
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

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
