##Megan Da Costa
##Supervisor: Tallulah Andrews

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

###Function to run Spliced_Processing over all droplets in the dataset
##@param: inputFile - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile - the attachedment after inputFile$, numeric 1 or 2 representing data ambient1_spliced or ambient2_spliced, respectively.
##@return: aud - A large list, where the first entry is the names of all droplets tested, and the second contains lists each of the droplets corresponding ambient genes that were identified
AllSplicedDroplets <- function(inputFile, ambientFile){
  print("Process started......")
  Sys.time()
  droplet_thresholds <- c()
  goi <- c()
  for(x in 1:ncol(inputFile$spliced_cells_matrix)){
    gois <- Spliced_Processing(inputFile, ambientFile, x)
    droplet_thresholds <- c(droplet_thresholds, gois[[1]])
    goi <- c(goi, gois[2])
    print(x)
  }
  print("Process ended......")
  Sys.time()
  aud <- list(droplet_thresholds, goi)
  return(aud)
}


###Function to run Unspliced_Processing over all droplets in the dataset
##@param: inputFile - of form readRDS('xyz_Cells__inputtrimmed.rds')
##@param: ambientFile - the attachedment after inputFile$, numeric 1 or 2 representing data ambient1_spliced or ambient2_spliced, respectively.
##@return: aud - A large list, where the first entry is the names of all droplets tested, and the second contains lists each of the droplets corresponding ambient genes that were identified
AllUnsplicedDroplets <- function(inputFile, ambientFile){
  print("Process started......")
  Sys.time()
  droplet_thresholds <- c()
  goi <- c()
  for(x in 1:ncol(inputFile$unspliced_cells_matrix)){
    gois <- Unspliced_Processing(inputFile, ambientFile, x)
    droplet_thresholds <- c(droplet_thresholds, gois[[1]])
    goi <- c(goi, gois[2])
    print(x)
  }
  print("Process ended......")
  Sys.time()
  aud <- list(droplet_thresholds, goi)
  return(aud)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
