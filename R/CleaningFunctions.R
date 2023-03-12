##Megan Da Costa
##Supervisor: Tallulah Andrews

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

###Function to take output of droplets and their corresponding ambient genes from both spliced and unspliced samples and create a readable version for calculation functions to use
##@param: sgoi - A list of all droplets corresponding ambient genes from spliced data, can be taken from the second list value from AllSplicedDroplets function or a basic run of Spliced_Processing function
##This parameter should be of form splicedData[[2]] if splicedData is a list of lists and splicedData[2] for a singular list, or else this function will not be able to read the data correctly
##@param: ugoi - A list of all droplets corresponding ambient genes from unspliced data, can be taken from the second list value from AllUnsplicedDroplets function or a basic run of Unspliced_Processing function
##This parameter should be of form unsplicedData[[2]] if unsplicedData is a list of lists and unsplicedData[2] for a singular list, or else this function will not be able to read the data correctly
##@param: drops - A list of all droplets taken from the first list value from either AllSplicedDroplets or AllUnsplicedDroplets functions or from a basic run of either Spliced_Processing or Unspliced_Processing
##This parameter should be of form splicedData[[1]] no matter what type of file is being run
##@return: sAnduAmbient - A list where each entry is named after the droplet it corresponds to, and has the names of all identified ambient genes from the spliced and unspliced data
CombineAndClean <- function(sgoi, ugoi, drops){
  sAnduAmbient <- list()
  for(x in 1:length(sgoi)){
    sNames <- c(names(sgoi[[x]]))
    uNames <- c(names(ugoi[[x]]))
    vec <- c(sNames, uNames)
    vec <- vec[!duplicated(vec)]
    sAnduAmbient <- append(sAnduAmbient, list(vec))
  }
  sAnduAmbient <- setNames(sAnduAmbient, names(drops))
  return(sAnduAmbient)
}


###Function to find all genes that were identified as ambient across the dataset
##@param: cleanedData - A list where each entry is named after the droplet it corresponds to, and has the names of all identified ambient genes from the spliced and unspliced data
##@return: genes - A vector of all genes that were identified as ambient in the dataset
AllAmbientGenes <- function(cleanedData){
  genes <- c()
  for(x in 1:length(cleanedData)){
    genes <- c(genes, cleanedData[[x]])
    genes <- genes[!duplicated(genes)]
  }
  return(genes)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
