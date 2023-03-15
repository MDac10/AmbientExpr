# AmbientExpr

In the R script where you wish to use this algorithm...
Load the contents of this algorithm, packages must be installed and are loaded:

  > install.packages("Matrix")
  > install.packages("proxy")
  > install.packages("mclust")
  > install.packages("vioplot")
  > library(Matrix) ...

The seurat package is also required to be able to use functions in this code package:

  > require(Seurat)

Either set the working directory to the folder that contains your input file or incorporate it into:

  > input <- readRDS('exampleFolder\\exampleInputFile.rds')

This input file should be preprocessed containing files, minimum:
* unspliced_cells_matrix
* spliced_cells_matrix
* ambient1_spliced
* ambient1_unspliced
* ambient2_spliced
* ambient2_unspliced

## Function Specifics

### Processing Functions

Inputs for Spliced_Processing and Unspliced_Processing needs an indication of which ambient files to use (1 or 2):
* Ambient file 1 contains all data including droplets with more RNA (>= 15 mRNA molecules detected) but not a valid cell - identified as ambient RNA
* Ambient file 2 contains all data including droplets with low RNA (< 15 mRNA molecules detected) - identified as ambient RNA
