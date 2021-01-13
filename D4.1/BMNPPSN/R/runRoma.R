#' Computes the eigenGene matrix using rROMA
#'
#' This function loads the expression matrix in the input directory, and the .Rdata object
#' produced by the run_gprofiler function, storing a list of gene sets used by rRoma.R function of
#' the rRoma package.
#' It stores the eigenGene table as a .csv file at the input directory.
#'
#'
#' @param directory Path to the expression matrix
#' @param outDirectory Path to the .RData gene set list
#' @return Stores the eigenGene table as a .csv objet
#' @export
runRoma <- function(directory, outDirectory){

  file = list.files(directory, pattern = "*.csv$", full.names = T)
  if(length(file)>1){
    file <- file[1]
    message("Running ROMA for:", basename(file),", the firts data matrix found in the folder")
  }

  cluster.files = list.files(outDirectory, pattern = "gprofiler.Rdata$", full.names = T)


  ### Attention here the matrix is assumed to have genes on rows and samples on columns
  dat.matrix <- read.csv(file, header = TRUE, row.names = 1)

  for (cluster.list in cluster.files){
    load(cluster.list)
    filename <- gsub("\\_cluster_gprofiler.Rdata$", "", basename(cluster.list))
    rRoma.results <- rRoma.R(ExpressionMatrix = dat.matrix,
                             ModuleList = list.genes.in.cluster,
                             FixedCenter = TRUE, UseParallel = FALSE,
                             ClusType = "FORK", nCores = 15,
                             centerData = TRUE,
                             MinGenes = 11, MaxGenes=3000)
    write.csv2(rRoma.results$SampleMatrix, file=file.path(directory, paste0(filename, "_eigenGenes.csv")),
               quote=FALSE)
  }
}
