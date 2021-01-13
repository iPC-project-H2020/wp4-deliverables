## MCL CLUSTERING AND FUNCTIONAL REPRESENTATION

#' Performs MCL clustering with inference of the inflation parameter
#'
#' This function loads .sif igraph networks present in the input directory,
#' after inflation parameter tuning for the MCL clustering (such that the number of clusters
#' with more than 10 genes is the max), it stores the clustering for the optimized inflation value
#' as a vertex attributes for the igraph object, the updated igraph object in a .sif format, and the
#' session as an Rdata object.
#'
#' @param directory Path to the input file
#' @param outDirectory Path to the output directory
#' @return Stores the MCL results as an attribute of the igraph network and the updated igraph network
#' @export
mcl_clustering <- function(directory, outDirectory){

    if (length(args) == 0)  stop("Must set input file directory", call. = F)
    directory = args[1]
    if (length(args) == 1)  stop("Must set output file directory", call. = F)
    outDirectory = args[2]

    files = list.files(directory, pattern = "*\\.sif$", full.names = T)


    for (file in files){
        ## perform comunity detection on the subnetwork
        filename = gsub("*\\.sif$", "", basename(file))
        load(paste0(gsub("*\\.sif$", "", file), ".Rdata"))

        cl.gsm <- tunning.inflation(file, outDirectory)
        n.cl10 <- cl.gsm$n10
        cl.gsm <- cl.gsm$mcl.clustering

        message("clustering for method ", filename, ": COMPLETED")
        message("Number of communities with n.genes > 10:", n.cl10)

        V(largest.conn.comp)$community <- cl.gsm$Cluster
        BioNet::saveNetwork(largest.conn.comp, file=file.path(outDirectory, filename), name="highest connected component", type="sif")
        save.image(paste0(gsub("*\\.sif$", "", file), ".Rdata"))
    }
}
