#' Creates a subnetwork from a COSIFER network
#'
#' This function loads a network file produced using COSIFER ( with extension .csv.gz) as a dataframe.
#' It takes the network file path directory as an input, and stores the subnetwork in the output directory.
#' If multiple .csv.gz files are present, the user can select a specific file using the idx parameter.
#' Additional parameters are passed to the function "filter_edges" for the creation of the subnetwork.
#' The function stores the subnetwork, the largest connected component of the network as well as
#' its coexpression matrix as an R session (.Rdata extension); the largest connected component edge weighet dataframe
#' in a .csv.gz file and as an igraph object in a .sif format (also to easy visualization into cytoscape)
#'
#'
#' @param file Path to the input file
#' @param outDirectory Path to the output directory
#' @param idx Optional: index of the network file when multiple are present
#' @return Stores the subnetwork and its largest connected component into different formats for further use
#' @export
create_subnetwork <- function(directory, outDirectory, idx=NULL, threshold=NULL, n.edges=10^5){

  ## create output directory
  if (!dir.exists(outDirectory))
    dir.create(outDirectory, recursive = T)

  # optionally, get file index (for multithreading)
  fileIdx = NULL
  if (length(args) >= 3)  fileIdx = as.numeric(args[3])
  # optionally, get desired weight threshold (for edge filtering)
  weight_threshold = NULL
  if (length(args) >= 4)  weight_threshold = as.numeric(args[4])

  # read files
  files = list.files(directory, pattern = "*.gz", full.names = T)
  if (!is.null(fileIdx)) files = files[fileIdx]

  for (file in files){
   filename = gsub("\\.csv.*$", "", basename(file))
   message("loading network for method ", filename)
   output = file.path(outDirectory, paste0(filename, ".Rdata"))
   if (!file.exists(output)) {
     tic()
     net = read_big_files(file)

     ## create subnetwork
     sub.net <- filter_edges(net, n.edges=n.edges)
     plot_density(net, sub.net, filename = filename, dir=outDirectory)

     sub.net_igraph <- graph_from_data_frame(d = sub.net, directed = FALSE)
     dis.components <- components(sub.net_igraph)
     gp.components <- groups(dis.components)
     largest.conn.comp <- subgraph(sub.net_igraph, gp.components[[which.max(dis.components$csize)]])
     coexpr <- as_adjacency_matrix(largest.conn.comp, type = c("both"),
                                    names = TRUE, attr = "weight",
                                    sparse = igraph_opt("sparsematrices"))
     coexpr <- as.matrix(coexpr)

     save(coexpr, sub.net_igraph, largest.conn.comp, gp.components, file = output)
     BioNet::saveNetwork(largest.conn.comp, file=file.path(outDirectory, filename), name="largest connected component", type="sif")
     write.csv.gz(sub.net, file= paste0(outDirectory, "reduced_", basename(file)), row.names = TRUE)
     toc()
    }
  }
}


