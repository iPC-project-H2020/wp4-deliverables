#' Load a big network
#'
#' This function loads network files produced using COSIFER as a matrix.
#' The column names of such networks are "e1", "e2", "intensity".
#' The function uses this characteristic to import column by column the network.
#' Finally, columns are concatenated and a dataframe is returned.
#'
#' @param file Path to the input file
#' @return A dataframe of the file
#' @export
read_big_files <- function(file){
  net.e1 <- fread(file, select = "e1", data.table = FALSE)
  net.e2 <- fread(file, select = "e2", data.table = FALSE)
  net.intensity <- fread(file, select = "intensity", data.table = FALSE)

  net <- cbind(net.e1, net.e2, net.intensity)
  colnames(net) <- c("from", "to", "weight")
  return(net)
}



#' Filter/remove edges from a network
#'
#' The function takes as input an edge weighted dataframe encoded in its most typical form:
#' with column names "from", "to", "weight".
#' It filters and removes edges either based on a threshold (w.threshold)
#' either based on the number of edges to keep (n.edges),
#' selecting for the ones that have the highest weights; and
#' it returns the reduced edge weighted dataframe.
#'
#' @param df.net An edge weighted dataframe
#' @param w.threshold A weight threshold
#' @param n.edges A number of edges to retain
#' @return A reduced edge weighted dataframe
#' @export
filter_edges <- function(df.net, w.threshold = NULL, n.edges = 10^5){ # Complexity for mcl O(n.nodes^3)
   if (is.null(w.threshold)){
      if (!(n.edges>dim(df.net)[1])) df.net <- df.net[order(df.net$weight, decreasing=TRUE)[1:n.edges], ]
    } else {
      df.net <- df.net[which(df.net$weight>w.threshold), ]
    }
  return(df.net)
}


#' Plot a network
#'
#' The function takes as input an igraph network.
#' If communities are present as vertex attributes of the network
#' nodes are colored based on the community they belong to.
#' If communities are not present, nodes are colored using the default igraph color.
#'
#' @param net An igraph network
#' @param title A title for the plot
#' @return A network plot
#' @export
plot_network <- function(net, title=NULL){
  if (!is.null(V(net)$community)){
    ## define colors
    nb.cols <- unique(V(net)$community)
    colrs <- adjustcolor(colorRampPalette(brewer.pal(8, "Set2"))(nb.cols), alpha=1)
    names(colrs) <- nb.cols
    community.colors <- colrs[as.character(V(net)$community)]
    X11
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    coords <- layout_with_fr(net, weights = E(net)$weight)
    plot(net, rescale=T, asp=0, vertex.size=2, vertex.color=community.colors,
         vertex.label=NA, vertex.frame.color=NA, edge.width= E(net)$weight,
         layout=coords, main=paste0("Method: ", title))
    legend("topright",legend = unique(V(net)$community)[order(unique(V(net)$community))],
           fill = community.colors[order(unique(V(net)$community))], box.col = "white",
           cex=0.6, ncol=3, inset = c(-0.3,0))
  } else {
    coords <- layout_with_fr(net, weights = E(net)$weight)
    plot(net, rescale=T, asp=0, vertex.size=2,
         vertex.label=NA, vertex.frame.color=NA, edge.width= E(net)$weight,
         layout=coords, main=paste0("Method: ", title))
  }
}

#' Plot edge weight density for two edge weighted dataframes
#'
#' The function takes as input two edge weighted dataframes encoded in its most typical form:
#' with column names "from", "to", "weight". It produces the edge weight density plots
#' one beside the other and store it as a pdf file to the given directory.
#'
#' @param df.net An edge weighted dataframe
#' @param df.sub.net An edge weighteddataframe
#' @param filename A pdf filename
#' @param dir An output directory
#' @param main.1 A title for the first density plot
#' @param main.2 A title for the second density plot
#' @return A PDF plot
#' @export
plot_density <- function(df.net, df.sub.net, filename="density", dir = "", main.1="density", main.2="reduced network"){
  pdf(paste0(dir, filename, ".pdf"), width=14, height=7)
  par(mfrow=c(1,2))
  #den.net <- density(df.net$weight)
  #den.sub.net <- density(df.sub.net$weight)
  plot(density(df.net$weight), main=main.1, pch=19, cex=0.3, xlab="Edge weights", ylab="Density")
  plot(density(df.sub.net$weight), main=main.2, pch=19, cex=0.3, xlab="Edge weights", ylab="Density")
  dev.off()
}


#' Tunes the inflation parameter of MCL clustering
#'
#' The function takes as input the directory of a R working space
#' where a coexpression matrix of a network has been stored;
#' and an output directory where to store the clustering results
#' computed with different values of the inflation parameter.
#' Selection of the inflation parameter is performed selecting for the value
#' leading to the maximum number of clusters with more than 10 nodes.
#' The function returns a list with the MCL clustering for the selected
#' inflation parameter and the number of clusters with more than 10 nodes
#'
#' @param file An input directory
#' @param outDirectory An output directory
#' @return A list with the mcl.clustering results and the number of communities with more than 10 genes
#' @export
tunning.inflation <- function(file, outDirectory){
  cl.gsm <- NULL
  filename = gsub("*\\.sif$", "", basename(file))
  load(paste0(gsub("*\\.sif$", "", file), ".Rdata"))

  for (ii in 1:25){
    cl.gsm[[ii]] <- mcl(coexpr, addLoops = FALSE, expansion = 2,
                        inflation = 1+ ii/10 , allow1 = FALSE,
                        max.iter = 10^10, ESM = FALSE)
    message("Interaction ", ii, ": COMPLETED")
  }
  save(cl.gsm, file=paste0(outDirectory, "tune_inflation_",filename, ".RData"))
  lc <- NULL
  for (ii in 1:length(cl.gsm)){
    if(length(cl.gsm[[ii]])==1) {cl.n <- 0
    } else {
      cl.n <- length(which(table(cl.gsm[[ii]]$Cluster)>10))
    }
    lc <- c(lc, cl.n)
  }
  sel.inflation <- which(lc==max(lc))[1]
  return(list("mcl.clustering" = cl.gsm[[sel.inflation]], "n10" = max(lc)))
}
