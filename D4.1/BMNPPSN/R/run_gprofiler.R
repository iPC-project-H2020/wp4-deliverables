#' Performs enrichment analysis
#'
#' This function loads the Rdata session produced using the mcl_clustering function:
#' There mcl clustering was performed and stored as vertex attribute 'community' of an igraph object.
#' For each community (or set of nodes), the function performs gene set enrichment analysis using gprofiler2 package.
#' Results are stored in a .csv file with columns cluster number, function, and gene set.
#' The detailed gprofiler results is stored as a list in a Rdata working session.
#'
#' @param directory Path to the .Rdata object
#' @return Stores the enrichment results in .csv and .Rdata objets
#' @export
run_gprofiler <- function(directory){
  mydir <- directory
  myfiles = list.files(mydir, pattern = "*\\.Rdata$", full.names = T)

  for (file in myfiles){
    filename <-  gsub("\\.Rdata$", "", basename(file))
    load(file)
    if (is.null(V(largest.conn.comp)$community)) next
    outDirectory <- file.path(mydir, paste0(filename,"_cluster_gprofiler.csv"))
    if (!file.exists(outDirectory)) {
      nodedata <- cbind(Genes= V(largest.conn.comp)$name, mclCluster= V(largest.conn.comp)$community)
      #nodedata <- as.data.frame(nodedata)
      clusters <- unique(nodedata[,"mclCluster"])
      list.genes.in.cluster <- lapply(clusters, function(cl){
        return(list(Name = paste0("Cluster ",cl),
                Genes = nodedata[which(cl==nodedata[,"mclCluster"]), "Genes"]))
      })
      names(list.genes.in.cluster) <- as.character(clusters)
      group.table <- NULL
      for (cl.net in clusters){
        gost.res <- gost(list.genes.in.cluster[[as.character(cl.net)]]$Genes,
                   organism = "hsapiens", ordered_query= FALSE,
                   correction_method="bonferroni",
                   significant=TRUE, user_threshold=0.05)
        group.table <- rbind(group.table,
                       c(cl.net, ifelse(is.null(gost.res$result$term_name),
                                        ".",
                                        paste(gost.res$result$term_name[1:min(length(gost.res$result$term_name), 5)],collapse=", ")
                       ),
                       paste0(list.genes.in.cluster[[as.character(cl.net)]]$Genes, collapse = ", ")))
        list.genes.in.cluster[[as.character(cl.net)]]$Function <- paste(gost.res$result$term_name[1:min(length(gost.res$result$term_name), 5)],collapse=", ")
        list.genes.in.cluster[[as.character(cl.net)]]$gprofiler <- gost.res$result
      }
      write.csv2(group.table, file = outDirectory, row.names = FALSE, quote = FALSE)
      save(list.genes.in.cluster, file= paste0(mydir,"/", filename,"_cluster_gprofiler.Rdata"))
    }
  }
}
