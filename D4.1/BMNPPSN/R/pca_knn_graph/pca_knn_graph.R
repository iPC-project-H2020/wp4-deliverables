### MEDULLOBLASTOMA CAVALLI
mycols <- c("red", "green", "gold", "lightblue")
names(mycols) <- unique(cavalli_sample_description$Subgroup)

cols.subtypes <- c("brown4","lightgreen","darkgoldenrod1","darkgoldenrod4","yellow","brown1","darkolivegreen", 
                   "brown3", "cyan", "blue", "red", "aquamarine4")
names(cols.subtypes) <- unique(cavalli_sample_description$Subtype_formatted)

### Select a matrix here 
sel.matrix <- summa_eigenGenes_concatenated_cavalli
### Select an output directory
outDirectory <- getwd()

## PCA 
library("factoextra")
pca.res <- prcomp(sel.matrix)
fviz_pca_ind(pca.res, repel=TRUE, col.ind = mycols[cavalli_sample_description$Subgroup], 
             legend.title="Subtypes",
             habillage=cavalli_sample_description$Subgroup)
fviz_pca_var(pca.res, repel = TRUE)
fviz_pca_biplot(pca.res, label = 'none')

## construct knn graph and plot
library(cccd)
X11()
par(mfrow=c(3,3))
par(mar=c(0,0,0,0)+2)
dm <- as.matrix(dist(sel.matrix))       # full distance matrix
k.set <- seq(10, 90, by = 10)
for (sel.k in k.set){
  if (sel.k==min(k.set)){ title <- "Concatenated \n kNN k="
  }else{title <- "k="}
  knn.mat <- nng(x = sel.matrix, k=sel.k, method = "euclidean",
                 mutual = TRUE)
  V(knn.mat)$name <- rownames(sel.matrix)
  V(knn.mat)$color.groups <- mycols[cavalli_sample_description$Subgroup]
  V(knn.mat)$color.subtypes <- cols.subtypes[cavalli_sample_description$Subtype_formatted]
  l <- layout_with_fr(knn.mat)
  plot(knn.mat,
       main=paste(title, sel.k, sep=""), 
       vertex.label=NA, layout=l, 
       vertex.size=3, vertex.color=V(knn.mat)$color.subtypes)
  distance <- apply(get.edgelist(knn.mat,1:ecount(knn.mat)),1,function(x)dm[x[1],x[2]])
  compg.edges <- as.data.frame(get.edgelist(knn.mat))
  pp.net <- cbind(compg.edges, distance)
  colnames(pp.net) <- c("e1", "e2", "distance")
  write.csv.gz(pp.net, file= paste0(outDirectory, "knn_k", sel.k, ".csv.gz"), row.names = TRUE)
}


### MEDULLOBLASTOMA FORGET
## Select an eigenGene matrix
sel.matrix <- summa_eigenGenes_concatenated_forget
## Select a directory
outDirectory <- getwd()
  
## PCA 
library("factoextra")
mycols <- c("red", "green", "gold", "lightblue")
names(mycols) <- unique(forget_sample_description)

pca.res <- prcomp(sel.matrix)
fviz_pca_ind(pca.res, repel=TRUE, col.ind = mycols[forget_sample_description], legend.title="Subtypes",
             habillage=forget_sample_description)
## knn graph
library(cccd)
X11()
par(mfrow=c(3,3))
par(mar=c(0,0,0,0)+2)
dm <- as.matrix(dist(sel.matrix))       # full distance matrix
k.set <- seq(3, 12, 1)
for (sel.k in k.set){
  if (sel.k==min(k.set)){ title <- "Concatenated \n kNN k="
  }else{title <- "k="}
  knn.mat <- nng(x = sel.matrix, k=sel.k, method = "euclidean",
                 mutual = TRUE)
  V(knn.mat)$name <- rownames(sel.matrix)
  V(knn.mat)$color <- mycols[forget_sample_description]
  l <- layout_with_fr(knn.mat)
  plot(knn.mat,
       main=paste(title, sel.k, sep=""), 
       vertex.label=NA, layout=l, 
       vertex.size=8)
  distance <- apply(get.edgelist(knn.mat,1:ecount(knn.mat)),1,function(x)dm[x[1],x[2]])
  compg.edges <- as.data.frame(get.edgelist(knn.mat))
  pp.net <- cbind(compg.edges, distance)
  colnames(pp.net) <- c("e1", "e2", "distance")
  write.csv.gz(pp.net, file= paste0(outDirectory, "knn_k", sel.k, ".csv.gz"), row.names = TRUE)
}



### NEUROBLASTOMA or EWING SARCOMA (just select the matrix you prefer)
sel.matrix <- summa_eigenGenes_concatenated_henrich
outDirectory <- getwd()

## PCA 
library("factoextra")
pca.res <- prcomp(sel.matrix)
fviz_pca_ind(pca.res, repel=TRUE)
fviz_pca_var(pca.res, repel=TRUE)
## knn graph
library(cccd)
X11()
par(mfrow=c(3,3))
par(mar=c(0,0,0,0)+2)
dm <- as.matrix(dist(sel.matrix))       # full distance matrix
k.set <- seq(5, 50, 5)
for (sel.k in seq(5, 50, 5)){
  knn.mat <- nng(x = sel.matrix, k=sel.k, method = "euclidean",
                 mutual = TRUE)
  V(knn.mat)$name <- rownames(sel.matrix)
  l <- layout_with_fr(knn.mat)
  if (sel.k==min(k.set)){ title <- "Concatenated \n kNN k="
  }else{title <- "k="}
  plot(knn.mat,
       main=paste(title, sel.k, sep=""), 
       vertex.label=NA, layout=l, 
       vertex.size=5)
  distance <- apply(get.edgelist(knn.mat,1:ecount(knn.mat)),1,function(x)dm[x[1],x[2]])
  compg.edges <- as.data.frame(get.edgelist(knn.mat))
  pp.net <- cbind(compg.edges, distance)
  colnames(pp.net) <- c("e1", "e2", "distance")
  write.csv.gz(pp.net, file= paste0(outDirectory, "knn_k", sel.k, ".csv.gz"), row.names = TRUE)
}



