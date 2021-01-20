-   [Intalling BMNPPSN](#intalling-bmnppsn)
-   [Using BMNPPSN](#using-bmnppsn)
-   [Create a subnetwork](#create-a-subnetwork)
-   [Tunning the inflation parameter for clustering and
    clustering](#tunning-the-inflation-parameter-for-clustering-and-clustering)
-   [Enrichment analysis of the communities identified through
    MCL](#enrichment-analysis-of-the-communities-identified-through-mcl)
-   [Construction of eigenGene matrix using
    Roma](#construction-of-eigengene-matrix-using-roma)
-   [Knn patient similarity networks](#knn-patient-similarity-networks)
-   [Knn patient for WP4.1 multi-omics
    datasets](#knn-patient-for-wp4.1-multi-omics-datasets)

This package provides the codes to Build Molecular Networks and
Patient-Patient Similarity Networks (BMNPPSN) described in details in
deliverable 4.1. The package is currently being tested.

Intalling BMNPPSN
=================

The BMNPPSN package relies on R version 4.0.3, thus BioConductor must be
updated to version 3.12.

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(version = "3.12")

    ## Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installation path not writeable, unable to update packages: codetools, foreign,
    ##   KernSmooth, Matrix, nlme

    ## Old packages: 'backports', 'BH', 'caTools', 'cpp11', 'crosstalk', 'data.table',
    ##   'DBI', 'DelayedMatrixStats', 'deldir', 'diffobj', 'dplyr', 'DT', 'expm',
    ##   'fansi', 'gdtools', 'ggplot2', 'ggrepel', 'gmp', 'hexbin', 'hms',
    ##   'htmltools', 'httpcache', 'httpuv', 'parallelly', 'plotly', 'preprocessCore',
    ##   'prettydoc', 'pROC', 'quantreg', 'ragg', 'Rcpp', 'RcppArmadillo', 'rgl',
    ##   'rlang', 'RMariaDB', 'robustbase', 'RSQLite', 'rtracklayer', 'scRNAseq',
    ##   'sf', 'sp', 'spatstat.utils', 'speedglm', 'tibble',
    ##   'TxDb.Hsapiens.UCSC.hg38.knownGene', 'VGAM', 'VIM', 'withr', 'xfun'

In addition, BMNPPSN package relies on package `gprofiler2` that can be
installed using `BiocManager`, and `rROMA` that can be installed using
`devtools`:

    if(!"gprofiler2" %in% installed.packages()){
      BiocManager::install("gprofiler2")
    }
    if(!requireNamespace("devtools")){
      install.packages("devtools")
    }
    if(!requireNamespace("rRoma")){
      devtools::install_github("Albluca/rROMA")
    }

Similarly, it is possible to install the BMNPPS package using devtools
by typing:

    if(!requireNamespace("BMNPPSN")){
      devtools::install_github("iPC-project-H2020/wp4-deliverables/BMNPPSN")
    }

Using BMNPPSN
=============

Upload the BMNPPSN package:

    library(BMNPPSN)

    ## Loading required package: tictoc

    ## Loading required package: crunch

    ## 
    ## Attaching package: 'crunch'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, rstandard, setNames

    ## The following object is masked from 'package:graphics':
    ## 
    ##     title

    ## The following object is masked from 'package:utils':
    ## 
    ##     write.csv

    ## The following object is masked from 'package:base':
    ## 
    ##     table

    ## Loading required package: data.table

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:crunch':
    ## 
    ##     copy, cube, rollup

    ## Loading required package: igraph

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    ## Loading required package: readr

    ## Loading required package: gprofiler2

    ## Loading required package: MCL

    ## Loading required package: Matrix

    ## Loading required package: scater

    ## Loading required package: SingleCellExperiment

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:igraph':
    ## 
    ##     normalize, path, union

    ## The following objects are masked from 'package:crunch':
    ## 
    ##     combine, duplicated, sd, table, type, type<-

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:Matrix':
    ## 
    ##     expand

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     first, second

    ## The following objects are masked from 'package:crunch':
    ## 
    ##     values, values<-

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:data.table':
    ## 
    ##     shift

    ## The following objects are masked from 'package:crunch':
    ## 
    ##     members, trim

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## The following objects are masked from 'package:crunch':
    ## 
    ##     description, description<-, notes, notes<-

    ## Loading required package: ggplot2

    ## 
    ## Attaching package: 'ggplot2'

    ## The following object is masked from 'package:crunch':
    ## 
    ##     resolution

    ## Loading required package: devtools

    ## Loading required package: usethis

    ## Loading required package: BiocManager

    ## Bioconductor version 3.12 (BiocManager 1.30.10), ?BiocManager::install for help

    ## 
    ## Attaching package: 'BiocManager'

    ## The following object is masked from 'package:devtools':
    ## 
    ##     install

    ## Loading required package: cccd

BMNPPS requires an expression matrix - with column names indicating
samples and row names indicating gene names - and a network file
produced by COSIFER - with columns e1, e2 and intensity, in which each
row indicates the intensity of the edge between node e1 and node e2.

Create a subnetwork
===================

It is possible to create a subnetwork simply by typing
`create_subnetwork(directory, outDirectory, n.edges=10^5)` where
`directory` is the directory of the original network in .csv.gz format,
and `outDirectory` the user defined output directory.

Tunning the inflation parameter for clustering and clustering
=============================================================

Tunning and clustering is as simple as typing
`mcl_clustering(resDirectory, outDirectory)` where now `resDirectory` is
the directory where subnetworks have been saved and `outDirectory` the
user defined output directory. MCL clustering is performed with
inflation parameters from 1 to 3.5 with a step of 0.1. The inflation for
which the max number of clusters with more than 10 genes is attaint, is
retained and stored.

Enrichment analysis of the communities identified through MCL
=============================================================

Now that we have identify communities in the network,
`run_gprofiler(directory)` performs gene enrichment analysis of the
communities (saved in `directory`) and stores the results in
`directory`.

Construction of eigenGene matrix using Roma
===========================================

To construct the eigenGene matrix we need now an input `directory` where
the expression matrix is stored and the `outDirectory` containing the
results from the previous analysis. Thus,
`runRoma(directory, outDirectory)`computes the eigenGene matrix and
stores it in the `outDirectory`.

Knn patient similarity networks
===============================

Once we have the eigenGene matrix, construction of knn graph can be
performed just by
typing:`knn.mat <- nng(x = eigenGene.matrix, k=sel.k, method = "euclidean", mutual = TRUE)`,
where `sel.k` is the user selected k in the knn graph.

Knn patient for WP4.1 multi-omics datasets
==========================================

The construction of knn graph for the multiomics datasets considered in
WP4.1 is made available in folder `R\pca_knn_graph.R` and the eigenGene
matrices in `data`
