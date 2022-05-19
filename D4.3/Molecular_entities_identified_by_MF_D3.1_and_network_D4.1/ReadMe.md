Crossing the results of matrix factorization and network-based analysis of multi-omics data

We report here the code and approach used for the analyses performed in WP4 D4.3. The aim is to investigate if the molecular entities defined in D3.1 and D4.1 are similar. For that, we compare the results of the applications of matrix factorization (D3.1) with the applications of network inference (D4.1) obtained on solid pediatric cancers studied in iPC.

As a reminder:
In D3.1, we performed unsupervised deconvolution using matrix factorization (more precisely stabilized Independent Component Analysis) of gene expression data from the 4 solid tumor types of interest in iPC. We used the results of the meta-analysis performed on gene expression of the solid tumor datasets, using the generated meta-weighted metagenes as defined in D3.1.

In D4.1, the aim was to build cancer type-specific multi-layered molecular and patient similarity networks. In the context of D4.3, we used the results obtained from the networks derived using Spearman correlation only, to be comparable with what was done in D3.1.

For this comparison, we investigate if some meta-weighted metagenes obtained in D3.1 are related to some subnetworks/modules defined in D4.1. More specifically, we test if the meta-weighted metagenes are enriched in the networks generated by omics layer for each dataset. To this end, the modules deciphered from network analyses in D4.1 were gathered in a gmt file.


Inputs:

The inputs needed to run this analysis are:

* the meta weighted metagenes, previously computed in D3.1 (we provide here a subset of them, due to size limitation)

* the community composition obtained after the computation of meta weighted metagenes, computed in D3.1

* the information about the datasets used in D3.1

* a gmt file containing the infered networks, obtained from the networks we derived using BMNPPSN R package

* a csv file containing the infered networks

They are provided here in the input_files folder.


Outputs:

The analysis provides different informations:

* a report "Report_*molecular_*entities_*identified_*by_*MF_*D3.1_*network_*D4.1.html" describing step by step how we selected the communities that we think to be of interest.(The folder ReportmolecularentitiesidentifiedbyMFD3.1networkD4.1_files contains all the files and images needed to build the html report)

* 2 tables describing these communities "community*stats.txt" and "community_complete_*description.txt"

* a table "tabletopcontributing_genes.txt" describing the 10 top contributing genes of these communities of interest

Usage recommendations:
The main document is Report_molecular_entities_identified_by_MF_D3.1_network_D4.1.Rmd, which calls template_enrichment_in_modules.Rmd and template_enrichment_in_GO.Rmd.
To launch this analysis, one has just to launch the main Rmd document, in a folder where the subfolder input_files is included.