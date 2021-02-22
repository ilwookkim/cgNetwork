<!-- badges: start -->
  [![R-CMD-check](https://github.com/ilwookkim/cgNetwork/workflows/R-CMD-check/badge.svg)](https://github.com/ilwookkim/cgNetwork/actions)
  <!-- badges: end -->

# **cgNetwork**
Weighted gene co-expression network analysis by mutation status

## **Approach**
There are number of novel gene co-expression networks that regulates by transcriptional regulation genes. This package aims to analyze these networks of interesting gene that possibly interacts with and regulated by transcriptional regulation genes. In order to figure out these network, we analyze RNAseq data in dynamic gene status model such as cancers that has significant deregulation of these processes. 

  1. Retrive RNAseq data from various cancer database.
  1. RNAseq by **mutation_gene** mutation.
  1. Build network around **gene_of_interest**
  1. Comparison Network of **gene_of_interest** by mutation of **mutation_gene** 

### Installation

The **development** version can be installed from GitHub using:

``` r
devtools::install_github("ilwookkim/cgNetwork")
```

### Usage

``` r
library(cgNetwork)

gene_of_interest = "CDKN1A"
mutation_gene = "TP53"
```

**Acquire Data**

For the tutorial we only use a subset of genes ([Transcriptional Regulation by TP53 : Reactome R-HSA-3700989](https://reactome.org/content/detail/R-HSA-3700989)).
``` r 
library(fgsea)
gmt.file <- system.file("extdata", "ReactomePathways.gmt", package="cgNetwork")
TP53_pathway <- gmtPathways(gmt.file)[["Transcriptional Regulation by TP53"]]
```

Option 1: Using the [cgdsr package](https://cran.r-project.org/web/packages/cgdsr/index.html)
We need to create a cgds object, which also gives us the IDs of available projects. We choose "laml_tcga". Then we acquire information about this project via the cgInfo function, listing the available data types. We choose "mRNA expression (RNA Seq V2 RSEM)". Finally we get the data via the cgData function. We also get mutation data via the cgMutation function.
```r
cgds <- cgBase() #lists the available studies
studyID <- "laml_tcga"

cgInfo(cgds, studyID) #lists the available profiles and caselists
profile_name <- "mRNA expression (RNA Seq V2 RSEM)"

countdata <- cgData(cgds, studyID, profile_name, genes=TP53_pathway)
mut_df <- cgMutation(cgds, studyID, genes="TP53")
```

Option 2: Using the [TCGABiolinks package](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) (Approximately 1 GB of data will be downloaded)
We obtain the RNA-Seq data via the TCGA_RNAseq_RSEM function and mutation information via mutation_info function. Somatic Mutation Calling Workflow as tumor-normal pairs. Variant calling is performed using four separate pipelines (muse, varscan2, somaticsniper, mutect2), see details : [Somatic Mutation Calling Workflow](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#somatic-variant-calling-workflow).
```r
studyID = "LAML"
countdata <- TCGA_RNAseq_RSEM(studyID)
countdata <- countdata[rownames(countdata) %in% TP53_pathway,]
mut_df <- mutation_info(countdata, studyID, gene = mutation_gene, pipeline = "mutect2")
```

**Neighbor genes finder**

Finding neighbor genes network around the interesting gene using weighted network function from [igraph package](https://igraph.org/r/).

``` r
common_neighbor <- neighbor_finder(countdata, gene=gene_of_interest,
                                   cor.cut.off=.39, weight.cut.off=.5)                            
```

**Creating the Network**

Calculate networks (one for each mutation status) around the gene of interest. Interactive networks are vizualized via shiny.
``` r
cgNetwork_list <- cgNetwork(countdata, mut_df, common_neighbor)
DiNetplot(cgNetwork_list)
```

<img src="inst/extdata/DiNetwork_example.png"/>

Example shiny server: https://ilwookkim.shinyapps.io/dinetplot/


**Network data export and differential network visualization using Cytoscape**

  1. Download and Install Cytoscape.
    https://cytoscape.org/download.html
  1. Open Cytoscape
  1. Install plugin **Diffany**
  1. Run below in R
  ``` r
  g1 <- cgNetwork_list[[1]]
  g2 <- cgNetwork_list[[2]]
  RCy3::createNetworkFromIgraph(g1,"network_wt")
  RCy3::createNetworkFromIgraph(g2,"network_mut")
  ```
  1. Save network file as .cys.
  1. Import network_wt.cys and network_mut.cys files
  1. Apps/Diffany > Run Diffany project
  1. Go to Diffany tab
  1. Input networks window : include - check Two networks, Reference - check wt network
  1. Options window : Comparison mode - here I used One to all, Cutoff - Here I used 0.5 Check Differntial networks and Consensus networks.
  1. Press Start
  1. export image

<img src="inst/extdata/cytoscape_Diffany_example.png"/>

  Red: Up-regulate in Mut - Green: Down-regulate in Mut
