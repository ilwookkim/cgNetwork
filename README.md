# **TCGANetwork**
TCGA gene co-expression network analysis

**Abstract**
There are number of gene networks that regulates by mutation by transcription target genes in various cancer. This package aims to analyze these networks under the certain mutation status.

**Approach**
In order to figure out these network, we analyze cancer RNAseq data from TCGA database. 
  1. Retrive RNAseq data from TCGA database.
  1. Seperate RNAseq data by mutation status of interesting gene.
  1. Check networks and neighbor genes under the mutation status

# Installation

The **development** version can be installed from GitHub using:

``` r
devtools::install_github("ilwookkim/TCGANetwork")
```

# Usage

``` r
library(TCGANetwork)
```

TCGA RNAseq data download

``` r
TCGA_study_name = "STAD"
countdata <- TCGA_RNAseq_RSEM(TCGA_study_name)
```

Mutation information

``` r
TCGA_study_name = "STAD"
gene = "TP53"
pipeline = "mutect2"
# There are four pipelines: muse, varscan2, somaticsniper, mutect2
mut_df <- mutation_info(countdata,TCGA_study_name, gene, pipeline)
```

Neighbor genes finder

``` r
common_neighbor <- neighbor_finder(countdata, gene="CDKN1A", cor_method = "spearman", cor.cut.off=.39, weight.cut.off=.2)
```

TCGA Network by mutation status of interesting gene

``` r
TCGANetwork_list <- TCGANetwork(countdata, mut_df, common_neighbor, cor_method = "spearman", weight.cut.off=.5)
```

Interactive clustered network plots by mutation status using shiny and visNetwork
``` r
Netplot(TCGANetwork_list)
```
