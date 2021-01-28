# **DDNetwork**
TCGA gene co-expression network analysis

**Abstract**
There are number of gene networks that regulates by mutation by transcription target genes in various cancer. This package aims to analyze these networks under the certain mutation status.

**Approach**
In order to figure out these network, we analyze cancer RNAseq data from TCGA database. 
  1. Retrive RNAseq data from TCGA database.
  1. Seperate RNAseq data by mutation status of interesting gene.
  1. Check networks and neighbor genes under the mutation status

## **Functions**
  **1.TCGA_RNAseq_RSEM(TCGA_study_name)**
    
  **1.mutation_info(TCGA_study_name, gene, pipeline)**
