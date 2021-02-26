#' cgNetwork
#'
#' Find gene networks around gene of interest, one for WT and one for mutated.
#'
#' @param countdata, output of the cgData or TCGA_RNAseq_RSEM function.
#' @param mut_df data.frame, output of the cgMutation or mutation_info function.
#' @param common_neighbor character vector, names of genes considered neighbors of the gene of interest. Output of the neighbor_finder function.
#' @param cor_method string, correlation coeffiecint method "spearman" or "pearson". Defaults to "spearman".
#' @param weight.cut.off numeric, cut off value of network weight to reduce edges. Defaults to 0.5.
#' @return list of two igraph objects
#' @examples
#' example.file <- system.file("extdata", "example.RData", package="cgNetwork")
#' load(example.file)
#' cgNetwork_list <- cgNetwork(countdata, mut_df, common_neighbor, "spearman", 0.5)
#' @export
#' @import dplyr
#' @importFrom igraph graph.adjacency delete.edges V E simplify delete.vertices degree
#' @importFrom stats cor na.omit sd

cgNetwork <- function(countdata, mut_df, common_neighbor, cor_method = "spearman", weight.cut.off=.5){
  # make sure that the mut_df only has genes that occur in the countdata
  mut_df <- subset(mut_df, rownames(mut_df) %in% colnames(countdata))

  wt  <- rownames(mut_df)[which(mut_df[,1]==0)]
  mut <- rownames(mut_df)[which(mut_df[,1]==1)]

  na_omit_df <- na.omit(countdata)
  rnaseq <- data.frame(t(na_omit_df))

  network <- list(
    wtNetwork  = .subNetwork(rnaseq=rnaseq, mut_status=wt,  common_neighbor=common_neighbor, cor_method=cor_method, weight.cut.off=weight.cut.off),
    mutNetwork = .subNetwork(rnaseq=rnaseq, mut_status=mut, common_neighbor=common_neighbor, cor_method=cor_method, weight.cut.off=weight.cut.off)
  )

  return(network)
}

#function to create the igraph network
.subNetwork <- function(rnaseq, mut_status, common_neighbor, cor_method, weight.cut.off){
  rnaseq <- rnaseq[mut_status, common_neighbor]
  rnaseq <- rnaseq[, !sapply(rnaseq, function(x) { sd(x) == 0} )]
  rnaseq_cor <- cor(rnaseq, method = cor_method)
  correl <- rnaseq_cor
  g1 <- graph.adjacency(correl, weighted=TRUE, mode="lower")
  g1 <- delete.edges(g1, E(g1)[ weight < weight.cut.off ])
  g1 <- simplify(g1, remove.multiple = TRUE, remove.loops = TRUE)
  g1 <- delete.vertices(g1, which(degree(g1)<1))
  return(g1)
}
