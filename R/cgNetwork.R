#' cgNetwork
#'
#' Imports:
#' igraph
#'
#' Find gene networks around gene of interest, one for WT and one for mutated.
#' 
#' @param data.frame, output of the cgData or TCGA_RNAseq_RSEM function.
#' @param mut_df data.frame, output of the cgMutation or mutation_info function.
#' @param common_neighbor character vector, names of genes considered neighbors of the gene of interest. Output of the neighbor_finder function.
#' @param cor_method string, correlation coeffiecint method "spearman" or "pearson". Defaults to "spearman".
#' @param weight.cut.off numeric, cut off value of network weight to reduce edges. Defaults to 0.5.
#' @return list of two igraph objects
#' @examples
#' cgNetwork_list <- cgNetwork(countdata, mut_df, common_neighbor, cor_method = "spearman", weight.cut.off=.5)
#' @export

cgNetwork <- function(countdata, mut_df, common_neighbor, cor_method = "spearman", weight.cut.off=.5){
  # make sure that the mut_df only has genes that occur in the countdata
  mut_df <- subset(mut_df, rownames(mut_df) %in% colnames(countdata))
  
  wt  <- rownames(mut_df)[which(mut_df[,1]==0)]
  mut <- rownames(mut_df)[which(mut_df[,1]==1)]

  na_omit_df <- na.omit(countdata)
  rnaseq <- data.frame(t(na_omit_df))

  network <- list(
    wtNetwork  = .subNetwork(rnaseq=rnaseq, mut_status=wt,  cor_method=cor_method),
    mutNetwork = .subNetwork(rnaseq=rnaseq, mut_status=mut, cor_method=cor_method)
  )
  
  return(network)
}

#function to create the igraph network
.subNetwork <- function(rnaseq, mut_status, cor_method){
  rnaseq_1 <- rnaseq[mut_status, common_neighbor]
  rnaseq_1 <- rnaseq[, !sapply(rnaseq_1, function(x) { sd(x) == 0} )]
  rnaseq_cor <- cor(rnaseq_1, method = cor_method)
  cor_1 <- rnaseq_cor
  rm(rnaseq_1, rnaseq_cor)
  g1 <-  igraph::graph.adjacency(cor_1, weighted=TRUE, mode="lower")
  g1 <- igraph::delete.edges(g1, igraph::E(g1)[ weight < weight.cut.off ])
  g1 <- igraph::simplify(g1, remove.multiple = TRUE, remove.loops = TRUE)
  g1 <- igraph::delete.vertices(g1, which(igraph::degree(g1)<1))
  return(g1)
}
