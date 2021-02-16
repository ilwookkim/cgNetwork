#' neighbor_finder
#'
#' Imports:
#' igraph,
#' bigmemory
#'
#' This function allows the user to find neighbor genes network of interested gene.
#' @param countdata RNAseq countdata
#' @param gene string, gene of interest to find a network for.
#' @param cor.cut.off cut off value of correlation coeffiecint to build network. Defaults to .39
#' @param weight.cut.off cut off value of network weight to reduce edges. Defaults to .5
#' @param t Matrix Transpose
#' @return character vector, genes that are considered neighbors of the gene of interest (including the gene of interest itself).
#' @examples
#' common_neighbor <- neighbor_finder(big_cor_matrix, gene="CDKN1A", cor_method = "spearman", cor.cut.off=.39, weight.cut.off=.5, t=TRUE)
#' @export

neighbor_finder <- function(countdata, gene=gene_of_interest, cor.cut.off=.39, weight.cut.off=.5, t=TRUE){
  if(!t){
    countdata_t <- data.frame(countdata)
  } else {
    countdata <- data.frame(na.omit(countdata))
    countdata_t <- data.frame(t(countdata))
  }
  countdata_t <- countdata_t[, !sapply(countdata_t, function(x) { sd(x) == 0} )]
  big_cor_matrix <- bigmemory::as.big.matrix(cor(countdata_t, method = "spearman"))

  #check if the gene is at all present. Otherwise there would be a calculation first (which might take long) and then the error.
  if(!any(rownames(big_cor_matrix[,]) %in% gene)) stop("Your gene of interest is not in the matrix.")

  big_cor_matrix[which(big_cor_matrix[,] == 1)] <- 0
  big_cor_matrix[which(big_cor_matrix[,] < cor.cut.off)] <- 0
  net <- igraph::graph_from_adjacency_matrix(big_cor_matrix[,], mode='undirected', weighted = T, diag=F)
  net <- igraph::simplify(net, remove.multiple = T, remove.loops = T)
  net.sp <- igraph::delete_edges(net, igraph::E(net)[weight<weight.cut.off])
  neigh.nodes <- igraph::neighbors(net.sp, igraph::V(net.sp)[gene])
  common_neig <- c(names(neigh.nodes),gene)
  return(common_neig)
}


