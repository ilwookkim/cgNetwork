#' neighbor_finder
#'
#' This function allows the user to find neighbor genes network of interested gene.
#' @param big_cor_matrix bigmatrix (bigmemory package) of Spearman's correlation coefficient from TCGA RNAseq countdata
#' @param gene Interesting gene to find network. Defaults to "CDKN1A".
#' @param cor.cut.off cut off value of correlation coeffiecint to build network. Defaults to .39
#' @param weight.cut.off cut off value of network weight to reduce edges. Defaults to .5
#' @examples
#' common_neighbor <- neighbor_finder(big_cor_matrix, gene="CDKN1A", cor_method = "spearman", cor.cut.off=.39, weight.cut.off=.5)
#' @export

neighbor_finder <- function(big_cor_matrix, gene="CDKN1A", cor.cut.off=.39, weight.cut.off=.5){
  big_cor_matrix[which(big_cor_matrix[,] == 1)] <- 0
  big_cor_matrix[which(abs(big_cor_matrix[,]) < cor.cut.off)] <- 0
  net <- igraph::graph_from_adjacency_matrix(big_cor_matrix[,], mode='undirected', weighted = T, diag=F)
  net <- igraph::simplify(net, remove.multiple = T, remove.loops = T)
  net.sp <- igraph::delete_edges(net, igraph::E(net)[weight<weight.cut.off])
  neigh.nodes <- igraph::neighbors(net.sp, igraph::V(net.sp)[gene])
  common_neig <- c(names(neigh.nodes),gene)
  return(common_neig)
}

