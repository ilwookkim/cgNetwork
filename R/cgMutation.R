#' convenient shortcut to find mutations
#'
#' @param cgds object, output from the cgBase() function
#' @param cancer_study_id string, the as given in the cgBase() printout
#' @param genes genes string vector with genes of interest. Symbols, EntrezIds and EnsemblIds work. If nothing (NA) is provided, all genes in the GO database will be used. This will take some time, so if you only need a specific set of genes, specify it here.
#' @param caselist_name string, should be chosen so that only samples with mutation data are used.
#' @return data.frame, each row is a patient, each column is a gene
#' @examples
#' cgds <- cgBbase()
#' mystudy <- "laml_tcga"
#' myMutations <- cgData(cgds, mystudy, genes=c("FLT3","TP53"))
#' @export
#' @import cgdsr
cgMutation <- function(cgds, cancer_study_id, genes=NA, caselist_name="Samples with mutation data"){
  outtable <- cgData(cgds=cgds, cancer_study_id=cancer_study_id, profile_name="Mutations", caselist_name=caselist_name, genes=genes, dropNApatients=F)
  outtable <- t(outtable)
  outtable <- apply(outtable, 2, function(x) ifelse(is.na(x), 0, 1))
  colnames(outtable) <- paste0(colnames(outtable), "_mutation_status")
  #==output==#
  as.data.frame(outtable)
}
