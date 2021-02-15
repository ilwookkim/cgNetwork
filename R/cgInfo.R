#' get info about a study
#'
#' Imports:
#' cgdsr
#'
#' @param cgds object, output from the cgBase() function
#' @param cancer_study_id string, the ID as given in the first column of the cgBase() printout
#' @return nothing, but prints out two tables with information for the follow-up function
#' @details uses cgdsr::getCancerStudies(cgds)
#' @examples
#' cgds <- cgBbase()
#' mystudy <- "laml_tcga"
#' cgInfo(cgds, mystudy)
#' @export
cgInfo <- function(cgds, cancer_study_id){
  #==user info print==#
  cat("\n## Profiles:\n")
  print(cgdsr::getGeneticProfiles(cgds, cancer_study_id)["genetic_profile_name"])
  cat("## Samples:\n")
  print(cgdsr::getCaseLists      (cgds, cancer_study_id)["case_list_name"])
}
