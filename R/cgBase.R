#' create CGDS object
#'
#' Imports:
#' cgdsr
#'
#' @param address string, web address
#' @return cdgs object, and prints the first 2 columns of the cdgs object (the first column contains the IDs, to be used in later functions)
#' @details uses cgdsr::CGDS(address)
#' @examples
#' cgds <- cgBbase()
#' @export
cgBase <- function(address="https://www.cbioportal.org/"){
  #==cdgs object==#
  cgds <- cgdsr::CGDS("https://www.cbioportal.org/")
  #==user info print==#
  print( cgdsr::test(cgds) )
  print( cgdsr::getCancerStudies(cgds)[,1:2] )
  #==output==#
  cgds
}
