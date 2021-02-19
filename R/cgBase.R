#' create CGDS object
#'
#' @param address string, web address
#' @return cdgs object, and prints the first 2 columns of the cdgs object (the first column contains the IDs, to be used in later functions)
#' @details uses cgdsr::CGDS(address)
#' @examples
#' cgds <- cgBbase()
#' @export
#' @import cgdsr
cgBase <- function(address="https://www.cbioportal.org/"){
  #==cdgs object==#
  cgds <- CGDS("https://www.cbioportal.org/")
  #==user info print==#
  print( test(cgds) )
  print( getCancerStudies(cgds)[,1:2] )
  #==output==#
  cgds
}
