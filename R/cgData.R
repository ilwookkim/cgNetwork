#' output data of the specified type
#'
#' Imports:
#' cgdsr,
#' org.Hs.eg.db
#'
#' @param cgds object, output from the cgBase() function
#' @param cancer_study_id string, the ID as given in the first column of the cgBase() printout
#' @param profile_name string, profile name (not ID) to be used (obtained from the printout of the cgStudy() function)
#' @param caselist_name string, caselist name (not ID) to be used (obtained from the printout of the cgStudy() function)
#' @param genes string vector with genes of interest. Symbols, EntrezIds and EnsemblIds work. If nothing (NA) is provided, all genes in the GO database will be used. This will take some time, so if you only need a specific set of genes, specify it here.
#' @return data.frame, each row is a patient, each column is a gene
#' @examples
#' cgds <- cgBbase()
#' mystudy <- "laml_tcga"
#' cgStudy(cgds, mystudy)
#' myData <- cgData(cgds, mystudy, profile_name="mRNA expression (RNA Seq RPKM)", genes=c("FLT3","TP53"))
#' @export
cgData <- function(cgds, cancer_study_id, profile_name, caselist_name="All samples", genes=NA){
  #if no genes are specified, get all genes from GO
  if(is.na(genes)[1]) genes <- as.vector(unique(unlist(as.list(org.Hs.eg.db::org.Hs.egGO2EG))))
  nmath <- floor(length(genes)/400)
  nmath <- seq(0,nmath)*400

  #==get profile and caselist ID (type of experiment)==#
  profileTable  <- cgdsr::getGeneticProfiles(cgds, cancer_study_id)
  caselistTable <- cgdsr::getCaseLists(      cgds, cancer_study_id)
  profileId  <- subset(profileTable, genetic_profile_name %in% profile_name )$genetic_profile_id
  caselistId <- subset(caselistTable,      case_list_name %in% caselist_name)$case_list_id
  print( paste("profile ID:",  profileId ) )
  print( paste("caselist ID:", caselistId) )
  #==get data for the genes of interest==#
  data1 <- lapply(nmath, function(x) {
    if(x%%4000==0) print(paste("getting genes",x+1,"to",min(x+4000,length(genes))))
    genes400 <- as.vector(na.omit(genes[seq(x+1,x+400)]))
    cgtable <- cgdsr::getProfileData(x=cgds, genes=genes400, geneticProfiles=profileId, caseList=caselistId)
  })
  data1 <- Reduce(function(x,y) transform(merge(x,y,by=0), row.names=Row.names, Row.names=NULL), data1) #successively merge all tables
  data1[,1:ncol(data1)] <- apply(data1, 2, function(x) ifelse(x %in% "NaN",NA,x)) #change NaNs to NAs
  #==output==#
  data1 <- as.data.frame(data1)
  data1 <- data1[apply(data1,1,function(x) !all(is.na(x))),]
  as.data.frame(t(data1))
}
