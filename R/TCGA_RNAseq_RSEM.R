#' TCGA_RNAseq_RSEM Function
#'
#' @import TCGAbiolinks
#'
#' This function allows the user to download TCGA RNAseq data.
#' @param study_name name of TCGA study. Defaults to STAD (Stomach adenocarcinoma). Find more TCGA Study Abbreviations: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
#' @examples
#' countdata <- TCGA_RNAseq_RSEM("STAD")
#' countdata <- TCGA_RNAseq_RSEM("PAAD")
#' @export

TCGA_RNAseq_RSEM <- function(study_name="STAD"){
  query <- GDCquery(project = paste0("TCGA-",study_name),
                                  data.category = "Gene expression",
                                  data.type = "Gene expression quantification",
                                  platform = "Illumina HiSeq",
                                  file.type  = "results",
                                  experimental.strategy = "RNA-Seq",
                                  legacy = TRUE)
  GDCdownload(query, files.per.chunk = 10)
  expdat <- GDCprepare(query = query, save = FALSE)
  # Transform the data to get raw count table (RSEM)
  dataPrep <- TCGAanalyze_Preprocessing(object = expdat)
  countdata <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                                       geneInfo = geneInfo,
                                                       method = "geneLength")
  samples <- colnames(countdata)
  samples <- sub("^(.{15}).*", "\\1", samples)
  colnames(countdata)<- samples
  rm(dataPrep, query)
  return(countdata)
}

