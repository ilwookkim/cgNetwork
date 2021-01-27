#' TCGA_RNAseq_RSEM Function
#'
#' This function allows the user to download TCGA RNAseq data.
#' @param study_name name of TCGA study. Defaults to STAD (Stomach adenocarcinoma). Find more TCGA Study Abbreviations: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
#' @examples
#' countdata <- TCGA_RNAseq_RSEM(STAD)
#' countdata <- TCGA_RNAseq_RSEM(PAAD)
#' @export

TCGA_RNAseq_RSEM <- function(study_name="STAD"){
  library(TCGAbiolinks)
  query <- GDCquery(project = paste0("TCGA-",study_name),
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq",
                    file.type  = "results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  GDCdownload(query)
  expdat <- GDCprepare(query = query, save = FALSE)
  # Retrieve TP53_mutation status
  # CAUSION: It doesn't have all TCGA sample ID
  subtype <- TCGAquery_subtype(tumor = study_name)
  # Transform the data to get raw count table (RSEM)
  dataPrep <- TCGAanalyze_Preprocessing(object = expdat)
  countdata <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                         geneInfo = geneInfo,
                                         method = "geneLength")
  rm(dataPrep, query)
  # rename samples IDs of column to unifying samples IDs of subset
  # xx <- colnames(countdata)
  # xc <- gsub("-", ".", xx)
  # xz <- sub("^(.{12}).*", "\\1", xc)
  # colnames(countdata) <- xz
  # countdata <- as.data.frame(countdata)
  # wt.df <- data.frame(samples = subtype[subtype[, "TP53.mutation"] == "0",]$patient, condition = "WT")
  # mut.df <- data.frame(samples = subtype[subtype[, "TP53.mutation"] == "1",]$patient, condition = "MUT")
  # phenotypes <- rbind(wt.df, mut.df)
  # rm(wt.df, mut.df,subtype)
  # x <- phenotypes$samples
  # c <- gsub("-", ".", x)
  # phenotypes$samples <- c
  #
  # x <- data.frame(samples = colnames(countdata), condition = NA)
  # c <- merge(x=phenotypes, y= x, by = "samples")
  # c <- unique(c$samples)
  #
  # rownames(phenotypes) <- phenotypes$samples
  # phenotypes <- phenotypes[-1]
  # phenotypes <- data.frame( samples = c, condition = phenotypes[c,])
  # rownames(phenotypes) <- phenotypes$samples
  # phenotypes <- phenotypes[-1]
  # countdata <- countdata[,rownames(phenotypes)]
  # rm(x,c,xc,xx,xz)
  # And here we used only TP53 WT and Mut samples whcih is in total 275.
  # But count table has 450 samples, therefore we take only above 412 samples.
}
