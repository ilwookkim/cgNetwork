#' mutation_info Function
#'
#' This function allows the user to download mutation infomation.
#' @param study_name name of TCGA study. Defaults to "STAD" (Stomach adenocarcinoma). Find more TCGA Study Abbreviations: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
#' @param gene Hugo Symbol of gene Default to "TP53"
#' @param pipeline Variant calling is performed using several separate pipelines; muse, varscan2, somaticsniper, mutect2. Default to mutect2. Find details: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#somatic-variant-calling-workflow
#' @examples
#' mut_info <- mutation_info("STAD", "TP53", "mutect2")
#' mut_info <- mutation_info("PAAD", "KRAS","muse")
#' @export

mutation_info <- function(study_name="STAD", gene="TP53", pipeline = "metect2"){
  maf <- TCGAbiolinks::GDCquery_Maf(study_name, pipelines = pipeline)
  sel_maf <- maf[which(maf$Hugo_Symbol == gene),]
  xx <- sel_maf$Tumor_Sample_Barcode
  #xc <- gsub("-", ".", xx)
  xz <- sub("^(.{12}).*", "\\1", xx)
  mut_samples <- unique(xz)
  subtype <- TCGAbiolinks::TCGAquery_subtype(tumor = study_name)
  mut_df <- data.frame(row.names = subtype$patient)
  mut_df$mut_status <- 0
  mut_df$mut_status[which(rownames(mut_df) %in% mut_samples)] <- 1
  names(mut_df) <- paste0(gene,"_mut_status")
  return(mut_df)
}
