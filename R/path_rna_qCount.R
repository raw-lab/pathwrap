

#' Title
#'
#' @param genomeFile
#' @param geneAnnotation
#' @param aligned_proj
#' @param corenum
#' @param result.dir
#' @param txdb
#'
#' @return
#' @export
#'
#' @examples
run_qCount <- function(genomeFile, geneAnnotation, aligned_proj,corenum, result.dir, txdb){
  ##
  #for mapping
  library(Rsamtools) #scanFaIndex
  library(GenomicFeatures)
  library(parallel)
  cl2 <- makeCluster(corenum)
  geneLevels <- QuasR::qCount(aligned_proj, txdb, reportLevel ="gene", clObj=cl2)
  saveRDS(geneLevels, file.path(result.dir, "combinedcount.trimmed.RDS"))
  return(geneLevels)
}
