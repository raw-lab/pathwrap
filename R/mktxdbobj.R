#' Title
#'
#' @param geneAnnotation
#' @param corenum
#' @param genomeFile
#' @param entity
#'
#' @return
#' @export
#'
#' @examples
make_txdbobj <- function(geneAnnotation, corenum, genomeFile, entity){
txdb <- try(loadDb(geneAnnotation), silent = T)
cl2 <- makeCluster(corenum)
if (class(txdb)==  "TxDb"){
  if(!grepl("chr", seqlevels(txdb)[1])){ #check if this is necessary
    newSeqNames <- paste('Chr', seqlevels(txdb), sep = '')
    names(newSeqNames) <- seqlevels(txdb)
    txdb <- renameSeqlevels( txdb, newSeqNames )
    seqlevels(txdb)
  }
}else{
  library(GenomicFeatures)
  library(Rsamtools)
  chrLen <- Rsamtools::scanFaIndex(genomeFile)
  chrominfo <- data.frame(chrom = as.character(seqnames(chrLen)),
                          length = width(chrLen),
                          is_circular = rep(FALSE, length(chrLen)))
  txdb <- makeTxDbFromGFF(file = geneAnnotation, format = "gtf",
                          chrominfo = chrominfo,
                          dataSource = "Ensembl",
                          organism = entity)
}
return(txdb)
}
