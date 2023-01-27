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
run_qCount <- function(genomeFile, geneAnnotation, aligned_proj,corenum, result.dir, txdb, entity, coldata){
  ##
  #for mapping
  library(Rsamtools) #scanFaIndex
  library(GenomicFeatures)
  library(parallel)
  cl2 <- makeCluster(corenum)
  geneLevels <- QuasR::qCount(aligned_proj, txdb, reportLevel ="gene", clObj=cl2)

  library(DESeq2)
  library(edgeR)
  #####################
  #post processing for count
  library(gage)
  library(pathview)
  cnts <- geneLevels[,-1]
  kegg.gs.species <- kegg.gsets(entity)
  orgcode<- kegg.species.code(entity)
  data(bods)
  #if(!all(rownames(cnts)%in% unlist(unname(kegg.gs.species$kg.sets)))){ #check if the use of "all" is appropriate
  if(sum(rownames(cnts)%in% unlist(unname(kegg.gs.species$kg.sets)) ) < 10){
    rownames(cnts)<- str_remove(rownames(cnts),"\\.[0-9]+$" )
    cnts<- mol.sum(cnts, id.map = "ENSEMBL", gene.annotpkg =bods[ which(bods[,3]==orgcode)]) #converting to entrez # what if gene id is not ensembl and what if arabidopsis thaliana id.map might be ath or else thing
  }
  
  cnts <- cnts[, rownames(coldata)]
  saveRDS(cnts, file.path(result.dir, "combinedcount.trimmed.RDS"))
  
  if(  all(rownames(coldata) == colnames(cnts)) ){#if this then proceed
    ref <- which(coldata[, 2] ==  levels(coldata[, 2])[1])
    samp <- which(coldata[, 2] ==  levels(coldata[, 2])[2])
    grp.idx <-NULL
    grp.idx[ref] <- "reference"
    grp.idx[samp] <- "sample"
  }
  else{
    message("make sure pheno file have only samples analysed")
  }
  #print("this is the counts data")
 #print( head(cnts))
  countreturnlist <- list("cnts" = cnts , "grp.idx"= grp.idx)
  return(countreturnlist )
}
