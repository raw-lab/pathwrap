
#' Title
#' function to make good data frame for deseq2 and edgeR 
#' @param diff.tool
#' @param result.dir
#' @param grp.idx
#' @param geneLevels
#' @param entity
#' @param deseq2.dir
#'
#' @return
#' @export
#' 
#' 
#'
#' @examples
run_difftool <- function(diff.tool, result.dir,grp.idx, geneLevels, entity, deseq2.dir){ #make sure correct dir is chosen
  #geneData_my <- as.data.frame(readRDS(file.path(result.dir, "combinedcount.trimmed.RDS"))) #TO DO maybe I can just use geneLevels variable
  library(DESeq2)
  library(edgeR)
    #setwd("/scratch/edhungel/Rsubread/script")
  geneData_my <- geneLevels
  library(gage)
  library(pathview)
  cnts <- geneData_my[,-1]
  kegg.gs.species <- kegg.gsets(entity)
  orgcode<- kegg.species.code(entity)
  data(bods)
  #if(!all(rownames(cnts)%in% unlist(unname(kegg.gs.species$kg.sets)))){ #check if the use of "all" is appropriate
  if(sum(rownames(cnts)%in% unlist(unname(kegg.gs.species$kg.sets)) ) < 10){
    rownames(cnts)<- str_remove(rownames(cnts),"\\.[0-9]+$" )
    cnts<- mol.sum(cnts, id.map = "ENSEMBL", gene.annotpkg =bods[ which(bods[,3]==orgcode)]) #converting to entrez # what if gene id is not ensembl and what if arabidopsis thaliana id.map might be ath or else thing
  }

  #ref <- which(grp.idx == "reference")
  #samp <- which(grp.idx == "samples")
  #grp.idx[ref] <- "reference"
  #grp.idx[samp] <- "sample"

  if(diff.tool=="edgeR"){
   exp.fc <- run_edgeR(cnts, grp.idx, edgeR.dir)
  }
  else
  {
  exp.fc  <-run_deseq2(cnts,grp.idx, deseq2.dir)
  }

  ######
  #include bayseq

  return(c(exp.fc, cnts))
}
