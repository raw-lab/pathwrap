#############################################################################
#2. RUN THE ANALYSIS # Alignmnet and counting
#############################################################################

#' Title
#'
#' @param corenum
#' @param endness
#' @param sampleFile
#' @param genomeFile
#' @param geneAnnotation
#' @param ref.dir
#'
#' @importFrom QuasR qAlign
#'
#' @return
#' @export
#'
#' @examples
run_qAlign <- function(corenum, endness, sampleFile, genomeFile,geneAnnotation, ref.dir,cacheDir){
  #does ref.dir also have ref index, if not make indexes
  if( !is.na(ref.dir)){
      if( length(list.files(ref.dir , ".Rhisat2$", full.names = T)) !=1){
    sampleFiletmp <- read.table(sampleFile, "\t", header = T)[1,]
    sampleFiletmp_name <- paste0(gsub("sampleFile.txt","sampleFiletmp.txt", sampleFile))
    write.table(sampleFiletmp, sep=  "\t", col.names = T, row.names = F, file = sampleFiletmp_name)
    cl2 <- makeCluster(corenum)
    if(endness == "SE"){
      aligned_proj <-  QuasR::qAlign(sampleFiletmp_name, paired ="no", clObj=cl2, alignmentsDir =aligned_bam ,
                                     genome=genomeFile,geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner ="Rhisat2",cacheDir )
    }else{
      aligned_proj <-  QuasR::qAlign(sampleFiletmp_name, paired ="fr", clObj=cl2, alignmentsDir =aligned_bam ,
                                     genome=genomeFile,geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner ="Rhisat2" ,cacheDir)
    }# this will form the reference index

    #the program check for aligned bam before running so we dont really need to remove this sample from our sampleFile
    unlink(sampleFiletmp_name)
    print("We made tmp file, and made index and one alignment")
      }
}
  cl2 <- makeCluster(corenum)
  print("Alignment is running")
  if (endness=="PE"){
    aligned_proj <- QuasR::qAlign(sampleFile,paired ="fr", clObj=cl2, alignmentsDir =aligned_bam , genome=genomeFile, geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner ="Rhisat2" ,cacheDir)
  } else {

    aligned_proj <- QuasR::qAlign(sampleFile,paired ="no", clObj=cl2, alignmentsDir =aligned_bam ,  genome=genomeFile,geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner ="Rhisat2",cacheDir )
  print("done")
    }
  saveRDS(aligned_proj, file.path(aligned_bam , "alltrimmedalignedobj.RDS"))
  stopCluster(cl2)
  return(aligned_proj)
}
