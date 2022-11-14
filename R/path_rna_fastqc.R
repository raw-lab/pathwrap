
#' Title
#'
#' @param fq.dir
#' @param qc.dir
#' @param corenum
#'
#' @return
#' @export
#'
#' @examples
run_qc <- function(fq.dir, qc.dir, corenum){
  library(fastqcr)
  library(ggplot2)
  fastqcr::fastqc(fq.dir, qc.dir, threads =8,fastqc.path =system( "which fastqc", intern = T))
   qc <- qc_aggregate(qc.dir)
   tiff(file.path(qc.dir,"total_seq.tiff") ,units="in", width=15, height=15, res=300)
   ggplot(qc, aes(x=sample, y=tot.seq))+ geom_bar(stat="identity")  + coord_flip()
   dev.off()
   tiff(file.path(qc.dir,"qc_heatmap.tiff"), units="in", width=15, height=15, res=300)
   ggplot(qc, aes(x=sample, y =status)) + geom_tile()
   dev.off()
  print("Multiqc report done")
}
