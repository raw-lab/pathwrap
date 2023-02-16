
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
  on.exit(closeAllConnections())
library(fastqcr)
library(ggplot2)
#install fastqc if system( "which fastqc", intern = T) fails
if (Sys.which("fastqc")=="" & !file.exists(paste0(qc.dir, "/FastQC/fastqc"))){
  ## work here
  fastqcr::fastqc_install(dest.dir = qc.dir)
  fastqc.path <- paste0(qc.dir, "/FastQC/fastqc")
}
else{
  fastqc.path <- Sys.which("fastqc")
}

fastqcr::fastqc(fq.dir, qc.dir,fastqc.path =fastqc.path) #what if threads is removed

qc <- qc_aggregate(qc.dir)
print ("this line is printing, we are making tiff")
#pdf(file.path(qc.dir,"total_seq.pdf"))# ,units="in", width=15, height=15, res=300)
tiff( file.path(qc.dir,"total_seq.tiff") ,units="in", width=15, height=15, res=300)
print ("this line is printing")

g<- ggplot(qc, aes(x=sample, y=tot.seq))+ geom_bar(stat="identity")  + coord_flip()
plot(g)
dev.off()
# pdf(file.path(qc.dir,"qc_heatmap.pdf"), width=15, height=15, res=300)
tiff(file.path(qc.dir,"qc_heatmap.tiff"), width=15, height=15, res=300)

plot(ggplot(qc, aes(x=sample, y =status)) + geom_tile())
dev.off()
print("Multiqc report done")
}

