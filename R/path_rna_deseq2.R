#' Title
#'
#' @param cnts
#' @param grp.idx
#' @param deseq2.dir
#'
#' @import EnhancedVolcano EnhancedVolcano
#' @import DESeq2
#'
#' @return
#' @export
#'
#' @examples
run_deseq2 <- function(cnts,grp.idx, deseq2.dir){
  library(DESeq2)
  coldat=DataFrame(grp=factor(grp.idx))
  dds <- DESeqDataSetFromMatrix(cnts, colData=coldat, design =~ grp)
  dds <- DESeq(dds)
  deseq2.res <- results(dds)
  #direction of fc, depends on levels(coldat$grp), the first level
  #taken as reference (or control) and the second one as experiment.
  deseq2.fc=deseq2.res$log2FoldChange
  names(deseq2.fc)=rownames(deseq2.res)
  exp.fc<- deseq2.fc
  table(is.na(deseq2.res$padj))
  write.table(deseq2.res , file.path(deseq2.dir, "DESEQ2_logfoldchange.txt"),sep = "\t" ,col.names =NA, row.names =TRUE, quote =FALSE)
  tiff(file.path(deseq2.dir, "Volcano_deseq2.tiff"), units="in", width=15, height=15, res=300)
  plot(EnhancedVolcano::EnhancedVolcano(deseq2.res, x ='log2FoldChange', y ='pvalue', lab =rownames(deseq2.res)))
  dev.off()
  return( exp.fc)
}
