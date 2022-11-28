path_wgs_markdups <- function(){
  setwd(align.dir)
  cmd <- paste0("java -Xmx8G -jar $PICARD/picard-2.8.0.jar MarkDuplicates INPUT=samplename_to_sed_sorted.sam OUTPUT=samplename_to_sed_sorted_marked.bam METRICS_FILE=metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT")
  cmd <- stringr::str_replace_all(cmd, "samplename_to_sed", samplename)

  system(cmd)
}
