path_wgs_picard <- function(){
  setwd(align.dir)
  cmd <- paste0("java -Xmx8G -jar $PICARD/picard-2.8.0.jar SortSam INPUT= samplename_to_sed.sam OUTPUT=samplename_to_sed_sorted.sam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT")
  cmd <- stringr::str_replace_all(cmd, "samplename_to_sed", samplename)

  system(cmd)
}

