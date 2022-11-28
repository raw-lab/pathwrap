prebuild <- function(ref_genome, ref.dir){
  setwd(ref_dir)
  system(paste0("bwa index -p ref_genome", ref_genome))

}
