#run fastp
#' to run trimming
#'
#' @param samplename
#'
#' @return
#' @export
#'
#' @examples 
run_fastp <-function(samplename){

  if (seq_tech == "PacBio" | seq_tech == "Nanopore" ){ #use custom adapters

    if (endness=="PE"){

      cmd <- paste0("fastp -i " ,file.path(fq.dir , "samplename_to_sed_1.fastq.gz") ,
                    " -I ", file.path(fq.dir , "samplename_to_sed_2.fastq.gz"),
                    " -o ", file.path(trim.dir , "samplename_to_sed_paired_1.fastq.gz"),
                    " -O ", file.path(trim.dir , "samplename_to_sed_paired_2.fastq.gz"),
                    "--adapter_fasta", "data/adapters.fna",
                    " -h " , file.path(trim.dir, "samplename_to_sed.html") ,
                    " -j " , file.path(trim.dir,  "samplename_to_sed.json"))

    } else {
      cmd <- paste0("fastp -i ",
                    file.path(fq.dir , "samplename_to_sed.fastq.gz"),
                    " -o ",  file.path(trim.dir , "samplename_to_sed_paired.fastq.gz"),
                    " --adapter_fasta ", "data/adapters.fna", " -h " ,
                    file.path(trim.dir,  "samplename_to_sed.html"),
                    " -j " , file.path( "samplename_to_sed.json"))
    }




  } else {
    if (endness=="PE"){

      cmd <- paste0("fastp -i " ,
                    file.path(fq.dir , "samplename_to_sed_1.fastq.gz") , " -I ",
                    file.path(fq.dir , "samplename_to_sed_2.fastq.gz"), " -o ",
                    file.path(trim.dir , "samplename_to_sed_paired_1.fastq.gz"), " -O ",
                    file.path(trim.dir , "samplename_to_sed_paired_2.fastq.gz"),   " -h " ,
                    file.path(trim.dir, "samplename_to_sed.html") , " -j " ,
                    file.path(trim.dir ,  "samplename_to_sed.json"))

    } else {
      cmd <- paste0("fastp -i ",
                    file.path(fq.dir , "samplename_to_sed.fastq.gz"), " -o ",
                    file.path(trim.dir , "samplename_to_sed_paired.fastq.gz"), " -h " ,
                    file.path(trim.dir,  "samplename_to_sed.html"),  " -j " ,
                    file.path(trim.dir , "samplename_to_sed.json"))

    }
  }
  cmd <- stringr::str_replace_all(cmd, "samplename_to_sed", samplename)
  print(cmd)
  if(!file.exists(file.path(trim.dir , "samplename_to_sed.json"))){
    print(paste0(file.path(trim.dir , "samplename_to_sed.json"), "does not exit")
    system(cmd)
  }

}
