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
  library(stringr)

  if (seq_tech == "PacBio" | seq_tech == "Nanopore" ){ #use custom adapters

    if (endness=="PE"){

      cmd <- paste0("fastp -i " ,file.path(fq.dir ,paste0( "samplename_to_sed" , filenamepattern)) , 
                    " -I ", file.path(fq.dir , paste0("samplename_to_sed", str_replace(filenamepattern, "1", "2"))),
                    " -o ", file.path(trim.dir , paste0("samplename_to_sed_trimmed", filenamepattern )),
                    " -O ", file.path(trim.dir , paste0("samplename_to_sed_trimmed", str_replace(filenamepattern, "1", "2"))),
                    "--adapter_fasta", "data/adapters.fna",
                    " -h " , file.path(trim.dir, "samplename_to_sed.html") ,
                    " -j " , file.path(trim.dir,  "samplename_to_sed.json"))

    } else {
      cmd <- paste0("fastp -i ",
                    file.path(fq.dir , paste0( "samplename_to_sed" , filenamepattern)) ,
                    " -o ",  file.path(trim.dir ,  paste0("samplename_to_sed_trimmed", filenamepattern )),
                    " --adapter_fasta ", "data/adapters.fna", " -h " ,
                    file.path(trim.dir,  "samplename_to_sed.html"),
                    " -j " , file.path( "samplename_to_sed.json"))
    }
  
    
  } else {
    if (endness=="PE"){

      cmd <- paste0("fastp -i " ,
                    file.path(fq.dir , paste0( "samplename_to_sed" , filenamepattern)) , " -I ",
                    file.path(fq.dir , paste0("samplename_to_sed", str_replace(filenamepattern, "1", "2"))), " -o ",
                    file.path(trim.dir , paste0("samplename_to_sed_trimmed", filenamepattern )), " -O ",
                    file.path(trim.dir , paste0("samplename_to_sed_trimmed", str_replace(filenamepattern, "1", "2"))),   " -h " ,
                    file.path(trim.dir, "samplename_to_sed.html") , " -j " ,
                    file.path(trim.dir ,  "samplename_to_sed.json"))

    } else {
      print("fastp running for single end illumina running")
      cmd <- paste0("fastp -i ",
                    file.path(fq.dir , paste0( "samplename_to_sed" , filenamepattern)) , " -o ",
                    file.path(trim.dir , paste0("samplename_to_sed_trimmed", filenamepattern )) , " -h " ,
                    file.path(trim.dir,  "samplename_to_sed.html"),  " -j " ,
                    file.path(trim.dir , "samplename_to_sed.json"))

    }
  }
  cmd <- stringr::str_replace_all(cmd, "samplename_to_sed", samplename)

  if(!file.exists(file.path(trim.dir , str_replace_all("samplename_to_sed.json", "samplename_to_sed", samplename)))){
    print("STEP 1b : running fastp")
    print(cmd)
        system(cmd)
  }

}
