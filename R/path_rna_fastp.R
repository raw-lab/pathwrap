#run fastp
#' to run trimming
#'
#' @param sampleName
#'
#' @return
#' @export
#'
#' @examples
run_fastp <-function(sampleName){
  intformatch <- grep(sampleName, FileName[,1], value = F, fixed = T) #integer
  trimmedoutfile <-  FileName[intformatch,] #trimmedoutfile$FileName1 , trimmedoutfile$Filename2
  infile <- filenames[intformatch,]  #infile$FileName1 , infile$FileName2 #infile for SE
  if (endness=="PE"){
    infileoutfile <- paste0("-i ", infile$FileName1 ,   " -I ", infile$FileName2,   " -o ", trimmedoutfile["FileName1"],  " -O ", trimmedoutfile["FileName1"])
  } else {
    infileoutfile <- paste0("-i ", infile, " -o ", trimmedoutfile)
  }

 logfile <- paste0( " -h " , file.path(trim.dir, paste0(sampleName, ".html")) ,  " -j " , file.path(trim.dir,  paste0(sampleName, ".json")))

 #actual command
  if (seq_tech == "PacBio" | seq_tech == "Nanopore" ){ #use custom adapters
    cmd <- paste0("fastp " , infileoutfile,  "--adapter_fasta", "data/adapters.fna", logfile)
  } else {
    cmd <- paste0("fastp " , infileoutfile, logfile)
  }
 #file check before running command
  if (length(list.files(trim.dir, pattern ="json")) <=0 | length(list.files(trim.dir, pattern ="trimmed")) != length(list.files(trim.dir, pattern ="json"))){
     print("STEP 1b : running fastp")
    print(cmd)
    system(cmd)
  }
 }
