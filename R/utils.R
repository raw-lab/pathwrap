#make sure force = T is good
#' Title
#'
#' @param fq.dir
#' @param ref.dir
#' @param phenofile
#' @param outdir
#' @param endness
#' @param entity
#' @param corenum
#' @param compare
#'
#' @return
#' @export
#'
#' @examples
sanity_check <- function(fq.dir, ref.dir , phenofile, outdir, endness,  entity , corenum , compare, rerun){

  library(stringr)

  #################################################################
  ##
  #Check if files/folders  exists and create if not
  ##
  #################################################################
  if (file.exists(outdir) & rerun == T){
    unlink(outdir, recursive = T)
  }
  
  if (!file.exists(outdir)){
    # default output file
    dir.create(outdir)
  }
  result.dir <- outdir
  #print("The results will be organized in ",result.dir)

  # make sure the second column is class and first column is sample name
  # make sure file is tab seperated
  if (!file.exists(phenofile)){ ###TO DO make sure reference is first aplhanumerically#
    print("Please provide phenofile with Class information")
  }
  coldata <- read.table(phenofile, sep = "\t", header = T)
  if(colnames(coldata)[2]!="Class"){
    print("Please make sure class information is in cloumn 2 with colname 'Class' . ")
  }
  coldata$Class <- as.factor(coldata$Class)
  # ref <- which(coldata[, 2] ==  levels(coldata[, 2])[1])
  # samp <- which(coldata[, 2] ==  levels(coldata[, 2])[2])
  # grp.idx <-NULL
  # grp.idx[ref] <- "reference"
  # grp.idx[samp] <- "sample"
  ##TO DO write something to automatically determine paired information, rev/fr etc
  #print("this is first grp.idx")
  #print(grp.idx)
  #check and create dir for organizing results
  checkcretdir <- function(parentname, dirname){
    if(!file.exists(file.path(parentname, dirname))) {
      dir.create(file.path(parentname, dirname))
    }
    assign(dirname,value = file.path(parentname, dirname), envir = .GlobalEnv)
  }

  folder_to_create<- list("fastqc_results", "fastp_results","gage_results", "differential_analysis","aligned_bam","pathway_analysis" )
  trim_dir <-list ( "fastp_log", "unpaired")
  diff_dir <-list ("DESeq2","edgeR")
  pathway_types <- list("KEGG", "GO")
  kegg_types <- list("signalling", "metabolism", "disease", "sig_n_met")
  go_types <- list("biological_process", "molecular_function", "cellular_component")
  lapply(folder_to_create, checkcretdir, parentname= result.dir  )
  lapply(trim_dir, checkcretdir, parentname= file.path(result.dir ,"fastp_results")  )
  lapply(diff_dir, checkcretdir, parentname= file.path(result.dir , "differential_analysis")  )
  lapply(pathway_types, checkcretdir, parentname= file.path(result.dir , "gage_results")  )
  lapply(kegg_types, checkcretdir, parentname= file.path(result.dir ,"gage_results","KEGG" )  )
  lapply(go_types, checkcretdir, parentname= file.path(result.dir ,"gage_results","GO" )  )

  #just to make sure rest of codes are same
  qc.dir <- fastqc_results
  diff.dir <- differential_analysis
  trim.dir <- fastp_results
  gage.dir <- gage_results
  trim.log <- fastp_log
  pathway.dir <- pathway_analysis
  edger.dir <- edgeR
  deseq2.dir <- DESeq2
  kegg.dir <- KEGG
  go.dir <- GO


  ### To run qAlign we need samplefile
  ##############################################################################

  if( endness== "SE"){
    #pinfo_string <- ".fastq"
    pinfo_string <- ".fastq.gz"
  }else{
   # pinfo_string <- "_1.fastq"
    pinfo_string <- "_1.fastq.gz"
  }
  library(stringr)
  FileName <- grep(pinfo_string,list.files(fq.dir, full.names=T) ,value =T)
  FileName <- str_replace_all(file.path(trim.dir,  basename(FileName)),pinfo_string, paste0("_paired", pinfo_string))
  sampleFile <- file.path(result.dir, "sampleFile.txt")
  SampleName <-  str_remove_all(basename(FileName), paste0("_paired",pinfo_string))
  if(endness == "SE") {
    write.table(file =sampleFile,sep = "\t", as.data.frame( cbind(FileName, SampleName)) ,quote =F ,  col.names=T, row.names=F)
  } else{
    FileName1 <- FileName
    FileName2 <- str_replace_all(FileName1, "_1.fastq.gz", "_2.fastq.gz")
    
    #FileName2 <- str_replace_all(FileName1, "_1.fastq", "_2.fastq")
    sampleFile <- file.path(result.dir, "sampleFile.txt")
    write.table(file =sampleFile,sep = "\t", as.data.frame( cbind(FileName1, FileName2, SampleName)) ,quote =F ,  col.names=T, row.names=F)
  }



  #References
  #if only species name is given and both geneAnnotation and genome is NULL
  if( is.na(ref.dir)){
    #ref_info <- read.table("data/species_genome_annotation_pkg", sep = "\t", header = T, na.strings=c(""," ","NA")) #this file is supplied with script
    #ref_info <- as.data.frame(readRDS("data/anntpkglist.RDS"))#, sep = "\t", header = T, na.strings=c(""," ","NA")) #this file is supplied with script
    ref_info <- anntpkglist

    species_no <- which(ref_info$species==entity)
    annotate_pkg <- ref_info$annotation[species_no]
    genome_pkg <- ref_info$genome[species_no]

    ###make sure both annot and genome package is installed for the species
    # (set of genome and annotation pkg come from developers list)

    pkg.on = require(annotate_pkg, character.only = TRUE, lib.loc = .libPaths()[1])
    if (!pkg.on) {
      if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
      BiocManager::install(annotate_pkg,force = T, suppressUpdates =TRUE, lib.loc = .libPaths()[1] )
      pkg.on = require(annotate_pkg, character.only = TRUE, lib.loc = .libPaths()[1])
      if (!pkg.on)
        stop(paste("Fail to install/load gene annotation package ",annotate_pkg, "!", sep = ""))
    }
    geneAnnotation <-  file.path(.libPaths()[1],annotate_pkg, "extdata", paste0(annotate_pkg, ".sqlite" ) )
    genomeFile <- genome_pkg
  } else {
    genomeFile <- list.files(ref.dir, ".fa$", full.names= T)
    geneAnnotation <- list.files(ref.dir, ".gtf$", full.names = T) #could be changed to include one of gtf, gff etc, check with quasR package

  }
  return (c(qc.dir,trim.dir,sampleFile, genomeFile, geneAnnotation, deseq2.dir, edger.dir, gage.dir, coldata))
}

