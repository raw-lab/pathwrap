

#' Title
#'
#' @param fq.dir : path of folder where raw fastq files are : one fasta per sample
#' @param ref.dir : path to reference directory which contain reference file(*.fa) and annotation file(*.gtf), can be NULL
#' @param phenofile : path to phenofile ; see note on test data to see the format of phenofile 
#' @param outdir  : give the name of parent result dir , this can be existing or not , rest of directory will be formed by program for organization
#' @param endness : can be “PE” for paired and “SE” for unpaired data
#' @param entity  : Scientific name of species whose RNA is being analyzed
#' @param corenum : number of cores available for analysis #defaut 2
#' @param diff.tool :  what differential tool to to use, “DESEQ2” or “edgeR” available
#' @param compare : what is sample/experimental design you have, paired or unpaired
#' @param seq_tech : Illumina, pacbio or nanopore
#'
#' @return
#' @export
#'
#' @examples pathviewwrap(fq.dir = "~/Documents/Research/UNCC/old/mouse/fresh/mouse_raw", ref.dir = "~/Documents/Research/UNCC/old/Data/Reference/mouse", phenofile = "~/Documents/Research/UNCC/old/mouse/fresh/mouse_raw/pheno.txt", outdir="~/Documents/Research/UNCC/Research_rotation-II/newtmp/pathviewwrap/results_withoutref", endness = "SE", entity = "Mus musculus", corenum = 8, diff.tool  = "DESEQ2", compare = "unpaired", seq_tech = "Illumina")


pathviewwrap <- function(fq.dir="mouse_raw", ref.dir = NA, phenofile= NA, outdir="results", endness="SE",  entity="Mus musculus", 
                         corenum = 8,  compare="unpaired",diff.tool = "DESeq2", seq_tech="Illumina", keep_tmp = FALSE,rerun = FALSE ){
  
    setwd(outdir)
    dirlist <- sanity_check(fq.dir, ref.dir , phenofile, outdir, endness,  entity , corenum , compare, rerun)
    
    coldata <- as.data.frame(dirlist[9:10])
    rownames(coldata) <- str_remove(coldata$Sample, pattern=".fastq.gz")

    dirlist <- unlist(dirlist)
    qc.dir <- dirlist[1]
    print(qc.dir)
    print("this is qc dir")
    trim.dir <- dirlist[2]
    sampleFile <- dirlist[3]
    genomeFile <- dirlist[4]
    geneAnnotation <- dirlist[5]
    deseq2.dir <- dirlist[6]
    edger.dir <- dirlist[7]
    gage.dir <- dirlist[8]
   
    if (!file.exists(file.path(qc.dir,"qc_heatmap.tiff"))){
      print("STEP 1 ; running fastqc")
      print("this is qc.dir")
      print(qc.dir)
      run_qc(fq.dir, qc.dir, corenum)
    }

   
    print("calling fastp")
    cl <- makeCluster(corenum)
    seq_tech = seq_tech
    clusterExport(cl,c("fq.dir","endness","seq_tech", "trim.dir"), envir = environment())#.GlobalEnv)
    ans <- parSapply(cl , read.csv( sampleFile , header =T, sep ="\t")$SampleName  ,run_fastp )
    print("the trim run is complete")
    stopCluster(cl)
    
    #make txdb from annotation
    if(!file.exists(paste0(outdir, gsub(" ", "", entity), "_txdbobj"))){
      print("STEP 2; making txdb obj")
      txdb <- make_txdbobj(geneAnnotation, corenum, genomeFile, entity, outdir)
    } else{
      txdb <- AnnotationDbi::loadDb(paste0(outdir,"/", gsub(" ", "", entity), "_txdbobj"))
    }
    
    if(!file.exists(file.path(aligned_bam , "alltrimmedalignedobj.RDS"))){
      setwd(outdir)#make sure you delete this file before rerunning can be better
      print("STEP 3 : aligning the sequence")
      aligned_proj <- run_qAlign(corenum, endness, sampleFile, genomeFile,geneAnnotation, ref.dir) #can be better?? 
    }
    else{
      aligned_proj <- readRDS(file.path(aligned_bam , "alltrimmedalignedobj.RDS"))
    }
    
     if(!file.exists(file.path(outdir, "combinedcount.trimmed.RDS")  ))  { 
    
      # why is this not recognizing alignemtn already present
      print("STEP 4: counting aligned sequences")
      returnlistcntsindx <-run_qCount(genomeFile, geneAnnotation, aligned_proj, corenum, outdir, txdb, entity, coldata)
      cnts <- returnlistcntsindx$cnts
      grp.idx <- returnlistcntsindx$grp.idx
    }
    else{
      cnts <- as.data.frame(readRDS(file.path(outdir, "combinedcount.trimmed.RDS") )) #check
      if(  all(rownames(coldata) == colnames(cnts)) ){#if this then proceed
        ref <- which(coldata[, 2] ==  levels(coldata[, 2])[1])
        samp <- which(coldata[, 2] ==  levels(coldata[, 2])[2])
        grp.idx <-NULL
        grp.idx[ref] <- "reference"
        grp.idx[samp] <- "sample"
      
      }
    }
   
    if (keep_tmp == FALSE){
      print("deleting aligned bam files, bam file index and log files")
      unlink(file.path(outdir, "aligned_bam", "*bam*"))
    }

    if(!file.exists(paste0(deseq2.dir, "/Volcano_deseq2.tiff"))){
      print("STEP 5a ; running differential analysis using DESeq2")
      exp.fcncnts.deseq2 <- run_deseq2(cnts,grp.idx, deseq2.dir)
      print(head(exp.fcncnts.deseq2))
    }
    else{
      deseq2.res.df  <- read.table(file.path(deseq2.dir, "DESEQ2_logfoldchange.txt"), header = T, sep = "\t", row.names = 1) #works with gage
      exp.fcncnts.deseq2 <- deseq2.res.df  $log2FoldChange
      names( exp.fcncnts.deseq2) <-  rownames(deseq2.res.df )
     
      }
    if(!file.exists(paste0(edger.dir, "Volcano_edgeR.tiff"))){
      print("STEP 5b ; running differential analysis using edgeR")
      exp.fcncnts.edger <- run_edgeR(cnts,grp.idx, edger.dir)
    } else{
      edger.res.df  <- read.table(file.path(edger.dir, "edgeR_logfoldchange.txt"), header = T, sep = "\t", row.names = 1) #works with gage
      exp.fcncnts.deseq2 <-edger.res.df  $log2FC
      names( exp.fcncnts.deseq2) <-  rownames(edger.res.df )
    }
    
    setwd(gage.dir)
    #chosing to use deseq2 result or edger result for gage
    if(diff.tool == "DESeq2"){
      exp.fc <- exp.fcncnts.deseq2
    } else{
      exp.fc <- exp.fcncnts.edger
    }
    if(!file.exists("*.txt")){
      print("STEP 6 : running pathway analysis using GAGE")
      print(paste0(compare, "this is from pathviewwrap"))
      run_pathway(entity,exp.fc, compare, gage.dir, cnts, grp.idx) 
      }
}
