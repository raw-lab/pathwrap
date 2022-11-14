

#' Title
#'
#' @param fq.dir
#' @param ref.dir
#' @param phenofile
#' @param outdir
#' @param endness
#' @param entity
#' @param corenum
#' @param diff.tool
#' @param compare
#' @param seq_tech
#'
#' @return
#' @export
#'
#' @examples
pathviewwrap <- function(fq.dir="mouse_raw", ref.dir = NA, phenofile= NA, outdir="results", endness="SE",  entity="Mus musculus", corenum = 8, diff.tool="DESEQ2", compare="unpaired", seq_tech="Illumina"){
 dirlist <- sanity_check(fq.dir, ref.dir , phenofile, outdir, endness,  entity , corenum , diff.tool, compare)

 qc.dir <- dirlist[1]
    trim.dir <- dirlist[2]
    sampleFile <- dirlist[3]
    genomeFile<- dirlist[4]
    geneAnnotation <- dirlist[5]

    deseq2.dir <- dirlist[6]
    gage.dir <- dirlist[7]
    grp.idx <- dirlist[8:length(dirlist)]

    run_qc(fq.dir, qc.dir, corenum)

      #call function for quality trimming
    library(parallel)
    #setwd(fq.dir)
    print("calling fastp")
    #samaplenamelist <- read.csv( sampleFile , header =T, sep ="\t")$SampleName
    # for (names in samaplenamelist){
    #   rfastp(read1 = paste0(names, ".fastq"), read2 = "", outputFastq = paste0(names, "_outfastq"), merge = FALSE, adapterFasta= "adapters.fasta", thread = corenum)
    #         }
    cl <- makeCluster(corenum)
    seq_tech = seq_tech
    clusterExport(cl,c("fq.dir","endness","seq_tech", "trim.dir"), envir = environment())#.GlobalEnv)
    ans <- parSapply(cl , read.csv( sampleFile , header =T, sep ="\t")$SampleName  ,run_fastp )
    print("the trim run is complete")
    stopCluster(cl)
    #make txdb from annotation
    txdb <- make_txdbobj(geneAnnotation, corenum, genomeFile, entity)
    print("the cluster are done" )
    aligned_proj <- run_qAlign(corenum, endness, sampleFile, genomeFile,geneAnnotation, ref.dir) #can be better
    geneLevels <-run_qCount(genomeFile, geneAnnotation, aligned_proj, corenum, outdir, txdb)
    exp.fcncnts <-run_difftool(diff.tool = "DESEQ2",outdir, grp.idx, geneLevels, entity, deseq2.dir)
    run_pathway(entity,exp.fcncnts [1] , compare, gage.dir, exp.fcncnts [2], grp.idx)
}
