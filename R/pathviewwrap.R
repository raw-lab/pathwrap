

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
                         corenum = 8, diff.tool="DESEQ2", compare="unpaired", seq_tech="Illumina"){
    dirlist <- unlist(sanity_check(fq.dir, ref.dir , phenofile, outdir, endness,  entity , corenum , diff.tool, compare))

    qc.dir <- dirlist[1]
    trim.dir <- dirlist[2]
    sampleFile <- dirlist[3]
    genomeFile<- dirlist[4]
    geneAnnotation <- dirlist[5]

    deseq2.dir <- dirlist[6]
    gage.dir <- dirlist[7]
    coldata <- dirlist[8:9]
    #grp.idx <- dirlist[8:length(dirlist)]

    run_qc(fq.dir, qc.dir, corenum)

     #call function for quality trimming
    library(parallel)
    #setwd(fq.dir)
    print("calling fastp")

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
    exp.fcncnts <- run_difftool(diff.tool = "DESEQ2",outdir,coldata, geneLevels, entity, deseq2.dir)
    run_pathway(entity,exp.fcncnts [1] , compare, gage.dir, exp.fcncnts [2], exp.fcncnts [2]) # see if you can use grp.idx
}
