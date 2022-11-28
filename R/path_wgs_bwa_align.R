run_alignment <- function(ref.dir, sample.dir, corenum, samplename)
  {
              cmd <- paste0("bwa mem -M -t", corenum, " " , ref.dir, "/ref_genome ",
                sample.dir, "/sample_name_to_sed_1.fq " ,
              sample.dir, "/sample_name_to_sed_2.fq " ,
                "| samtools view -Sb sample_name_to_sed.bam"  )
            cmd <- stringr::str_replace_all(cmd, "samplename_to_sed", samplename)
              system(cmd)

}
