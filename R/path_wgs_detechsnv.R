#Calling copy number and structural variants


run_alignment <- function(ref.dir, sample.dir, corenum, samplename)
{
  cmd1 <- paste0("bcftools mpileup -O b -o",  results/bcf/sample_name_tosed_raw.bcf, "data/ref_genome/ecoli_rel606.fasta",
                results/bam/sample_name_tosed.aligned.sorted.bam   )
  cmd1 <- stringr::str_replace_all(cmd1, "samplename_to_sed", samplename)


cmd2 <-  paste0("bcftools call --ploidy 1 -m -v -o results/vcf/" ,sample_name_tosed_variants.vcf ," results/bcf/sample_name_tosed_raw.bcf ")
cmd2 <- stringr::str_replace_all(cmd2, "samplename_to_sed", samplename)
cmd3 <- paste0("vcfutils.pl varFilter results/vcf/sample_name_tosed_variants.vcf  > results/vcf/sample_name_tosed_final_variants.vcf")
cmd3 <- stringr::str_replace_all(cmd3, "samplename_to_sed", samplename)
system(cmd1)
system(cmd2)
system(cmd3)
}
