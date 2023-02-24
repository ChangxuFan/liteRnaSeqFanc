mini.align <- function(fastq.vec, R1.pattern="1.fastq.gz", R2.pattern = "2.fastq.gz", k = NULL, a = F,
                       bowtie.genome, map.bed = NULL, out.bam, bowtie.log.file=NULL, thread = 16, run = T,
                       bowtie2 = "/opt/apps/bowtie2/2.3.4.1/bowtie2", no.unal = F) {
  system(paste0("mkdir -p ", dirname(out.bam)))
  
  R1 <- fastq.vec[grepl(R1.pattern, fastq.vec)][1]
  R2 <- fastq.vec[grepl(R2.pattern, fastq.vec)][1]
  # cmd <- paste0("/bar/cfan/scripts/dna-seq/bowtie2_4_22_20.sh ", " -p ", thread, 
  #               " -x ", bowtie.genome, " -i ", R1, " -I ", R2, " -o ", out.bam, " -k ", 50, 
  #               " -s ", map.bed[1,2])
  cmd <- paste0(bowtie2, " -X2000 --reorder --very-sensitive --xeq --seed 42 ",
                " -p ", thread, " -x ", bowtie.genome, " -1 ", R1, " -2 ", R2)
  
  # if (!is.null(bowtie.log.file))
  #   cmd <- paste0(cmd, " --met-file ", bowtie.log.file)
  
  if (!is.null(k)) 
    cmd <- paste0(cmd, " -k ", k)
  if (a == T)
    cmd <- paste0(cmd, " -a ")
  if (no.unal == T)
    cmd <- paste0(cmd, " --no-unal")
  
  if (!is.null(map.bed)) {
      shift <- map.bed[1,2] 
      cmd <- paste0(cmd , " | ", "awk -F \"\\t\" 'BEGIN {OFS = \"\\t\"} {$4 = $4+",shift,"; print $0}' ")
  }

  cmd <- paste0(cmd, " | /bar/cfan/anaconda2/envs/jupyter/bin/samtools sort -O bam -m 4G -@ ",thread," -  > ", out.bam)
  
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run)
  if (!file.exists(out.bam))
    stop(paste0(out.bam, " was not successfully generated"))
  
  cmd <- paste0("samtools index ", out.bam)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run, stdout.file = bowtie.log.file)
  if (!file.exists(out.bam %>% paste0(".bai")))
    stop(paste0(out.bam, ".bai was not successfully generated"))
  
  return(out.bam)
}