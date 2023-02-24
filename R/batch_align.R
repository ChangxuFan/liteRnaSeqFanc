bwa.batch <- function(fastq.vec, pe.regex,
                      genome.fa.vec, genome.names = NULL, genome.master.dir = NULL,
                      mapped.only = T, bwa,
                      master.dir, threads.genome = 3, 
                      threads.fastq = 3, threads.bwa = 2, add.params = "") {
  if (!is.null(genome.master.dir) && !is.null(genome.names))
    genome.fa.vec <- paste0(genome.master.dir, "/", genome.names, "/", genome.names,".fa")
  
  if (is.null(genome.names))
    genome.names <- names(genome.fa.vec)
  if (is.null(genome.names))
    stop("genome names must be specified")
  if (is.null(names(genome.fa.vec)))
    names(genome.fa.vec) <- genome.names
  
  fastq.groups <- liteRnaSeqFanc::fastq.group.gen(fastqs = fastq.vec, pe.regex = pe.regex, return.name = T)
  names(fastq.vec) <- fastq.groups

  bams <- mclapply(genome.names, function(genome.name) {
    fa <- genome.fa.vec[genome.name]
    bams.fastq <- fastq.vec %>% split(f = fastq.groups) %>% mclapply(function(fastq) {
      if (length(fastq) != 2)
        stop("expecting pair end!!")
      out.bam <- paste0(master.dir, "/", genome.name, "/", names(fastq)[1], ".bam")
      bam.log <- paste0(master.dir, "/", genome.name, "/", names(fastq)[1], ".log")
      liteRnaSeqFanc::bwa.mem(fastq.vec = fastq, ref.fa = fa, out.bam = out.bam,
                              mapped.only = mapped.only,
                              bwa = bwa, threads = threads.bwa, log.file = bam.log,
                              add.params = add.params)
      if (!file.exists(out.bam))
        stop(paste0(out.bam, " failed to generate"))
      return(out.bam)
    }, mc.cores = threads.fastq, mc.cleanup = T)
    
    return(bams.fastq %>% unlist())
  }, mc.cores = threads.genome, mc.cleanup = T)
  names(bams) <- genome.names
  return(bams)
}


