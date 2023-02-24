# the logic: pair end qc must start with a sample info sheet. one additional field is required: "pe_flag".
##it can be set as either se or pe. if se, the sample will be processed just as single end samples,
##if pe, the sample will be processed as pair end samples.
##this mode also supports split files. which means that GSMs with split SRRs will also be accommodated.

qc.pe <- function(fastq.pair, fastqc.dir=NULL, fastp.out.dir=NULL, thread) {
  fastqc.dir.sample <- multi.fastqc(fastq = fastq.pair, fastqc.dir = fastqc.dir, thread=thread)
  adapter.dir.pair <- sapply(fastqc.dir.sample, function (x) fastqc.parse(fastqc.dir.sample=x, over.only = T))
  fastp.out <- fastp.pe(adapter.dir.pair = adapter.dir.pair, fastq = fastq.pair, thread = thread, fastp.out.dir = fastp.out.dir)
  return(fastp.out %>% normalizePath(mustWork = F))
}

