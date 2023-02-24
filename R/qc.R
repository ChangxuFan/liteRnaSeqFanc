#
trim.fastq <- function (fastq, trim_dir=F) {
  if (trim_dir == T) fastq <- basename(fastq)
  trim <- sub("(.*)(.fastq.gz)", "\\1", fastq)
  return(trim)
}

trim.bam <- function (bam, trim_dir=F) {
  if (trim_dir == T) bam <- basename(bam)
  trim <- sub("(.*)(.bam)", "\\1", bam)
  return(trim)
}


fastqc <- function(fastq, fastqc.dir=NULL, thread) {
  if (is.null(fastqc.dir)) fastqc.dir="./fastqc"
  system(paste0("mkdir -p ", fastqc.dir))
  system(paste0("/opt/apps/FastQC/0.11.9/fastqc -q --extract -o ", fastqc.dir, " -t ", thread,
                " -a ~/scripts/rna-seq/adapters_and_indices.txt ",
                fastq))

  return(paste0(fastqc.dir, "/", trim.fastq(fastq = fastq, trim_dir = T), "_fastqc/") %>% normalizePath(mustWork = F))
}

multi.fastqc <- function(fastq.vector, fastqc.dir=NULL, thread ) {
  if (is.null(fastqc.dir)) fastqc.dir="./fastqc"
  system(paste0("mkdir -p ", fastqc.dir))
  file.list <- paste0(fastq.vector, collapse = " ")
  system(paste0("/opt/apps/FastQC/0.11.9/fastqc -q --extract -o ", fastqc.dir, " -t ", thread,
                " -a ~/scripts/rna-seq/adapters_and_indices.txt ",
                file.list))
  fastqc.out <- sapply(fastq.vector, function (x) {
    paste0(fastqc.dir, "/", trim.fastq(fastq = x, trim_dir = T), "_fastqc/") %>% normalizePath(mustWork = F) %>%
      return()
  })
  return(fastqc.out)
}

fastqc.parse <- function(fastqc.dir.sample, over.only = F) {
  system(paste0("~/scripts/rna-seq/parse_fastqc_overrepresented.sh ", fastqc.dir.sample, "/fastqc_data.txt", " ",
                fastqc.dir.sample,"/over.txt"))
  over <- read.table(paste0(fastqc.dir.sample,"/over.txt"), as.is = TRUE, sep="\t")

  over.df <- data.frame(head=paste0(">add_", runif(nrow(over), 1, 1000) %>% round()), seq = over$V1)
  over.dir <- paste0(fastqc.dir.sample, "/", "over.fasta")
  write.table(over.df, file = over.dir, sep = "\n", quote = F, row.names = F, col.names = F)
  if (over.only == T) {
    return(over.dir)
  }
  adapter.dir <- paste0(fastqc.dir.sample, "/fastp_adapters.fasta") %>% normalizePath(mustWork = F)
  system(paste0("cat ", "~/scripts/rna-seq/adapter.fasta ", over.dir, " >", adapter.dir))
  return(adapter.dir)
}




fastp.se <- function(adapter.dir, fastq, thread, fastp.out.dir=NULL) {
  if (is.null(fastp.out.dir)) fastp.out.dir <- dirname(fastq)

  fastp.out <- paste0(fastp.out.dir, "/", trim.fastq(fastq = fastq, trim_dir = T), "_fastp_fastq.gz") %>% normalizePath(mustWork = F)

  system(paste0("fastp -i ", fastq, " -o ", fastp.out, " -x A -x T -w ", thread, " --adapter_fasta ", adapter.dir))
  return(fastp.out)
}

fastp.pe <- function (adapter.dir.pair, fastq.pair, thread, fastp.out.dir = NULL) {
  if (is.null(fastp.out.dir)) fastp.out.dir <- dirname(fastq.pair[1])
  fastp.out <- sapply(fastq.pair, function (fastq) {
    paste0(fastp.out.dir, "/", trim.fastq(fastq = fastq, trim_dir = T), "_fastp_fastq.gz") %>% normalizePath(mustWork = F) %>%
      return()
  })
  # cat 2 over.txt into one
  adapter.dir <- paste0(fastp.out.dir, "/adapters/", trim.fastq(fastq.pair[1], trim_dir = T), "_adapters.fasta")
  system(paste0("mkdir -p ",paste0(fastp.out.dir, "/adapters/") ))
  system(paste0("cat ", "~/scripts/rna-seq/adapter.fasta ", adapter.dir.pair[1], " ", adapter.dir.pair[2], " > ", adapter.dir))
  system(paste0("fastp -i ", fastq.pair[1], " -o ", fastp.out[1]," -I ", fastq.pair[2], " -O ", fastp.out[2],
                " -x A -x T -w ", thread, " --adapter_fasta ", adapter.dir))
  return(fastp.out)
}



qc.se <- function(fastq, fastqc.dir, fastp.out.dir, thread) {
  fastqc.dir.sample <- fastqc(fastq = fastq, fastqc.dir = fastqc.dir, thread=thread)
  adapter.dir <- fastqc.parse(fastqc.dir.sample = fastqc.dir.sample)
  fastp.out <- fastp.se(adapter.dir = adapter.dir, fastq = fastq, thread = thread, fastp.out.dir = fastp.out.dir)
  return(fastp.out %>% normalizePath(mustWork = F))
}



# qc.se(fastq = fastqs[1],
#       fastqc.dir = "~/rna-seq/pipeline/snakemake/NK_rep1_fastqc",
#       fastp.out.dir = "~/rna-seq/pipeline/snakemake/NK_rep1_fastqc",
#       thread = 8)

