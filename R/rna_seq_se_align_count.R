insert.name <- function(name, insert, ext, trim.dir=F) {
  if (trim.dir == T) name <- basename(name)
  inserted <- sub(paste0("(.*)(", ext, ")"), paste0("\\1_", insert, "\\2"), name)
  return(inserted)
}

rna.se.align <- function(fastq, genome.index, align.out.dir, thread) {
  align.out.dir <- normalizePath(align.out.dir, mustWork = F)
  root.name <- trim.fastq(fastq=fastq, trim_dir=T)
  bam <- paste0(align.out.dir, "/", root.name, "_star_genome.bam")
  trash <- system(paste0("~/scripts/rna-seq/encode_se_star.sh ",genome.index,
                         " ", fastq, " ",root.name, " ", thread, " ",  align.out.dir, " ", root.name ))
  return(bam)
}

# t.align <- rna.se.align("mini.fastq.gz", "~/genomes/mm10/STAR_gencode_vM24", "mini_star", 2)


rna.se.rmdup <- function(bam, pseudo=F) {
  # if pseudo ==T, the function will not perform dedup and will be just a file name formatter.
  bam.nodup <- insert.name(name = bam, insert = "rmdup", ext=".bam", trim.dir = F)
  trash <- system(paste0("samtools rmdup -s ", bam, " ", bam.nodup))
  if (file.exists(bam.nodup)) return(bam.nodup)
  else stop("something wrong at samtools rmdup; output file was not generated")
}


# t.rmdup <- rna.se.rmdup("star/mini.1_star/mini.1_fastp_star_genome.bam")

# rna.se.featureCounts <- function(bam=NULL, annot, count.out.name=NULL, thread, bam.info) {
#   # if bam is null: file names of the bam files will be read from bam.info
#   #bam.info: a sample name ~ bamfile name correspondance file. tsv format.
#   #must have a colum named "sample" and another column named "bamfile"
#   #bam.info can also be a dataframe mapping bam file names to sample names.
#   if (is.null(bam))
#     bam <- read.table(bam.info, as.is = T, header = T, sep="\t") %>% pull(bamfile)
#   bam.expand <- paste0(bam, collapse = " ")
#   trash <- system(paste0("~/scripts/rna-seq/featureCounts_se.sh ", "-a ", annot,
#                          " -o ", count.out.name, " -t ", thread, " ", bam.expand))
#   count.formatted <- count.table.change.name(count.table = count.out.name, bam.info = bam.info)
#   if(file.exists(count.formatted)) return(count.formatted)
#   else return("error: count table was not successfully generated")
# 
# }

# t.count <- t.align %>% rna.se.featureCounts(annot = "~/genomes/mm10/gencode/gencode.vM24.annotation.gtf",
#                                             count.out.name = "mini_count", thread = 2,
#                                             bam.info = data.frame(sample="mini_0", bamfile=t.align))

count.table.change.name <- function(count.table, bam.info) {
  if (is.character(bam.info))
    bam.info <- read.table(bam.info, as.is = T, header = T, sep="\t")
  if (!is.data.frame(bam.info)) stop("bam.info is neither a file path nor a dataframe")
  # print(head(count.table))
  # print(bam.info)
  bam.info$bamfile <- sapply(bam.info$bamfile, function (x) normalizePath(x, mustWork = F))
  ori.count.table <- count.table
  count.table <- read.table(count.table, as.is=T, header = F)
  count.table[, paste0("V", 2:6)] <- NULL
  colnames(count.table) <- paste0("V", 1:ncol(count.table))
  bams <- as.character(count.table[1,])
  bams[c(2:length(bams))] <- normalizePath(bams[c(2:length(bams))], mustWork = T)
  tmpt.df <- data.frame(bamfile = bams)
  converted <- left_join(tmpt.df, bam.info) %>% pull(sample)
  converted[1] <- "gene"
  count.table[1,] <- converted
  write.table(count.table, paste0(ori.count.table, ".renamed"), quote = F, col.names = F, row.names = F, sep = "\t")
  return(paste0(ori.count.table, ".renamed") %>% normalizePath(mustWork = F))
}




rna.se.qc.align.count <- function(sample.info=NULL, qc=T, fastq=NULL, rmdup = F,  genome.index=NULL, annot=NULL, count.out.name,
                                  fastqc.dir="./qc", fastp.out.dir="./", thread=8,
                                   align.out.dir="./star") {
# this function accepts infomation given separately via fastq, genome.index and annot.
##alternatively you can feed a sample info tsv file or sample info dataframe. the corresponding column names should be:
##sample, fastqfile, genomeIndex, annot.
  if (is.character(sample.info))
    sample.info <- read.table(sample.info, as.is = T, header = T, sep="\t")
  fastq <- sample.info$fastqfile
  samples <- sample.info$sample
  genome.index <- sample.info$genomeIndex %>% unique()
  if (!length(genome.index) ==1 ) stop("all the files must be aligned to the same genome and thus must have the same index file")
  annot <- sample.info$annot %>% unique()
  if (!length(annot) ==1 ) stop("all the files must use the same annotation file")

  bams <- sapply(fastq, function (x) {
    if (qc==T)
      fastp <- qc.se(fastq = x, fastqc.dir = fastqc.dir, fastp.out.dir = fastp.out.dir, thread = thread)
    else
      fastp <- x
    bam <- rna.se.align(fastq = fastp, genome.index = genome.index,
                 align.out.dir = paste0(align.out.dir, "/", trim.fastq(fastp, trim_dir = T), "_star/"), thread = thread)
    return(bam)
  })
  bam.info <- data.frame(sample=samples, bamfile=bams)
  write.table(bam.info, paste0(count.out.name, "_bam.info"), col.names = T, row.names = F, quote = F, sep="\t")
  count.table <- rna.se.featureCounts(bam=bams, annot = annot, count.out.name = count.out.name, thread = thread,
                                      bam.info = bam.info)
  bams.dedup <- sapply(bams, function (x) return(rna.se.rmdup(bam=x, pseudo = F)))
  bam.info.dedup <- data.frame(sample=samples, bamfile=bams.dedup)
  write.table(bam.info.dedup, paste0(count.out.name, "_bam_dedup.info"), col.names = T, row.names = F, quote = F, sep="\t")
  count.table.dedup <- rna.se.featureCounts(bam=bams.dedup, annot = annot, count.out.name = count.out.name, thread = thread,
                                            bam.info = bam.info.dedup)

  if (rmdup == T) return(count.table.dedup)
  else return(count.table)

}

# t.pipeline <- rna.se.qc.align.count(sample.info = "sample_info.tsv", count.out.name = "multi_mini", rmdup = T)
# t.pipeline
