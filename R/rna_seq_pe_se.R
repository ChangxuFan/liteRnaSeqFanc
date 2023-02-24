#
rna.se.pe.featureCounts <- function(bam.info, annot, count.out.name=NULL, thread,
                                    field = "gene_name", 
                                    count.frac.overlap, other.count.options) {
  # if bam is null: file names of the bam files will be read from bam.info
  #bam.info: a sample name ~ bamfile name correspondance file. tsv format.
  #must have a colum named "sample" and another column named "bamfile"
  #bam.info can also be a dataframe mapping bam file names to sample names.
  if (grepl("vM26", annot)) {
    stop("gencode vM26 doesn't work well with featureCounts. It gives you
    very low assignment rates for unknown reason")
  }
  if (is.character(bam.info))
    bam.info <- read.table(bam.info, as.is = T, header = T, sep="\t")
  bam.expand <- paste0(bam.info$bamfile, collapse = " ")
  count.cmd <- paste0("featureCounts ", bam.expand, " -a ", annot, " -F GTF -g ",field,
                      " -o ", count.out.name, " ",
                      " -t exon ", " -T ", thread, " --fracOverlap ", count.frac.overlap, " ", other.count.options)
  print(count.cmd)
  system(count.cmd)
  count.formatted <- count.table.change.name(count.table = count.out.name, bam.info = bam.info)
  if(file.exists(count.formatted)) return(count.formatted)
  else stop("error: count table was not successfully generated")

}


bam.merge.gsm <- function(in.bam.vector, out.bam=NULL, align.out.dir=NULL, thread) {
  if (is.null(out.bam)) {
    srr <- grab.SRR(in.bam.vector) %>% unique() %>% paste0(collapse = "_")
    dedup <- ifelse(T %in% grepl("dedup", in.bam.vector), "dedup", "")
    mkdup <- ifelse(T %in% grepl("mkdup", in.bam.vector), "mkdup", "")
    dup <- paste0(dedup,"_", mkdup)
    dup <- sub("_$", "", dup)
    root.name <- basename(in.bam.vector[1])

    root.name <- sub("SRR\\d+", srr, root.name)
    root.name <- sub("[md][ek]dup", dup, root.name)
    if (! is.null(align.out.dir)) out.bam <- paste0(align.out.dir, "/", root.name)
    else
      stop("error in bam merge: either bam file name for output directory should be offered ")
  }
  merge.cmd <- paste0("samtools merge -f -@ ", thread, " ", out.bam, " ", paste0(in.bam.vector, collapse = " "))
  print(merge.cmd)
  system(merge.cmd)
  if(!file.exists(out.bam))
    stop(paste0("error [merging bam files]: ", out.bam, " failed to generate"))

  return(out.bam)
}

star.align.pe.se <- function(fastq, genome.index, align.out.dir, thread, sepe) {
  if (sepe == "se") fastq[2] <- "single"
  root.name <- trim.fastq(fastq[1], T) %>% pair.end.rename()
  align.out.dir <- paste0(normalizePath(align.out.dir, mustWork = F), "/", root.name)
  star.out.bam <- paste0(align.out.dir,"/",root.name, "_star_genome.bam")
  align.cmd <- paste0("~/scripts/rna-seq/encode_se_pe_star.sh ",genome.index,
                      " ", paste0(fastq, collapse= " "), " ",root.name, " ", thread, " ",  align.out.dir, " ", root.name )
  print(align.cmd)
  system(align.cmd)

  if(!file.exists(star.out.bam)) stop(paste0("error: ",star.out.bam, " failed to generate"))
  return(star.out.bam)
}



star.align.pe.se.2 <- function(fastq, genome.index, align.out.dir, thread, sepe, pseudo=F) {
  if (sepe == "se") fastq[2] <- "single"
  root.name <- trim.fastq(fastq[1], T) %>% pair.end.rename()
  align.out.dir <- paste0(normalizePath(align.out.dir, mustWork = F), "/", root.name)
  star.out.bam <- paste0(align.out.dir,"/",root.name, "_star_genome.bam")
  align.cmd <- paste0("~/scripts/rna-seq/encode_se_pe_star.sh ",genome.index,
                      " ", paste0(fastq, collapse= " "), " ",root.name, " ", thread, " ",  align.out.dir, " ", root.name )
  utilsFanc::cmd.exec.fanc(cmd = align.cmd, run = !pseudo)

  if(!file.exists(star.out.bam) && pseudo == F) stop(paste0("error: ",star.out.bam, " failed to generate"))
  return(star.out.bam)
}


# rna.seq <- function(fastq.info, genome.index, annot, qc=T, pe.rmdup = T, se.rmdup = F,
#                     filter.expression.se, filter.expression.pe, count.out.name, count.frac.overlap,
#                     other.count.options.pe, other.count.options.se,
#                     fastqc.dir="./qc", fastp.out.dir="./", generate.bw=T, bw.dir = "./bigwig/" ,thread,
#                     align.out.dir="./star2_aligned/", other.dedup.options = " --hash-table-size=3000000 --overflow-list-size=3000000 ",
#                     normalization="RPKM", other.bw.options = "", mc.cores.GSM = 4, mc.cores.SRR = 2) {
#   # first create folders for all the output data
#   lapply(list(fastqc.dir, fastp.out.dir, bw.dir, align.out.dir), function (dir) {
#     system(paste0("mkdir -p ",dir))
#   })
#   # then read in fastq.info as a dataframe or external file that can be read in to be a dataframe. This part is adopted from
#   # my single end pipeline mainly because I need to later map bam file names to sample names after featureCounts
# 
#   if (is.character(fastq.info))
#     fastq.info <- read.table(fastq.info, as.is = T, header = T, sep="\t")
# 
#   fastq.gsm <- data.frame(fastq= fastq.info$fastqfile, sample =  fastq.info$sample, GSM = grab.GSM(fastq.info$fastqfile))
# 
#   fastq.by.gsm <- split(fastq.gsm, f=fastq.gsm$GSM)
#   print(fastq.by.gsm)
#   bam.by.gsm <- mclapply(fastq.by.gsm, function(x) {
#     fastq.df <- data.frame(fastq = x$fastq, SRR = grab.SRR(x$fa))
#     fastq.by.SRR <- split(fastq.df, f=fastq.df$SRR)
#     #print(fastq.by.SRR)
# 
#     # note: now I tell se/pe at gsm level. Under the assumption that for a fiven GSM, it's not supposed to have both se and pe.
#     # because one GSM is defined by one SRX (usually, at least the way I handle it.)
#     # logic: every SRR should correspond to 1 file or 2 files.
#     print(sapply(fastq.by.SRR, nrow))
#     if (sum(sapply(fastq.by.SRR, nrow) != 2) == 0) {
#       sepe = "pe"
#     } else if (sum(sapply(fastq.by.SRR, nrow) != 1) == 0) {
#       sepe = "se"
#     } else stop ("Error: each GSM should contain only pairend SRRs or only single end SRRs")
#     print(sepe)
#     bam.by.SRR <- mclapply(fastq.by.SRR, function (df) {
#       print(df)
#       fastq <- df$fastq
# 
#       if (!length(fastq) %in% c(1,2))
#         stop(paste0("something is wrong with ", df$SRR[1], ": there is no files or more than 2 files associated with it."))
#       if (qc == T) {
#         if (sepe == "se")
#           fastq <- qc.se(fastq, fastqc.dir = fastqc.dir, fastp.out.dir = fastp.out.dir, thread = thread)
#         else
#           fastq <- qc.pe(fastq, fastqc.dir = fastqc.dir, fastp.out.dir = fastp.out.dir, thread = thread)
#       }
#       if (sepe == "se") {
#         rmdup <- se.rmdup
#         filter.expression <- filter.expression.se
#       }
#       else {
#         rmdup <- pe.rmdup
#         filter.expression <- filter.expression.pe
#       }
# 
#       bam <- fastq %>% star.align.pe.se(genome.index = genome.index, thread = thread, align.out.dir = align.out.dir, sepe = sepe)
# 
#       if (!is.null(filter.expression)) bam <- bam.filter(in.bam = bam, filter.expression = filter.expression, thread = thread)
# 
# 
#       bam <- bam.dedup(in.bam = bam,thread = thread, remove = rmdup, other.options = other.dedup.options)
# 
#       return(bam)
# 
#     }, mc.cores = mc.cores.SRR) %>% unlist()
#     bam <- bam.merge.gsm(in.bam.vector = bam.by.SRR, align.out.dir = align.out.dir,thread = thread)
#     if (generate.bw == T)
#       try(bam.browser(bam = bam, bw.dir=bw.dir, thread=thread, normalization=normalization, other.bw.options = other.bw.options))
#     # prepare bam info for featureCounts
#     bam.info <- data.frame(bamfile = bam, sample=x$sample[1], sepe = sepe)
#     return(bam.info)
#   }, mc.cores = mc.cores.GSM)
#   bam.by.gsm <- Reduce(rbind, bam.by.gsm)
#   # split into single end and pair end for featureCounts:
#   bam.by.gsm.se <- bam.by.gsm %>% filter(sepe == "se")
#   if (nrow(bam.by.gsm.se) > 0)
#     count.se <- rna.se.pe.featureCounts(bam.info = bam.by.gsm.se, annot = annot, count.out.name = paste0(count.out.name, ".se"),
#                                         thread = thread, count.frac.overlap = count.frac.overlap,
#                                         other.count.options = other.count.options.se)
#   else count.se <- "no single end files"
#   bam.by.gsm.pe <- bam.by.gsm %>% filter(sepe == "pe")
#   if (nrow(bam.by.gsm.pe) > 0)
#     count.pe <- rna.se.pe.featureCounts(bam.info = bam.by.gsm.pe, annot = annot, count.out.name = paste0(count.out.name, ".pe"),
#                                         thread = thread, count.frac.overlap = count.frac.overlap,
#                                         other.count.options = other.count.options.pe)
#   else count.pe <- "no pair end files"
#   return(list(count.se = count.se, count.pe = count.pe))
# 
# }

star.parse.log <- function(log.vec, out.file=NULL) {
  df <- lapply(seq_along(log.vec), function(i) {
    tf <- tempfile()
    log <- log.vec[i]
    log.name <- names(log.vec)[i]
    system(paste0("grep -v \":$\" ", log, " > ", tf))
    df <- read.table(tf, as.is = T, header = F, sep = "\t", quote = "")
    df$V1 <- df$V1 %>% sub("^ +", "", .) %>% sub(" +\\|", "",.)
    df$V1[7:20] <- paste0("UNIQUE READS::", df$V1[7:20])
    df$V1[21:24] <- paste0("MULTI-MAPPING READS::", df$V1[21:24])
    df$V1[25:27] <- paste0("UNMAPPED READS::", df$V1[25:27])
    df$V1[28:29] <- paste0("CHIMERIC READS::", df$V1[28:29])
    names(df) <- c("qc", log.name)
    return(df)
  }) %>% Reduce(left_join, .)
  if (!is.null(out.file)) {
    write.table(df, out.file, sep = "\t", row.names = F, col.names = T, quote = F)
  }
  return(df)
}

# star.index.gen <- function() {
#   
# }
