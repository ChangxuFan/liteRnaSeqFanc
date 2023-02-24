dna.rna.seq <- function(fastq.info, type, genome.index, annot, qc=T, pe.rmdup = T, se.rmdup = F,
                        skip.qc = F, skip.align = F, skip.dedup=F,generate.bw=T,
                        bowtie2.script = "~/scripts/dna-seq/bowtie2_pe_se.sh",
                        filter.expression.se = NULL, filter.expression.pe=NULL, count.out.name, count.frac.overlap,
                        other.count.options.pe, other.count.options.se,
                        baminfo.file = "./bam_info.tsv",
                        fastqc.dir="./qc", fastp.out.dir="./fastp", bw.dir = "./bw/" ,
                        align.out.dir="./aligned/", other.dedup.options = " --hash-table-size=3000000 --overflow-list-size=3000000 ",
                        normalization="RPKM", other.bw.options = "",
                        thread, mc.cores.GSM = 4, mc.cores.SRR = 2) {
  if (type == "rna")
    stop("the rna part has not been tested")

  # first create folders for all the output data
  lapply(list(fastqc.dir, fastp.out.dir, bw.dir, align.out.dir), function (dir) {
    system(paste0("mkdir -p ",dir))
  })
  # then read in fastq.info as a dataframe or external file that can be read in to be a dataframe. This part is adopted from
  # my single end pipeline mainly because I need to later map bam file names to sample names after featureCounts

  if (is.character(fastq.info))
    fastq.info <- read.table(fastq.info, as.is = T, header = T, sep="\t")

  fastq.gsm <- data.frame(fastq= fastq.info$fastqfile, sample =  fastq.info$sample, GSM = grab.GSM(fastq.info$fastqfile))

  fastq.by.gsm <- split(fastq.gsm, f=fastq.gsm$GSM)
  print(fastq.by.gsm)
  bam.by.gsm <- mclapply(fastq.by.gsm, function(x) {
    fastq.df <- data.frame(fastq = x$fastq, SRR = grab.SRR(x$fa))
    fastq.by.SRR <- split(fastq.df, f=fastq.df$SRR)
    #print(fastq.by.SRR)

    # note: now I tell se/pe at gsm level. Under the assumption that for a fiven GSM, it's not supposed to have both se and pe.
    # because one GSM is defined by one SRX (usually, at least the way I handle it.)
    # logic: every SRR should correspond to 1 file or 2 files.
    print(sapply(fastq.by.SRR, nrow))
    if (sum(sapply(fastq.by.SRR, nrow) != 2) == 0) {
      sepe = "pe"
    } else if (sum(sapply(fastq.by.SRR, nrow) != 1) == 0) {
      sepe = "se"
    } else stop ("Error: each GSM should contain only pairend SRRs or only single end SRRs")
    print(sepe)
    bam.by.SRR <- mclapply(fastq.by.SRR, function (df) {
      print(df)
      fastq <- df$fastq

      if (!length(fastq) %in% c(1,2))
        stop(paste0("something is wrong with ", df$SRR[1], ": there is no files or more than 2 files associated with it."))

      if (qc == T) {
        if (skip.qc == T) {
          fastq <- sapply(fastq, function (fastq) {
            paste0(fastp.out.dir, "/", trim.fastq(fastq = fastq, trim_dir = T), "_fastp_fastq.gz") %>% normalizePath(mustWork = F) %>%
              return()
          })

        } else {
          if (sepe == "se")
            fastq <- qc.se(fastq, fastqc.dir = fastqc.dir, fastp.out.dir = fastp.out.dir, thread = thread)
          else
            fastq <- qc.pe(fastq, fastqc.dir = fastqc.dir, fastp.out.dir = fastp.out.dir, thread = thread)
        }
      }

      if (sepe == "se") {
        rmdup <- se.rmdup
        filter.expression <- filter.expression.se
      }
      else {
        rmdup <- pe.rmdup
        filter.expression <- filter.expression.pe
      }

      if ( type == "rna") {
        bam <- fastq %>% star.align.pe.se.2(genome.index = genome.index, thread = thread,
                                          align.out.dir = align.out.dir, sepe = sepe, pseudo = skip.align)
      }
      if ( type == "dna") {
        bam <- fastq %>% bowtie2.align.pe.se.2(genome.index = genome.index, thread = thread, 
                                               bowtie2.script = bowtie2.script,
                                             filter.expression = NULL, align.out.dir = align.out.dir, pseudo = skip.align)
      }

      # here rna and dna behaves differently. For RNA: each SRR in one directory, because star generates too much random stuff.
      # for bowtie2, all samples will be in the same dir to make downstream stuff easier.
      if (!is.null(filter.expression)) 
        bam <- bam.filter(in.bam = bam, filter.expression = filter.expression, thread = thread, creat.index = T)
      
      if (skip.dedup == F) {
        bam <- bam.dedup(in.bam = bam,thread = thread, remove = rmdup, other.options = other.dedup.options)
      }

      return(bam)

    }, mc.cores = mc.cores.SRR) %>% unlist()

    if (length(bam.by.SRR) > 1) {
      bam <- bam.merge.gsm(in.bam.vector = bam.by.SRR, align.out.dir = align.out.dir,thread = thread)
    } else {
      bam <- bam.by.SRR
    }

    if (generate.bw == T)
      try(bam.browser(bam = bam, bw.dir=bw.dir, thread=thread, normalization=normalization, other.bw.options = other.bw.options))
    # prepare bam info for featureCounts
    bam.info <- data.frame(bamfile = bam, sample=x$sample[1], sepe = sepe)
    return(bam.info)
  }, mc.cores = mc.cores.GSM)

  bam.by.gsm <- Reduce(rbind, bam.by.gsm)

  write.table(bam.by.gsm, baminfo.file, quote = F, col.names = T, row.names = F, sep = "\t")
  json.df <- data.frame(name = basename(bam.by.gsm$bam),
                        url = utilsFanc::bash2ftp(bam.by.gsm$bam),
                        track_type = "bam",
                        metadata = "novalue",
                        options = "novalue")

  trash <- utilsFanc::jsongen(df = json.df, outfile = paste0(baminfo.file, ".json"))

  if (generate.bw == T) {
    system(paste0("Rscript --vanilla ~/R_for_bash/json_dir.R ", bw.dir, " bigwig"))
  }

  if (type == "RNA") {
    # split into single end and pair end for featureCounts:
    bam.by.gsm.se <- bam.by.gsm %>% filter(sepe == "se")
    if (nrow(bam.by.gsm.se) > 0)
      count.se <- rna.se.pe.featureCounts(bam.info = bam.by.gsm.se, annot = annot, count.out.name = paste0(count.out.name, ".se"),
                                          thread = thread, count.frac.overlap = count.frac.overlap,
                                          other.count.options = other.count.options.se)
    else count.se <- "no single end files"
    bam.by.gsm.pe <- bam.by.gsm %>% filter(sepe == "pe")
    if (nrow(bam.by.gsm.pe) > 0)
      count.pe <- rna.se.pe.featureCounts(bam.info = bam.by.gsm.pe, annot = annot, count.out.name = paste0(count.out.name, ".pe"),
                                          thread = thread, count.frac.overlap = count.frac.overlap,
                                          other.count.options = other.count.options.pe)
    else count.pe <- "no pair end files"
    return(list(count.se = count.se, count.pe = count.pe))
  }

  return()

}
