#
cutadapt.target <- function(in.fastq.vec, out.fastq.dir = NULL,
                            adapt1 = ATAC.ADAPTER.1, adapt2 = ATAC.ADAPTER.2,
                            cutadapt = "/bar/cfan/anaconda2/envs/jupyter/bin/cutadapt",
                            params = CUTADAPT.PARAMS.TARGET, threads.master = 4, 
                            threads.sub = 4, run = T) {
  if (is.null(out.fastq.dir))
    out.fastq.dir <- in.fastq.vec[1] %>% dirname() %>% paste0("/cutadapt/")
  log.dir <- paste0(out.fastq.dir, "/cutadapt.logs/")
  dir.create(out.fastq.dir, showWarnings = F, recursive = T)
  dir.create(log.dir, showWarnings = F, recursive = T)
  # automatically group reads based on R1 and R2. 
  
  fastq.groups <- in.fastq.vec %>% sub("_[12].fastq", "",.) %>% 
    sub("_R[12].001.fastq", "", .) %>% 
    sub("_R[12].fastq", "", .)
  out.fastqs <- in.fastq.vec %>% split(f = fastq.groups) %>% 
    mclapply(function(fastqs) {
      out.files <- paste0(out.fastq.dir, "/", basename(fastqs)) %>% 
        utilsFanc::insert.name.2(insert = "cutadapt_", ext = "R[12](.001)*.fastq") %>% 
        utilsFanc::insert.name.2(insert = "cutadapt_", ext = "[12].fastq", pre = "_")
      if (length(fastqs) == 1) {
        cat("single end\n")
        cat(fastqs)
        cat("\n")
        cmd <- paste0(cutadapt, " -a ", adapt1, " -j ", threads.sub, 
                      " ", params, " -o ", out.files[1], 
                      " ", fastqs[1])
        
      } else if (length(fastqs) == 2) {
        print("pairend")
        cat(fastqs)
        cat("\n")
        cmd <- paste0(cutadapt, " -a ", adapt1, " -A ", adapt2, " -j ", threads.sub, 
                      " ", params, " -o ", out.files[1], " -p ", out.files[2],
                      " ", fastqs[1], " ", fastqs[2])
      } else {
        stop(paste0(paste0(fastqs, collapse = " "), " is not single or paired end"))
      }
      
      utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = paste0(log.dir, "/", basename(out.files[1]), ".log"),
                               intern = F, run = run)
      for (out.file in out.files) {
        if (!file.exists(out.file))
          stop(paste0(out.file, " failed to generate"))
      }
      return(out.files)
    }, mc.cores = threads.master, mc.cleanup = T)
 
  return(out.fastqs)
}
bwa.mem <- function(fastq.vec, ref.fa, out.bam, threads,
                    samtools = SAMTOOLS,
                    bwa = BWA.TARGET,
                    mapped.only = F,
                    add.params = "",
                    log.file = NULL) {
  dir.create(dirname(out.bam), showWarnings = F,recursive = T)
  cmd <- paste0(bwa, " mem -t ", threads, " ", add.params, " ", ref.fa, " ", paste0(fastq.vec, collapse = " "))
  if (!is.null(log.file))
    cmd <- paste0(cmd, " 2>",log.file)
  if (mapped.only == T)
    cmd <- paste0(cmd, " | ", samtools, " view -h -@ ", threads, " -F 0x4 - ")
    # cmd <- paste0(cmd , " | ", sambamba, " view -t ", threads, " -f SAM -F ", "\"not (unmapped)\"", " /dev/stdin ")
  cmd <- paste0(cmd, " | ", samtools, " sort - -O bam -o ", out.bam, " -@ ", threads)
  utilsFanc::cmd.exec.fanc(cmd, stdout.file = NULL, run = T, intern = F)
  if (!file.exists(out.bam))
    stop(paste0(out.bam, " failed to generate"))
  cmd <- paste0(samtools, " index ", out.bam)
  utilsFanc::cmd.exec.fanc(cmd, run = T, intern = F)
  return(out.bam)
}

bowtie2.wrapper <- function(fastq.or.pair, genome.index, align.out.dir,
                            out.bam = NULL, threads = 6, 
                            a = F, k = 10, X = 10000, mm = F, report.unaligned = F,
                            preset = " --very-sensitive ", defaults = BOWTIE2.DEFAULTS,
                            add.params = "", run = T, log.file = NULL,
                            bowtie2 = BOWTIE2, samtools = SAMTOOLS) {
  if (is.null(out.bam)) {
    root.name <- trim.fastq(fastq.or.pair[1], T)
    align.out.dir <- normalizePath(align.out.dir, mustWork = F)
    out.bam <- paste0(align.out.dir, "/", root.name, ".bam")
  }
  
  if (length(fastq.or.pair) == 1) 
    cmd.fastq <- paste0(" -U ", fastq.or.pair)
  else if (length(fastq.or.pair) == 2) 
    cmd.fastq <- paste0(" -1 ", fastq.or.pair[1], " -2 ", fastq.or.pair[2])
  else stop("must provide one fastq or a pair of fastqs")
  
  if (a == T) {
    report <- " -a "
  } else if (!is.null(k)) {
    report <- paste0(" -k ", k)
  } else {
    report <- ""
  }
  
  if (mm == T) {
    mm <- " --mm "
  } else {
    mm <- ""
  }
  
  if (report.unaligned == F) {
    unal <- " --no-unal "
  } else {
    unal <- ""
  }
  dir.create(dirname(out.bam), showWarnings = F, recursive = T)
  cmd <- paste0(bowtie2, " -x ", genome.index, " ", cmd.fastq, " ", preset,
                " ", defaults, " ", report, " -p ", threads, " ", mm, " ", unal,
                " -X ", X, " ",add.params)
  if (!is.null(log.file)) {
    dir.create(dirname(log.file), showWarnings = F, recursive = T)
    cmd <- paste0(cmd, " 2>", log.file, " ")
  }
   
  
  cmd <- paste0(cmd, " | ", samtools, " sort - -O bam -m 2G -o ", out.bam, " -@ ", threads)
  
  utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = NULL, intern = F, run = run)
  
  system(paste0(samtools, " index ", out.bam))
  if (!file.exists(out.bam))
    stop(paste0(out.bam, " failed to generate"))
  else
    return(out.bam)
}

bowtie2.fast.realign <- function(in.bam.df, tmp.dir = tempdir(), 
                                 genome.index, 
                                 bowtie2 = BOWTIE2, ...) {
  
}

bowtie2.index <- function(fa, root.name = NULL, out.dir = NULL, threads,
                          bowtie2.build = "/opt/apps/bowtie2/2.3.4.1/bowtie2-build",
                          log.file = NULL) {
 
  if(is.null(root.name))
    root.name <- basename(fa) %>% sub(".fn*a.*$", "", .)
  if (is.null(out.dir))
    out.dir <- paste0(dirname(fa), "/bowtie2/")
  dir.create(out.dir, showWarnings = F, recursive = T)
  cmd <- paste0(bowtie2.build, " --threads ", threads, " ", fa, 
                " ", out.dir, "/", root.name)
  utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = log.file, intern = F, run = T)
  return(paste0(out.dir, "/", root.name))
}

bowtie2.align.pe.se <- function(fastq.or.pair, genome.index, thread, bowtie.out.bam, filter.expression=NULL) {
  if (length(fastq.or.pair) == 1)
    fastq.pair <- c(fastq.or.pair, "se")
  else
    if (length(fastq.or.pair) ==2 ) fastq.pair <- fastq.or.pair
    else stop("must provide one fastq or a pair of fastqs")
  cmd <- paste0("~/scripts/dna-seq/bowtie2_pe_se.sh ", " -p ", thread, " -x ", genome.index, " -i ", fastq.pair[1],
                     " -I ", fastq.pair[2], " -o ", bowtie.out.bam)
  if (!is.null(filter.expression)) cmd <- paste0(cmd, " -f ",'"', filter.expression, '"')
  system(cmd)
  if(file.exists(bowtie.out.bam)) return(bowtie.out.bam)
  else stop(paste0("bam file was not successfully generated for: ", fastq.pair[1]))
}

bowtie2.align.pe.se.2 <- function(fastq.or.pair, genome.index, thread, 
                                  align.out.dir, filter.expression=NULL, pseudo = F,
                                  bowtie2.script = "~/scripts/dna-seq/bowtie2_pe_se.sh") {
  root.name <- trim.fastq(fastq.or.pair[1], T) %>% pair.end.rename()
  align.out.dir <- normalizePath(align.out.dir, mustWork = F)
  bowtie.out.bam <- paste0(align.out.dir, "/", root.name, ".bam")
  bowtie.log <- sub(".bam$", ".log", bowtie.out.bam)
  if (length(fastq.or.pair) == 1)
    fastq.pair <- c(fastq.or.pair, "se")
  else
    if (length(fastq.or.pair) ==2 ) fastq.pair <- fastq.or.pair
  else stop("must provide one fastq or a pair of fastqs")
  cmd <- paste0("bash ", bowtie2.script, " ", " -p ", thread, " -x ", genome.index, " -i ", fastq.pair[1],
                " -I ", fastq.pair[2], " -o ", bowtie.out.bam)
  if (!is.null(filter.expression)) cmd <- paste0(cmd, " -f ",'"', filter.expression, '"')
  utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = bowtie.log, intern = F, run = !pseudo)

  if(file.exists(bowtie.out.bam) || pseudo == T) return(bowtie.out.bam)
  else stop(paste0("bam file was not successfully generated for: ", fastq.pair[1]))
}

# bowtie2.align.se <- function(fastq, genome.index, thread, bowtie.out.bam) {
#   system(paste0("~/scripts/dna-seq/bowtie2_pe_se.sh ", " -p ", thread, " -x ", genome.index, " -i ", fastq, " -o ", bowtie.out.bam))
#   if(file.exists(bowtie.out.bam)) return(bowtie.out.bam)
#   else stop(paste0("bam file was not successfully generated for: ", fastq))
# }


# bam.filter <- function(in.bam, out.bam=NULL, creat.index = F, filter.expression, thread,
#   samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
#   if (is.null(out.bam))
#     out.bam <- insert.name(in.bam, insert = "filtered", ext = ".bam", trim.dir=F)

#   filter.cmd <- paste0("~/software/sambamba/sambamba view -f bam -F ", '"',filter.expression, '"', " -t ", thread, " ", in.bam, " > ", out.bam )
#   print(filter.cmd)
#   system(filter.cmd)
#   if (creat.index == T) {
#     cmd <- paste0(samtools, " index ", out.bam)
#     print(cmd)
#     system(cmd)
#   }

#   if (file.exists(out.bam)) return(out.bam)
#   else (stop(paste0("bam file was not successfully filtered for: ", in.bam)))
#   return(out.bam)
# }

bam.filter <- function (in.bam, out.bam = NULL, creat.index = F, create.bw = F,
                        filter.expression, remove.mate=F, debug = F,
                        thread, samtools = SAMTOOLS,
                        sambamba = SAMBAMBA)
{
  if (is.null(out.bam))
    out.bam <- insert.name(in.bam, insert = "filtered", ext = ".bam",
                           trim.dir = F)
  out.bam.tmpt <- insert.name(out.bam, insert = "tempt", ext = ".bam",
                                         trim.dir = F)


  filter.cmd <- paste0(sambamba, " view -f bam -F ",
                       "\"", filter.expression, "\"", " -t ", thread, " ", in.bam,
                       " > ", out.bam.tmpt)
  print(filter.cmd)
  system(filter.cmd)

  if (remove.mate == T) {
    bam.nsort <- insert.name(name = out.bam.tmpt, insert = "nsort", ext = ".bam", trim.dir = F)
    cmd <- paste0(samtools, " sort -n ", " -@ ", thread,
                  " -o ", bam.nsort,
                  " ", out.bam.tmpt)
    print(cmd)
    system(cmd)

    bam.nsort.fixmate <- insert.name(name = bam.nsort, insert = "fixmate", ext = ".bam", trim.dir = F)
    cmd <- paste0(samtools, " fixmate ", bam.nsort, " ", bam.nsort.fixmate)
    print(cmd)
    system(cmd)

    cmd <- paste0(samtools, " sort ", " -@ ", thread, " ", bam.nsort.fixmate, " | ",
                  sambamba, " view -f bam -F ",
                  "\"", "paired", "\"", " -t ", thread, " ", "/dev/stdin",
                  " > ", out.bam.tmpt)
    print(cmd)
    system(cmd)

    if (debug == F) {
      system(paste0("rm ", bam.nsort, " ", bam.nsort.fixmate))
    }
  }

  cmd <- paste0("mv ", out.bam.tmpt, " ", out.bam)
  print(cmd)
  system(cmd)

  if (creat.index == T) {
    cmd <- paste0(samtools, " index ", out.bam)
    print(cmd)
    system(cmd)
  }
  if (create.bw == T) {
    bam.browser(bam = out.bam, bw.dir = dirname(out.bam), thread = thread)
  }

  if (file.exists(out.bam))
    return(out.bam)
  else (stop(paste0("bam file was not successfully filtered for: ",
                    in.bam)))
  return(out.bam)
}


bam.dedup <- function(in.bam, out.bam=NULL, thread, other.options ="", remove=T,
                      sambamba = SAMBAMBA) {
  if(is.null(out.bam) && remove == T)
    out.bam <- insert.name(in.bam, insert = "dedup", ext = ".bam", trim.dir=F)
  else if (is.null(out.bam) && remove == F)
    out.bam <- insert.name(in.bam, insert = "mkdup", ext = ".bam", trim.dir=F)

  if (remove == T)
    remove = "-r"
  else remove = ""
  cmd <- paste0(sambamba, " markdup -t ", thread, " ", remove, " ", other.options, " ", in.bam, " ", out.bam)
  system(cmd)
  if(file.exists(out.bam)) return(out.bam)
  else stop(paste0("bam file was not successfully dedupped for: ", in.bam))
  return(out.bam)
}

bam.browser <- function(bam, bw.dir=NULL, thread, normalization="RPKM", other.bw.options = "", stranded = 0,
  bamCoverage.path = BAMCOVERAGE) {
  if (is.null(bw.dir)) {
    bw.dir <- dirname(bam)
  }
  paste0("mkdir -p ", bw.dir) %>% system()
  if (!is.null(normalization)) normalization <- paste0(" --normalizeUsing ", normalization)
  if (!file.exists(paste0(bam, ".bai"))) {
      cmd <- paste0("samtools index -b ", bam)
      system(cmd)
  }

  bw.options = paste0(" -b ", bam, " -p ", thread, " ", normalization, " ", other.bw.options)

  if (stranded == 0) {
    bw <- paste0(bw.dir, "/", trim.bam(bam, trim_dir=T), ".bw")
    cmd <- paste0(bamCoverage.path, " -o ", bw , bw.options)
    print(cmd)
    system(cmd)
  } else {
    bw <- paste0(bw.dir, "/", trim.bam(bam, trim_dir=T), ".stranded_", stranded,".bedgraph")
    bw1 <- paste0(bw.dir, "/", trim.bam(bam, trim_dir=T), ".f.bedgraph")
    cmd1 <- paste0(bamCoverage.path, " -of bedgraph --filterRNAstrand forward -o ", bw1,  bw.options)
    print(cmd1)
    system(cmd1)
    bw2 <- paste0(bw.dir, "/", trim.bam(bam, trim_dir=T), ".r.bedgraph")
    cmd2 <- paste0(bamCoverage.path, " -of bedgraph --filterRNAstrand reverse -o ", bw2, bw.options)
    print(cmd2)
    system(cmd2)
    if (stranded == 1) {
      cmd3 <- paste0("~/scripts/stranded_bedgraph.sh ",bw1, " ", bw2, " ", bw)
    } else if (stranded == 2) {
      cmd3 <- paste0("~/scripts/stranded_bedgraph.sh ",bw2, " ", bw1, " ", bw)
    } else {
      stop("error in stranded bedgraph generation: 'stranded' parameter has to be 0, 1, or 2.")
    }
    system(cmd3)
    bw <- paste0(bw, ".sort.gz")
  }

  if(file.exists(bw)) return(bw)
  else stop(paste0("bw file was not successfully generated for: ", bam))
}

grab.SRR <- function(fastq) {
  sub(".+(SRR\\d+).+", "\\1", fastq) %>%
  return()
}
grab.GSM <- function(fastq) {
  sub(".+(GSM\\d+).+", "\\1", fastq) %>%
    return()
}

pair.end.rename <- function (fastq) {
  sub("(.+SRR\\d+)(_[12])(.+)", "\\1\\3", fastq) %>%
  return()
}

# fastq.tag.to.bam <- function (fastq) {
#   bam <- sub(".fastq.gz", ".bam", fastq)
#   return(bam)
# }

dna.seq <- function(fastq.files, genome.index, qc=T, pe.rmdup = T, se.rmdup = F,
                    filter.expression = "mapping_quality >= 30 and (not secondary_alignment)",
                    fastqc.dir="./qc", fastp.out.dir="./", bw.dir = "./bigwig/" ,thread,
                    align.out.dir="./bowtie2_aligned/", other.dedup.options = " --hash-table-size=5000000 --overflow-list-size=5000000 ",
                    normalization="RPKM", other.bw.options = "", mc.cores = 4) {
# fastq.files is a vector. pair end and single end can be mised. genome.index should be only 1
  options(stringsAsFactors = F)
  lapply(list(fastqc.dir, fastp.out.dir, bw.dir, align.out.dir), function (dir) {
    system(paste0("mkdir -p ",dir))
  })
  fastq.df <- data.frame(fastq = fastq.files, SRR = grab.SRR(fastq.files))
  print(fastq.df)
  fastq.by.SRR <- split(fastq.df, f=fastq.df$SRR)
  utilsFanc::safelapply(fastq.by.SRR, function (df) {
    print(df)
    fastq <- df$fastq

    if (!length(fastq) %in% c(1,2))
      stop(paste0("something is wrong with ", df$SRR[1], ": there is no files or more than 2 files associated with it."))
    if (qc == T) {
      if (length(fastq) == 1)
        fastq <- qc.se(fastq, fastqc.dir = fastqc.dir, fastp.out.dir = fastp.out.dir, thread = thread)
      else
        fastq <- qc.pe(fastq, fastqc.dir = fastqc.dir, fastp.out.dir = fastp.out.dir, thread = thread)
    }
    root.name <- trim.fastq(fastq[1], T) %>% pair.end.rename()
    bowtie.out.bam <- paste0(align.out.dir,"/",root.name, ".bam")
    if (!is.null(filter.expression)) bowtie.out.bam <- insert.name(bowtie.out.bam, insert = "filtered", ext = ".bam", trim.dir=F)

    if (length(fastq) == 1) rmdup <- se.rmdup
    else rmdup <- pe.rmdup

    if (rmdup == T) dedup.out.bam <- insert.name(bowtie.out.bam, insert = "dedup", ext = ".bam", trim.dir=F)
    else dedup.out.bam <- insert.name(bowtie.out.bam, insert = "mkdup", ext = ".bam", trim.dir=F)

    fastq %>% bowtie2.align.pe.se(genome.index = genome.index, thread = thread,
                                  bowtie.out.bam = bowtie.out.bam, filter.expression = filter.expression) %>%
      bam.dedup(out.bam = dedup.out.bam, thread = thread, remove = rmdup, other.options = other.dedup.options) %>%
      bam.browser(bw.dir=bw.dir, thread=thread, normalization=normalization, other.bw.options = other.bw.options) %>%
      return()
  }, threads = mc.cores) %>% unlist() %>%
  return()
}


