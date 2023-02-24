fastq.write <- function(seq.list, score = "F", out.root, out.suffix = ".fastq",
                        zip = T) {
  # format: seq.list has R1 and R2 2 elements. both are named vectors of strings. 
  # each string is a read. read names are the names of the vectors
  if (!identical(sort(names(seq.list)), c("R1", "R2"))) {
    stop("seq.list must be named R1 and R2")
  }
  # seq.list <- seq.list[order(names(seq.list))]
  if (grepl(".gz$", out.suffix)) {
    out.suffix <- out.suffix %>% sub(".gz", "", .)
    zip <- T
  }
  out.files <- lapply(1:2, function(i) {
    seqs <- seq.list[[paste0("R", i)]] %>% toupper()
    quals <- stringr::str_dup(score, nchar(seqs))
    df <- data.frame(rname = paste0("@", names(seqs)), seq = seqs, third = "+", score = quals)
    out.file <- paste0(out.root, i, out.suffix)
    trash <- utilsFanc::write.zip.fanc(df = df, out.file = out.file, zip = F,
                              bed.shift = F, col.names = F, row.names = F, sep = "\n")
    if (zip == T) {
      cmd <- paste0("gzip -f ", out.file)
      system(cmd)
      out.file <- paste0(out.file, ".gz")
    }
    return(out.file)
  })
  names(out.files) <- c("R1", "R2")
  return(out.files)
}

fastq.subsample <- function(in.fastqs, out.dir, n.reads, seed = 123, seqtk = SEQTK,
                            threads = 1, run = T, save.mem = F,
                            suffix.regex = "_R*[12](_001)*.fastq(.gz)*$") {
  dir.create(out.dir, showWarnings = F, recursive = T)
  pairs <- basename(in.fastqs) %>% sub(suffix.regex, "", . )
  n.reads.f <- utilsFanc::so.formatter(n.reads)
  out.fastqs <- in.fastqs %>% split(f = factor(pairs, levels = unique(pairs))) %>% 
    utilsFanc::safelapply(function(fastqs) {
      out.fastqs <- paste0(out.dir, "/", basename(fastqs)) %>% 
        utilsFanc::insert.name.2(insert = paste0("_", n.reads.f), ext = suffix.regex, pre = ".+")
      mem <- ifelse(save.mem, " -2 ", "")
      cmds <- paste0(seqtk, " sample -s ", seed, mem, " ", fastqs, " ", n.reads,
                     " | gzip -nc ", " > ", out.fastqs)
      print(cmds)
      if (run) {
        utilsFanc::safelapply(cmds, system, threads = 2)
        lapply(out.fastqs,function(x) {
          if (!file.exists(x)) {
            stop(paste0(x, " failed to generate"))
          }
        })
      }
      return(out.fastqs)
    },threads = threads) %>% unlist()
  return(out.fastqs)
}