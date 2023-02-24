
# library(jsonlite)
encode.dna.jsongen <- function (fastq.files, other.params.json, assay,
                                master.dir=NULL, run = F, pseudo = T,
                                thread = 3, normalize.path =T) {
  # pseudo: generate json file locations to return without really generating the file itself.
  # assay should be atac or chip
  # generate the json files used for encode.
  # note this script only takes care of the fastq part. other parameters.. you need to put into another json file and
  # read in.

  # it takes in a list of fastq files and works on GSM and rep. replicates must be labeled as "rep1, rep2, ..."
  # the fastq files must be named according to my specifications.

  # I do not accept different biological replicates having mixed pair-end/single-end.
  if (normalize.path == T)
    fastq.files <- normalizePath(fastq.files, mustWork = F)
  if (is.null(master.dir) && !jsonlite::validate(other.params.json))
    master.dir <- other.params.json %>% dirname()

  if (!jsonlite::validate(other.params.json)[1])
    other.json <- jsonlite::read_json(other.params.json)
  other.json <- lapply(other.json, jsonlite::unbox)

  title <- paste0(assay, ".title")
  description <- paste0(assay, ".description")
  fastqs <- paste0(assay, ".fastqs")

  sample.df <-  data.frame(fastq = fastq.files, sample = grab.sample(fastq.files),
                           rep = grab.rep(fastq.files),
                           read = grab.read(fastq.files))
  js.files <- lapply(sample.df %>% split(f = sample.df$sample), function(s) {
      s$read[s$read == "se"] <- "R1"

      s$read.cat <- paste0(fastqs, "_", s$rep, "_",s$read)
      js <- lapply(s %>% split(f=s$read.cat), function(x) x$fastq)

      js[[paste0(assay, ".paired_ends")]] <- lapply(s %>% split(f=s$rep), function(x) {
        if ("R2" %in% s$read)
          return(T)
        else
          return(F)
      }) %>% unlist ()

      js[[title]] <- s$sample[1] %>% jsonlite::unbox()
      js[[description]] <- js[[title]]

      js <- c(other.json, js)
      js <- js %>% jsonlite::toJSON() %>% jsonlite::prettify()
      outdir <- paste0(master.dir, "/", s$sample[1], "/")
      system(paste0("mkdir -p ", outdir))
      js.file <- paste0(outdir, s$sample[1], ".json")

      if (pseudo == F)
      write(js, js.file)

      return(js.file)
  })
 if (pseudo == F && run == T) {
   mclapply(js.files, function (x) {
     path <- read_file("/bar/cfan/PATHs/encode-atac-seq-pipeline.path") %>% sub("\\n", "", .)
     Sys.setenv(PATH=path)
     logfile <- sub("\\.json", ".log", x)
     cmd <- paste0("cd ", dirname(x),
                   " && caper run ", "/bar/cfan/software/atac/encode_pipeline/atac.wdl",
                   " -i ", x, " 1>", logfile, " 2>&1")
     paste0(cmd)
     system(cmd)
   }, mc.cores = thread)
 }
 return(js.files)

}





encode.parser <- function(js.files, assay, bam.dir=NA, unfiltered.bam=NA, ta.dir=NA, peak.dir=NA, idr.dir = NA,
                          qc.dir = NA, bw.dir = NA, frag.length.dir = NA, jsd.dir = NA, gc.dir = NA) {
  # note: when you add more pipe functions, below is the only line that you need to change (add 3 things)

  pipe.functions <- data.frame(filetype = c("filter", "bam2ta", "call_peak", "fraglen_stat_pe", "gc_bias", "align")[!is.na(c(bam.dir, ta.dir, peak.dir, frag.length.dir, gc.dir, unfiltered.bam))],
                          location =c(bam.dir, ta.dir, peak.dir, frag.length.dir, gc.dir, unfiltered.bam)) %>% na.omit()
  print(pipe.functions)

  lapply(pipe.functions$location,
         function (x) {
           print(paste0("mkdir -p ", x))
           system(paste0("mkdir -p ", x))})

  if (!is.null(idr.dir)) system(paste0("mkdir -p ", idr.dir))
  if (!is.null(qc.dir)) system(paste0("mkdir -p ", qc.dir))
  if (!is.null(jsd.dir)) system(paste0("mkdir -p ", jsd.dir))

  lapply(js.files, function (s) {
    sample <- basename(s) %>% sub(".json", "", .)
    rootdir <- paste0(dirname(s), "/atac/*/")
    # link bm files:
    reps <- system(paste0("ls ", rootdir, "call-filter/"), intern = T)

    lapply(pipe.functions %>% split(f=pipe.functions$filetype), function(pf) {
      lapply(seq_along(reps), function(i) {
        cmd <- paste0("rm -rf ", paste0(pf$location, "/", sample, "_rep",i),
                      " && ln -s ", paste0(rootdir, "call-", pf$filetype ,"/shard-",i-1,"/execution"),
                      " ", paste0(pf$location, "/", sample, "_rep",i))
        print(cmd)
        system(cmd)

        if (!is.null(bw.dir) && pf$filetype == "filter")
        {
          bam.name <- system(paste0("ls ", pf$location, "/", sample, "_rep",i, "/*.bam"), intern = T)[1]

          bam.browser(bam = bam.name,
                      bw.dir = bw.dir, thread = 16, normalization = "RPKM")
        }
      })
    })

    if (!is.null(idr.dir)) lapply(system(paste0("ls ", rootdir, "call-idr/"), intern = T), function(x) {
      a <- system(paste0("ls ", rootdir, "/call-idr/", x, "/execution"), intern = T)
      contrast <- a[grepl("rep\\d_vs_rep\\d", a)][1] %>% sub("(rep\\d_vs_rep\\d).+", "\\1", .)
      cmd <- paste0("rm -rf ", paste0(idr.dir, "/", sample, "_", contrast),
                    " && ln -s ", paste0(rootdir, "/call-idr/", x, "/execution"),
                    " ", paste0(idr.dir, "/", sample, "_", contrast))
      print(cmd)
      system(cmd)
    })

    if (!is.null(qc.dir)) {
      qc.report <- paste0(rootdir, "/call-qc_report/execution")
      cmd <- paste0("rm -rf ",paste0(qc.dir, "/", sample),
                    " && ln -s ", qc.report, " ", paste0(qc.dir, "/", sample))
      print(cmd)
      system(cmd)
    }

    if (!is.null(jsd.dir)) {
      jsd.report <- paste0(rootdir, "/call-jsd/execution")
      cmd <- paste0("rm -rf ",paste0(jsd.dir, "/", sample),
                    " && ln -s ", jsd.report, " ", paste0(jsd.dir, "/", sample))
      print(cmd)
      system(cmd)
    }

  })
}


# encode.bw.jsongen <- function(bw.dir.list, out.json) {
#   name <- lapply(bw.dir.list, function(x) Sys.glob(paste0(x, "/*.bw"))) %>% unlist()
#   url <- sub("/bar/cfan/", "https://wangftp.wustl.edu/~cfan/", name)
#   df <- data.frame(name = name, track_type="bw",
#                    url = url, metadata = "novalue", options = "novalue")
#   jsongen <-
#   write.table()
# }



# t <- read.table("~/4dn/nk/fanc/dna_sample_info.tsv", as.is = T, header = T)
#
# tt <- encode.dna.jsongen(fastq.files = t$fastqfile, assay = "atac",
#                          other.params.json = "~/4dn/nk/fanc/lodge/encode/json.template")
#
# ttt <- jsonlite::read_json("~/4dn/nk/fanc/lodge/encode/json.template")

grab.sample <- function(x, is.short.name=F) {
  # short.name: grab something like: "Sp_NK_rep1"
  if (is.short.name==F)
  y <- basename(x) %>% sub("(.+)_rep\\d.+", "\\1", .)
  else
    y <- basename(x) %>% sub("(.+)_rep\\d+$", "\\1", .)
}
grab.rep <- function(x, is.short.name=F) {
  if (is.short.name==F)
    y <- basename(x) %>% sub(".+(rep\\d+).+", "\\1", .)
  else
    y <- basename(x) %>% sub(".+(rep\\d+)$", "\\1", .)
}
grab.read <- function(fastq) {
  lapply(fastq, function (x) {
    if (grepl("SRR\\d+_[12]", x))
      return(basename(x) %>% sub(".+SRR\\d+_([12]).+", "\\1", .) %>% paste0("R", .))
    else
      return("se")
  }) %>% unlist()
}



