bam2ta <- function(bam, disable.shift, paired.end, outdir, thread) {
  rootname <- rnaSeqFanc::trim.bam(bam, trim_dir=T)
  cmd <- paste0("python3 /bar/cfan/software/atac/encode_pipeline/src/encode_task_bam2ta.py ",bam,
                " --out-dir ", outdir,
                " --nth ", thread)
  if (paired.end==T)
    cmd = paste0(cmd, " --paired-end ")
  if (disable.shift == T)
    cmd = paste0(cmd, " --disable-tn5-shift ")
  else rootname <- paste0(rootname, ".tn5")

  print(cmd)
  system(cmd)

  outfile <- paste0(outdir, "/", rootname, ".tagAlign.gz")
  if (!file.exists(outfile))
    stop(paste0("error: bam2tagAlign failed for bam file: ", bam))
  else return(outfile)
}



macs2.atac <- function (ta, genome, peak.pval, smooth.win, cap.num.peak=300000, outdir) {

  script <- "/bar/cfan/software/atac/encode_pipeline/src/encode_task_macs2_atac.py"

  interface <- tempfile()
  if (grepl("hg", genome))
    gensz <- "hs"
  else gensz <- "mm"
  cmd <- paste0("/bar/cfan/anaconda2/envs/jupyter/bin/python ", script, " ", ta,
                " --gensz ", gensz,
                " --pval-thresh ", peak.pval,
                " --cap-num-peak ", cap.num.peak,
                " --out-dir ", outdir,
                " --interface ", interface,
                " --smooth-win ", smooth.win)
  print(cmd)
  system(cmd)

  output <- read.table(interface, as.is = T)
  outfile <- output %>% filter(V1 == "narrowPeak") %>% pull(V2)
  if (!file.exists(outfile))
    stop(paste0("error: macs2 peak calling failed for ta file: ", ta))
  else return(outfile)
}

diffbind.fanc <- function(samples=NULL, dbo = NULL, rootname, outdir, ...) {
  # ...: additional commands to pass on to diffbind.
  if (is.null(dbo)) {
    if (!is.data.frame(samples))
      samples <- read.table(samples, as.is = T, header = T)
    else (stop("either samples (df or file location) or an already counted diffbind object needs to be offered"))
    dbo <- dba(sampleSheet = samples)
    system(paste0("mkdir -p ", outdir))

    png(file=paste0(outdir, "/", rootname, "_occupancy_hm.png"),
        height = 500, width = 500, units = "px")
    dba.plotHeatmap(dbo)
    dev.off()

    dbo <- dba.count(dbo, ...)
  }



  png(file=paste0(outdir, "/", rootname, "_affinity_hm.png"),
      height = 500, width = 500, units = "px")
  dba.plotHeatmap(dbo)
  dev.off()

  dbo <- dba.contrast(dbo, categories = DBA_FACTOR)
  dbo <- dba.analyze(dbo)

  png(file=paste0(outdir, "/", rootname, "_pca.png"),
      height = 500, width = 500, units = "px")
  dba.plotPCA(dbo, DBA_FACTOR,  label=DBA_FACTOR)
  dev.off()

  contrast <- mclapply(seq_along(dbo$contrasts), function(i) {
    options(scipen = 4)
    name1 <- dbo$contrasts[[i]]$name1
    name2 <- dbo$contrasts[[i]]$name2
    df <- dba.report(dbo, contrast = i) %>% as.data.frame()
    write.table(df, paste0(outdir, "/", rootname, "_",name1 , "_V_", name2, ".txt"),
                row.names = F, col.names = F, quote = F, sep = "\t")

    write.table(df[df[,6] > df[,7],], paste0(outdir, "/", rootname, "_",name1 , "_VV_", name2, ".txt"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(df[df[,7] > df[,6],], paste0(outdir, "/", rootname, "_",name2 , "_VV_", name1, ".txt"),
                row.names = F, col.names = F, quote = F, sep = "\t")

    return(df)
  }, mc.cores = 16)
  return(list(dba = dbo, contrast=contrast))
}

