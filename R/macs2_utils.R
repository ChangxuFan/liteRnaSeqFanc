
macs2.easy.fanc <- function(params, pseudo=F, zip.bedgraph.dir=T, force.zip = F) {
  cmd <- paste0("/opt/apps/python2/bin/macs2 ", params)
  print(cmd)
  if (pseudo == F)
    system(cmd)
  
  if (zip.bedgraph.dir == T) {
    patterns <- c("^.+outdir +([^ ]+).*$")
    dir <- lapply(patterns, function(x) {
      if (grepl(x, params))
        return(sub(x, "\\1", params))
      else return(NULL)
    }) %>% unlist()
    if (length(dir) != 1)
      stop(paste0("error in finding outdir: ", length(dir), " patterns were matched"))
    print(paste0("the dir grepped out is: ", dir))
    bedgraphs <- Sys.glob(paste0(dir, "/*.bdg"))
    lapply(bedgraphs, function(x) {
      if (file.exists(paste0(x, ".gz")) && force.zip != T)
        print(paste0("the zipped form of ", x, " already exists. use force.zip == T to overwrite"))
      else {
        cmd <- paste0("/bar/cfan/scripts/bed_browser_v2.sh ", x)
        print(cmd)
        if (pseudo == F) 
          system(cmd)
        return(NULL)
      }
    })
  }

  return(NULL)
  
}

zip.files.by.dir <- function(dir, patterns.vector, zip = F, as.bw=F, genome=NULL,thread) {
  files.to.zip <- lapply(patterns.vector, function(x) {
    Sys.glob(paste0(dir, "/", x))
  }) %>% unlist() 
  
  print("the files to zip are: ")
  print(files.to.zip)
  
  if (zip == T && length(files.to.zip) > 0) {
    print("zipping files ...")
    trash <- mclapply(files.to.zip, function(x) {
      if (as.bw == T) {
        if (!is.null(genome)) {
          cmd <- paste0("/opt/apps/kentUCSC/334/bedGraphToBigWig ", x, " ",
                        "/bar/cfan/genomes/", genome, "/", genome, ".chrom.sizes ", x, ".bw")
        } else
          stop("when as.bw is TRUE, genome must be supplied")
      } else 
        cmd <- paste0("/bar/cfan/scripts/bed_browser_v2.sh ", x)
      print(cmd)
      system(cmd)
      return(NULL)
    }, mc.cores = min(length(files.to.zip) + 1, thread))
  }
  return(NULL)
}



# untested:::
macs2.atac.callpeak <- function(infile, root.name, outdir, format, genome, macs2.path = "/opt/apps/python2/bin/macs2",
                           p.cutoff=NULL, q.cutoff=NULL, shift=NULL, ext=NULL, spmr=F, subsummit=F,
                           bdgcmp = T, run=T, zip = T, write.log = T, bdgcmp.method = NULL, thread = 5,
                           add.param = "") {
  cmd <- paste0(macs2.path, " callpeak ", 
                " --nomodel -B --keep-dup all ",
                " -t ", infile,
                " -n ", root.name,
                " --outdir ", outdir,
                " -f ", format,
                " -g ", genome)
  cutoff <- ""
  if (!is.null(p.cutoff))
    cutoff <- paste0(" -p ", p.cutoff)
  if (!is.null(q.cutoff))
    cutoff <- paste0(" -q ", q.cutoff) # this means that if p and q are both given, q overrides p.
  if (cutoff != "")
    cmd <- paste0(cmd, cutoff)
  
  if (!is.null(shift))
    cmd <- paste0(cmd, " --shift ", shift)
  if (!is.null(ext))
    cmd <- paste0(cmd, " --extsize ", ext)
  if (spmr == T)
    cmd <- paste0(cmd, " --SPMR")
  if (subsummit == T)
    cmd <- paste0(cmd, " --call-summits")
  
  logfile <- NULL
  if (write.log == T) {
    logfile <- paste0(outdir, "/", root.name, "_callpeak.log")
  }
  cmd <- paste0(cmd, " ", add.param)

  trash <- utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = logfile, intern = F, run = run)
  
  trash <- zip.files.by.dir(dir = outdir, 
                            patterns.vector = paste0(root.name,c("_control_lambda.bdg", "_treat_pileup.bdg",
                                                                 "_summits.bed", "_peaks.narrowPeak")),
                            zip = zip, thread = thread)
  if (bdgcmp == T) {
    if (!is.null(bdgcmp.method))
      trash <- macs2.atac.cmp(t.bdg = paste0(outdir, "/", root.name, "_treat_pileup.bdg"),
                              c.bdg = paste0(outdir, "/", root.name, "_control_lambda.bdg"), 
                              outdir = outdir, root.name = root.name, method = bdgcmp.method, run = run, zip = zip, 
                              macs2.path = macs2.path, write.log = write.log)
    else 
      stop("when cmp is called, method parameter has to be supplied")
  }

  
  return(paste0(outdir, "/", root.name))
}

macs2.bdg.count.total.fragments <- function(bdg, frag.length=NULL, run = T) {
  if (is.null(frag.length)) 
    frag.length <- 1
  cmd <- paste0("cat ", bdg, " | awk -F \"\\t\" 'BEGIN{sum=0} {sum = sum + ($3-$2)*$4} END {printf sum}'")
  if (grepl(".gz$", bdg)) {
    cmd <- paste0("z", cmd)
  }
  print(Sys.time())
  cat(cmd)
  cat("\n")
  if (run == T) {
    n.frag <- system(cmd, intern = T) %>% as.numeric()
    return(n.frag/frag.length)
  }
  else 
    return(NULL)
}



macs2.atac.cmp <- function(t.bdg, c.bdg, outdir, root.name, method, run = T, zip = T,
                           macs2.path = "/opt/apps/python2/bin/macs2", write.log=T, thread = 2) {
  if (!file.exists(t.bdg))
    stop(paste0(t.bdg, " does not exist"))
  if (!file.exists(c.bdg))
    stop(paste0(c.bdg, " does not exist"))
  
  cmd <- paste0(macs2.path, " bdgcmp ", " -t ", t.bdg, " -c ", c.bdg, " --outdir ", outdir, 
                " --o-prefix ", root.name, " -m ", method)
  
  logfile <- NULL
  if (write.log == T) {
    logfile <- paste0(outdir, "/", root.name, "_bdgcpm_", method,".log")
  }
  trash <- utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = logfile, intern = F, run = run)
  
  trash <- zip.files.by.dir(dir = outdir, patterns.vector = paste0(root.name, "_",method, ".bdg"),
                            zip = zip, thread = thread)
  
  return(paste0(outdir, "/", root.name))
  
}


# macs2.bdg.count.total.fragments("sth/macs2/ly49apos2_no_SPMR/ly49apos2_no_SPMR_treat_pileup.bdg", T)

macs2.atac.diff <- function(t1.bdg, t2.bdg, c1.bdg=NULL, c2.bdg = NULL, llr.cutoff, min.length, max.gap, depth.1=NULL,
                            depth.2=NULL, frag.length.1=NULL, frag.length.2=NULL,outdir, 
                            root.name, macs2.path = "/opt/apps/python2/bin/macs2", run = T,
                            zip=T, write.log = T, thread = 4) {
  # note thread is used for zipping only
  if (is.null(depth.1)) {
    depth.1 <- macs2.bdg.count.total.fragments(bdg = t1.bdg, frag.length = frag.length.1, run = T )
  }
  if (is.null(depth.2)) {
    depth.2 <- macs2.bdg.count.total.fragments(bdg = t2.bdg, frag.length = frag.length.2, run = T )
  }
  
  if (is.null(c1.bdg)) {
    c1.bdg <- sub("treat_pileup","control_lambda",t1.bdg)
  }
  
  if (is.null(c2.bdg)) {
    c2.bdg <- sub("treat_pileup","control_lambda",t2.bdg)
  }
  
  for (i in c(t1.bdg, t2.bdg, c1.bdg, c2.bdg)) {
    if (!file.exists(i))
      stop(paste0("miao~ file ", i, " does not exist~"))
  }
  
  cmd <- paste0(macs2.path, " bdgdiff ", 
                " --t1 ", t1.bdg, " --t2 ", t2.bdg, " --c1 ",c1.bdg, " --c2 ", c2.bdg, 
                " -C ", llr.cutoff, " -l ", min.length, " -g ", max.gap,
                " --d1 ", depth.1, " --d2 ", depth.2,
                " --outdir ", outdir, " --o-prefix ", root.name)
  
  logfile <- NULL
  if (write.log == T) {
    logfile <- paste0(outdir, "/", root.name, "_bdgdiff.log")
  }
  trash <- utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = logfile, intern = F, run = run)
  
  trash <- zip.files.by.dir(dir = outdir, patterns.vector = paste0(root.name, "*bed"),
                            zip = zip, thread = thread)
  
  return(paste0(outdir, "/", root.name))
}

# t <- system("cat sth/macs2/ly49apos2_no_SPMR/ly49apos2_no_SPMR_treat_pileup.bdg | awk -F \"\\t\" 'BEGIN{sum=0} {sum = sum + ($3-$2)*$4} END {printf sum}' & ", intern = T)

