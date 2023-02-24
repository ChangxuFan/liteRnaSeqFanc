plotAnnoPie.m <- function(anno.list, samples = NULL,
                          width = 10, height = 10,
                          sub.width = 800, sub.height = 500,
                          plot.out,
                          ...) {
  if (is.null(names(anno.list))) {
    names(anno.list) <- samples
  }
  if (is.null(names(anno.list)))
    stop("anno.list must be named. Or samples need to be provided")
  pl <- mapply(function(anno, name) {
    out.file <- tempfile()
    png(filename = out.file, width = sub.width, height = sub.height, res = 100, 
        ...)
    try((ChIPseeker::plotAnnoPie(anno, main = "\n\n\n" %>% paste0(name))))
    dev.off()
    p <- png::readPNG(out.file) %>% grid::rasterGrob()
    return(p)
  }, anno.list, names(anno.list), SIMPLIFY = F)
  p <- gridExtra::grid.arrange(grobs = pl)
  system(paste0("mkdir -p ", dirname(plot.out)))
  ggsave(plot.out, p, height = height, width = width)
  return(p)
}

wrap.png <- function(pngs, plot.out, height, width, ...) {
  stop("this function has been moved to utilsFanc::png.wrap()")
  # rl = lapply(pngs, png::readPNG)
  # gl = lapply(rl, grid::rasterGrob)
  # p <- gridExtra::grid.arrange(grobs=gl)
  # system(paste0("mkdir -p ", dirname(plot.out)))
  # ggsave(plot.out, p, width = width, height = height, ...)
  # return(p)
}
reads.distro <- function(bams, samples = NULL, genome, subsample.to = 200000, 
                         pre.subsample.frac = NULL, force = F,
                         plot.out, height = 10, width = 10,
                         threads = 1, 
                         samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  if (is.null(samples))
    samples <- names(bams)
  if (is.null(samples)) {
    stop("samples must be supplied or bams must be named") 
  }
  if (length(bams) != length(samples)) {
    stop("length(bams) != length(samples)")
  }
  bam.not.found <- bams %>% .[!file.exists(bams)]
  if (length(bam.not.found) > 0) {
    stop(paste0("bam files not found: ", paste0(bam.not.found, collapse = "\n")))
  }
  if (genome == "mm10") {
    TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    annoDb = "org.Mm.eg.db"
  } else if (genome == "hg38") {
    annoDb <- "org.Hs.eg.db"
    TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    
  } else {
    stop("only mm10 and hg38 genomes has been developed.")
  }
  tempt.dir <- paste0(dirname(plot.out), "/tempt/")
  dir.create(tempt.dir, showWarnings = F, recursive = T)
  res <- utilsFanc::safelapply(seq_along(bams), function(i) {
    bam <- bams[[i]]
    sample <- samples[[i]]
    if (!is.null(pre.subsample.frac)) {
      bam.ori <- bam
      bam <- bam.ori %>% utilsFanc::insert.name.before.ext(paste0("sub_", pre.subsample.frac), delim = "_")
      dirname(bam) <- dirname(bam) %>% paste0("/sub_", pre.subsample.frac, "/")
      dir.create(dirname(bam), showWarnings = F, recursive = T)
      if (file.exists(bam) && force == F) {
        print(paste0("bam file \n", bam, "\nalready exists, use force = T to regenerate it"))
      } else {
        cmd <- paste0(samtools, " view -s ", pre.subsample.frac, " -hbo ", bam, " ", bam.ori)
        print(cmd); system(cmd)
      }
    }
    
    gal.se <- readGAlignments(bam)
    subsample.to <- min(subsample.to, length(gal.se))
    gal.se <- gal.se %>% .[sample(1:length(.), size = subsample.to, replace = F)]
    
    gr.se <- gal.se %>% as("GRanges")
    anno <- ChIPseeker::annotatePeak(peak = gr.se, TxDb = TxDb, annoDb = annoDb)
    plot.sub <- paste0(tempt.dir, "/", sample, ".png")
    cageFanc::plotAnnoPie.fanc(anno = anno, out.file = plot.sub, 
                               main = paste0("\n\n",sample, "\nreads:", subsample.to %>% utilsFanc::so.formatter(),
                                             " pre.subsample: ", pre.subsample.frac))
    res <- list(anno = anno, plot.sub = plot.sub, gal.se = gal.se, subsample.to = subsample.to, 
                pre.subsample.frac = pre.subsample.frac)
    return(res)
  }, threads = threads)
  save.rds <- tools::file_path_sans_ext(plot.out) %>% paste0(".Rds")
  system(paste0("mkdir -p ", dirname(save.rds)))
  saveRDS(res, save.rds)
  p <- wrap.png(pngs = sapply(res,function(x) return(x$plot.sub)),
                plot.out = plot.out, height = height, width = width)
  return(p)
}

bam.2.gr.se <- function(bam, subsample.to = NULL) {
  gal.se <- readGAlignments(bam)
  if (!is.null(subsample.to)) {
    subsample.to <- min(subsample.to, length(gal.se))
    gal.se <- gal.se %>% .[sample(1:length(.), size = subsample.to, replace = F)]
  }
  gr.se <- gal.se %>% as("GRanges")
  return(gr.se)
}

import.narrowPeak.fanc <- function(file, return.gr = F, broadPeak = F) {
  df <- read.table(file, sep = "\t", header = F)
  if (broadPeak == T)
    n.col <- 9
  else
    n.col <- 10
  if (ncol(df) != n.col)
    stop("ncol(df) != n.col")
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "signal", "p", "q", "peak") %>% 
    .[1:n.col]
  if (return.gr == T)
    return(makeGRangesFromDataFrame(df = df, keep.extra.columns = T))
  else
    return(df)
}

motif.in.peaks <- function(peaksets, samples, motif.gr,
                           score.col = 9,
                           top.n.vec = c(10000, 30000, 50000),
                           out.file,
                           threads = 1) {
  # narrowPeak format is expected.
  if (length(peaksets) != length(samples))
    stop("length(peaksets) != length(samples)")
  if (is.data.frame(motif.gr))
    motif.gr <- makeGRangesFromDataFrame(df = motif.gr, keep.extra.columns = T)
  peaksets <- lapply(peaksets, function(peakset) {
    if (is.character(peakset))
      peakset <- read.table(peakset, header = F, sep = "\t")
    peakset <- peakset[rev(order(peakset[, score.col])),]
    if (!("chr" %in% colnames(peakset))) {
      colnames(peakset)[1:6] <- c("chr", "start", "end", "peak_id", "score", "strand")
    }
    peakset <- makeGRangesFromDataFrame(peakset, keep.extra.columns = T)
    return(peakset)
  })
  stat.df <- utilsFanc::safelapply(peaksets, function(peakset) {
    fracs <- sapply(top.n.vec, function(top.n) {
      peaks <- peakset %>% .[1:min(length(.), top.n)]
      n.peaks <- length(peaks)
      o <- findOverlaps(query = peaks, subject = motif.gr)
      n.w.motif <- queryHits(o) %>% unique() %>% length()
      frac <- round(n.w.motif/n.peaks, digits = 3)
      return(frac)
    })
    return(fracs)
  }, threads = threads) %>% Reduce(cbind, .) %>% as.data.frame()
  
  colnames(stat.df) <- samples
  stat.df <- cbind(data.frame(top.n = top.n.vec), stat.df)
  rownames(stat.df) <- NULL
  dir.create(dirname(out.file), showWarnings = F, recursive = T)
  write.table(stat.df, out.file, quote = F, col.names = T, row.names = F, sep = "\t")
  n.peaks <- sapply(peaksets, length)
  names(n.peaks) <- samples
  write(paste0("## ", samples, ": ", n.peaks), out.file, sep = "\n", append = T)
  res <- list(stat.df = stat.df, n.peaks = n.peaks)
  return(res)
}

frip.fanc <- function(sample.df, score.col = 9, subsample.to = NULL,
                      broadPeak = F,
                      top.n.vec = c(10000, 20000, 30000, 40000),
                      out.file) {
  # sample.df: sample, bam, peakset (narrowpeak format)
  if (is.character(sample.df))
    sample.df <- read.table(sample.df, sep = "\t", header = T)
  if (any(! c("sample", "bam", "peakset") %in% colnames(sample.df))) {
    stop("colnames of sample.df doesn't conform")
  }
  top.n.df <- data.frame(start = c(0, top.n.vec[-length(top.n.vec)]) + 1,
                         end = top.n.vec) %>% 
    mutate(interval = paste0(start, "-", end))
  n.peaks.df <- data.frame(interval = "n.peaks")
  frip.df <- lapply(1:nrow(sample.df), function(i) {
    x <- sample.df[i, ]
    reads <- bam.2.gr.se(bam = x$bam, subsample.to = subsample.to)
    n.reads <- length(reads)
    peakset <- import.narrowPeak.fanc(file = x$peakset, return.gr = T, broadPeak = broadPeak)
    n.peaks <- length(peakset)
    p <- parent.frame(2)
    p$n.peaks.df[, x$sample] <- n.peaks
    frips <- top.n.df %>% split(., f = 1:nrow(.)) %>% 
      sapply(function(y) {
        if (y$start > n.peaks)
          return(0)
        if (y$end > n.peaks) {
          int.length <- y$end - y$start
          n.peaks.left <- n.peaks - y$start
          if (n.peaks.left > 0.3 * int.length)
            scale.factor <- int.length/n.peaks.left
          else
            return(0)
        } else {
          scale.factor <- 1
        }
        peaks <- peakset[y$start:min(y$end, n.peaks)]
        o <- findOverlaps(query = reads, subject = peaks)
        n.in.peak <- queryHits(o) %>% unique() %>% length()
        frip <- scale.factor * n.in.peak/n.reads
        frip <- round(frip, digits = 3)
        return(frip)
      })
    return(frips)
  }) %>% Reduce(cbind, .) %>% as.data.frame()
  colnames(frip.df) <- sample.df$sample
  frip.df <- cbind(data.frame(interval = top.n.df$interval), frip.df)
  col.sum <- colSums(frip.df[, 2:ncol(frip.df)]) %>% c("sum", .)
  frip.df <- rbind(frip.df, col.sum)
  frip.df <- rbind(frip.df, n.peaks.df)
  dir.create(dirname(out.file), showWarnings = F,recursive = T)
  write.table(frip.df, file = out.file, quote = F, sep = "\t", col.names = T, row.names = F)
  return(frip.df)
}

macs2.cutrun <- function(sample.df, prefix = NULL,
                         sepe = NULL, p_q = NULL, thresh = NULL, outdir = NULL, genome = NULL,
                         no.model = F,
                         threads.each = 1, threads.master = NULL, run = T,
                         macs2 = "/opt/apps/python2/bin/macs2",
                         bedtools = "/opt/apps/python2/bin/bedtools") {
  # sample.df format: bam, control. no control: give NA. optional: sepe, p_q, thresh, outdir, genome
  if (is.character(sample.df))
    sample.df <- read.table(sample.df, header = T, sep = "\t")
  args <- c("sepe", "p_q", "thresh", "genome")
  for ( arg in args) {
    if (! arg %in% colnames(sample.df)) {
      arg.get <- get(arg, envir = environment(), inherits = F)
      if (is.null(arg.get)) {
        stop(paste0(arg, " needs to be specified"))
      } else {
        sample.df[, arg] <- arg.get
      }
    }
  }
  if (is.null(threads.master))
    threads.master <- nrow(sample.df)

  sample.df$bg.flag <- ""
  sample.df$bg.flag[! is.na(sample.df$control)] <- "_bg"
  sample.df <- sample.df %>% 
    dplyr::mutate(suffix = paste0(bg.flag, "_", p_q, thresh, "_", sepe))
  if (is.null(sample.df$outdir)) {
    sample.df <- sample.df %>% 
      dplyr::mutate(outdir = paste0(dirname(dirname(bam)), "/macs2", suffix))
  }
  sample.df <- sample.df %>% 
    dplyr::mutate(root.name = paste0(outdir, "/", tools::file_path_sans_ext(basename(bam)), suffix)) 
  
  if (!is.null(prefix)) {
    sample.df <- dplyr::mutate(sample.df, bam = paste0(prefix, "/bam/", bam),
                               control = paste0(prefix, "/bam/control"),
                               outdir = paste0(prefix, "/", outdir))
  }
  sample.df <- sample.df %>% filter(is.na(control) | bam != control)
  if (nrow(sample.df) < 1) {
    stop("nrow(sample.df) < 1")
  }
  lapply(c("bam", "control"), function(x) {
    file.to.test <- sample.df[, x] %>% .[!is.na(.)] %>% unique()
    if (length(file.to.test) > 0) {
      file.not.found <- file.to.test %>% .[!file.exists(.)]
      if (length(file.not.found) > 0)
        stop(paste0("files not found: \n", paste0(file.not.found, collapse = "\n")))
    }
  })
  out.df <- sample.df %>% split(., f = 1:nrow(.)) %>% 
    utilsFanc::safelapply(function(x) {
      cmd <- paste0(macs2, " callpeak -B -t ", x$bam, 
                    " -", x$p_q, " ", thresh)
      if (!is.na(x$control)) {
        cmd <- paste0(cmd, " -c ", x$control)
      } 
      
      if (no.model == T) {
        cmd <- paste0(cmd, " --nomodel")
      }
      
      if (x$sepe == "se") {
        cmd <- paste0(cmd, " -f BAM")
      } else {
        cmd <- paste0(cmd, " -f BAMPE")
      }
      cmd <- paste0(cmd, " -n ", x$root.name, " --keep-dup all")
      print(cmd)
      if (run == T) {
        dir.create(dirname(x$outdir), showWarnings = F, recursive = T)
        system(cmd)
      }
      x$np <- x$root.name %>% paste0("_peaks.narrowPeak")
      
      if (run == T && !file.exists(x$np)) {
        stop(paste0(x$root.name, " failed to generate"))
      }
      x$bl <- paste0("~/genomes/", genome, "/blacklist/", genome, ".blacklist.bed")
      if (!file.exists(x$bl)) {
        stop(paste0("blacklist file: ", x$bl, " not found"))
      }
      x$np_bl <- paste0(x$root.name, "_rmbl.narrowPeak")
      cmd <- paste0(bedtools, " intersect -v -a ", x$np, " -b ", x$bl, " > ", x$np_bl)
      print(cmd)
      if (run == T) {
        system(cmd)
        system(paste0("~/scripts/bb.sh ", x$np_bl))
      }
      return(x)
    }) %>% Reduce(rbind, .)
  return(out.df)
}


