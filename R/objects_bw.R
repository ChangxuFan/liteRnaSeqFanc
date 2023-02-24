bw.zero <- function(in.bw, out.bed, bgzip = BGZIP, tabix = TABIX) {
  in.gr <- rtracklayer::import(in.bw)
  # chrom.sizes <- seqinfo(in.gr) %>% as.data.frame() %>% .[, 1, drop = F] %>% 
  #   mutate(., chr = rownames(.), start = 1, end = seqlengths) %>% 
  #   makeGRangesFromDataFrame()
  zero.gr <- in.gr[in.gr$score == 0]
  zero.gr$score <- width(zero.gr)
  colnames(mcols(zero.gr)) <- "size"
  dir.create(dirname(out.bed), showWarnings = F, recursive = T)
  rtracklayer::export.bed(zero.gr, out.bed)
  system(paste0(bgzip, " -c -f ", out.bed, " > ", out.bed, ".gz"))
  system(paste0(tabix, " -f -p bed ", out.bed, ".gz"))
  return(out.bed)
}

bdg.average <- function(in.bdgs, regions.gr = NULL, out.bdg,
                        tmp.dir = NULL, bedtools = BEDTOOLS,
                        scale.to = "mean") {
  # scale.to: NULL: don't do scaling; "mean": normalize to the mean of all in.bdgs.
  if (is.null(tmp.dir)) {
    tmp.dir <- tempdir()
  }
  if (is.character(regions.gr)) {
    regions.gr <- utilsFanc::loci.2.gr(regions.gr)
  }
  dir.create(tmp.dir, showWarnings = F, recursive = T)
  union.bdg <- paste0(tmp.dir, "/", basename(out.bdg), "_union.bdg")
  cmd <- paste0(bedtools, " unionbedg -i ", paste0(in.bdgs, collapse = " "),
                " > ", union.bdg)
  print(cmd); system(cmd)
  gr <- rtracklayer::import.bedGraph(union.bdg)
  if (!is.null(regions.gr)) {
    gr <- subsetByOverlaps(gr, regions.gr)
  }
  mat <- mcols(gr) %>% as.data.frame() %>% as.matrix()
  if (!is.null(scale.to)) {
    if (scale.to == "mean") {
      scale.to <- mean(colSums(mat))
    }
    mat <- mat %*% diag(scale.to/colSums(mat))
  }
  mcols(gr) <- NULL
  gr$score <- rowMeans(mat)
  dir.create(dirname(out.bdg), showWarnings = F, recursive = T)
  rtracklayer::export.bedGraph(object = gr, out.bdg)
  cmd <- paste0("~/scripts/bb.sh ", out.bdg)
  print(cmd); system(cmd)
  return(out.bdg)
}

bdg.fill.zero <- function(bdg, out.dir = NULL, add = 0,
                          chrom.sizes = NULL, genome = NULL, seq.lengths = NULL) {
  suffix <- sub(".gz$", "", bdg) %>% tools::file_ext()
  root <- sub(paste0("\\.", suffix, ".*$"), "" , basename(bdg))
  if (is.null(out.dir)) {
    out.dir <- dirname(bdg)
  }
  dir.create(out.dir, showWarnings = F, recursive = T)
  if (tolower(suffix) %in% c("bdg", "bedgraph")) {
    gr <- rtracklayer::import.bedGraph(bdg)
    bdg.flag <- T
    if (is.null(seq.lengths)) {
      if (is.null(chrom.sizes)) {
        if (is.null(genome)) {
          seq.lengths <- utilsFanc::gr.get.seq.max(gr)
        } else {
          chrom.sizes <- paste0("~/genomes/", genome, "/", genome, ".chrom.sizes")
          seq.lengths.df <- read.table(chrom.sizes, sep = "\t", header = F)
          seq.lengths <- seq.lengths.df$V2
          names(seq.lengths) <- seq.lengths.df$V1
        }
      } else {
        seq.lengths.df <- read.table(chrom.sizes, sep = "\t", header = F)
        seq.lengths <- seq.lengths.df$V2
        names(seq.lengths) <- seq.lengths.df$V1
      }
    }
    avail.chroms <- seqlengths(gr) %>% names()
    utilsFanc::check.intersect(avail.chroms, x.name = "avail.chroms", 
                               names(seq.lengths), y.name = "names(seq.lengths)")
    seq.lengths <- seq.lengths[avail.chroms]
    seqlengths(gr) <- seq.lengths
  } else if (tolower(suffix) %in% c("bw", "bigwig")) {
    gr <- rtracklayer::import.bw(bdg)
    bdg.flag <- F
    seq.lengths <- seqlengths(gr)
  }
  
  cov <- coverage(gr, weight = "score") %>% as("GRanges")
  if (add != 0) {
    cov$score[cov$score == 0] <- add
  }
  
  out.bdg <- paste0(out.dir, "/", root, "_fill", add, ".", suffix)
  if (!bdg.flag) {
    seqlengths(cov) <- seq.lengths
    rtracklayer::export(cov, out.bdg)
  } else {
    rtracklayer::export.bedGraph(cov, out.bdg)
    system(paste0("~/scripts/bed_browser_ez.sh ", out.bdg))
    out.bdg <- paste0(out.bdg, ".gz")
  }
  cat(utilsFanc::bash2ftp(out.bdg))
  cat("\n")
  return(out.bdg)
}

bdg.ext.smooth <- function(in.bdg, out.dir = NULL, ext) {
  if (is.null(out.dir)) {
    out.dir <- dirname(in.bdg)
  }
  dir.create(out.dir, recursive = T, showWarnings = F)
  out.bdg <- paste0(out.dir, "/", 
                    utilsFanc::insert.name.before.ext(name = basename(in.bdg), 
                                                      insert = paste0("ext_", ext), delim = "_"))
  gr <- rtracklayer::import(in.bdg) 
  cov <- utilsFanc::gr.expand.smooth.bp(gr, ext = ext)
  rtracklayer::export(cov, sub(".gz$", "", out.bdg))
  if (grepl("(bedgraph|bdg)(.gz)*$", tolower(out.bdg))) {
    system(paste0("~/scripts/bed_browser_ez.sh ", sub(".gz$", "", out.bdg)))
  }
  if (file.exists(paste0(out.bdg, ".gz"))) {
    cat(utilsFanc::bash2ftp(paste0(out.bdg, ".gz")))
  }
  cat(utilsFanc::bash2ftp(out.bdg))
  cat("\n")
  return(out.bdg)
}


bw.reNorm <- function(in.vec, out.vec = NULL, factor, threads = 1) {
  if (length(factor) == 1 && grepl("[a-zA-Z]", factor)) {
    factor <- readLines(factor)
  }
  factor <- as.numeric(factor)
  if (any(is.na(factor))) {
    stop("any(is.na(factor))")
  }
  
  if (length(factor) != length(in.vec)) {
    stop("length(factor) != length(in.vec)")
  }
  
  if (is.null(out.vec)) {
    out.vec <- tools::file_path_sans_ext(in.vec) %>%
      paste0("_reNorm_", round(factor, digits = 2), ".bw")
  }
  utilsFanc::safelapply(seq_along(in.vec), function(i) {
    in.bw <- in.vec[i]
    out.bw <- out.vec[i]
    factor <- factor[i]
    in.gr <- rtracklayer::import.bw(in.bw)
    out.gr <- in.gr
    out.gr$score <- out.gr$score / factor
    dir.create(dirname(out.bw), showWarnings = F, recursive = T)
    rtracklayer::export.bw(out.gr, out.bw)
    return()
  }, threads = threads)
  
}

bw.count <- function(sample.info, bw.col = "bw", sample.col = "sample",
                     features.gr, sort = T,
                     threads = 1) {
  if (length(unique(sample.info[, sample.col])) != nrow(sample.info)) {
    stop("length(unique(sample.info[, sample.col])) != 1")
  }
  
  if (is.character(features.gr))
    features.gr <- utilsFanc::loci.2.gr(loci = features.gr)
  if (sort) {
    features.gr <- sort(features.gr)
  }
  features.gr$id <- utilsFanc::gr.get.loci(features.gr)
  strand(features.gr) <- "*"
  mat <- utilsFanc::safelapply(1:nrow(sample.info), function(n) {
    bw <- sample.info[n, bw.col]
    gr.bw <- rtracklayer::import(bw)
    j <- plyranges::join_overlap_left(features.gr, gr.bw)
    df <- mcols(j) %>% as.data.frame() %>% 
      group_by(id) %>% summarise(score = sum(score)) %>% 
      as.data.frame()
    df <- df[gtools::mixedorder(df$id),]
    if (!identical(df$id, features.gr$id)) {
      stop("!identical(df$id, features.gr$id)")
    }
    return(df$score)
  }, threads = threads) %>% do.call(cbind, .)
  rownames(mat) <- features.gr$id
  colnames(mat) <- sample.info[, sample.col]
  mat[is.na(mat)] <- 0
  return(mat)
}