
bam.count <- function(sample.info, bam.col = "bam", sample.col = "sample",
                      features.gr, bSingleEnd, sort = T,
                      return.mat = F, threads = 1) {
  # written with reference to sync_cage_to_dna_assays.R in cageFanc
  if (length(unique(sample.info[, sample.col])) != nrow(sample.info)) {
    stop("length(unique(sample.info[, sample.col])) != 1")
  }
  bpparam <- MulticoreParam(workers = threads)
  if (is.character(features.gr))
    features.gr <- utilsFanc::loci.2.gr(loci = features.gr)
  if (sort) {
    features.gr <- sort(features.gr)
  }
  se <- GenomicAlignments::summarizeOverlaps(features = features.gr, reads = sample.info[, bam.col],
                                             mode = "Union", inter.feature = T,
                                             singleEnd = bSingleEnd, fragments = F,
                                             ignore.strand = T, BPPARAM = bpparam)
  if (!identical(colnames(se), basename(sample.info[, bam.col]))) {
    stop("!identical(colnames(se), basename(sample.info[, bam.col]))")
  }
  colnames(se) <- sample.info[, sample.col]
  rownames(sample.info) <- sample.info[, sample.col]
  colData(se) <- DataFrame(sample.info)
  if (return.mat) {
    mat <- assay(se)
    rownames(mat) <- utilsFanc::gr.get.loci(rowRanges(se))
    return(mat)
  }
  return(se)
}


bam.promoter.pct <- function(bams, genomes, genome = NULL, mode = "ChIPseeker",
                             s = 42.001,  skip.subset = F, rm.subset = T,
                             bed.ext.1side = 1000, # only used in bed mode
                             scale.fac = 3,
                             out.file = NULL, 
                             threads = 1, samtools = SAMTOOLS) {
  modes.avail <- c("ChIPseeker", "liftover", "bed")
  # when using the bed mode, pass the bed file to the genomes argument 
  if (tolower(mode) == "liftover") {
    genomes.uniq <- unique(genomes)
    if (length(genomes.uniq) != 2) {
      stop("currently only supporting 2 genomes with the liftover mode")
    }
    liftover.dir <- paste0(genomes.uniq, collapse = "_")
    beds <- lapply(genomes.uniq, function(genome) {
      utilsFanc::import.bed.fanc(paste0("~/genomes/genes_shared/", liftover.dir, 
             "/shared_tss_", genome, ".bed"), return.gr = T) %>% return()
    })
    names(beds) <- genomes.uniq
  } 
  if (!is.null(genome)) {
    genomes <- rep(genome, length(bams))
    # back compatibility
  }
  if (length(bams) != length(genomes)) {
    if (length(genomes) == 1) {
      genomes <- rep(genomes, n = length(bams))
    } else {
      stop("length(bams) != length(genomes)")
    }
  }
  df <- utilsFanc::safelapply(seq_along(bams), function(n) {
    bam <- bams[n]
    print(paste0("processing file: ", bam))
    genome <- genomes[n]
    total.reads <- system(paste0(samtools, " view -c ", bam), intern = T) %>% 
      as.numeric()
    bam.s <- paste0(bam, "_", s, ".bam")
    cmd <- paste0(samtools, " view -hbo ", bam.s, " -s ", s, " ", bam)
    if (!skip.subset | !file.exists(bam.s)) {
      print(cmd); system(cmd)
    }
    gr <- suppressWarnings(GenomicAlignments::readGAlignments(bam.s) %>% GRanges())
    seqlengths(gr) <- NA
    strand(gr) <- "*"
    browser()
    if (tolower(mode) == "chipseeker") {
      gr <- utilsFanc::gr.fast.annot(gr, genome, use.klraps = F, anno.cols = "annotation")
      n.pro.1k <- sum(gr$annotation == "Promoter (<=1kb)")
    } else if (tolower(mode) == "liftover") {
      gr <- gr.annot.by.bed(gr= gr, bed = beds[[genome]], bed.ext.each.side = 1000, 
                            anno.col.name = "annotation", name.hits.as = "Pro1k")
      n.pro.1k <- sum(gr$annotation %in% "Pro1k")
    } else if (tolower(mode) == "bed") {
      bed <- genome
      gr <- gr.annot.by.bed(gr= gr, bed = bed, bed.ext.each.side = bed.ext.1side, 
                            anno.col.name = "annotation", name.hits.as = "Pro1k")
      n.pro.1k <- sum(gr$annotation %in% "Pro1k")
    } else {
      stop(paste0("mode ", mode, " is not in supported modes: ", 
                  paste0(modes.avail, collapse = ", ")))
    }
    # n.pro.2k <- sum(gr$annotation == "Promoter (1-2kb)")
    # n.pro.3k <- sum(gr$annotation == "Promoter (2-3kb)")
    # n.pro <- n.pro.1k + n.pro.2k + n.pro.3k
    # pct <- n.pro/length(gr)
    df <- data.frame(bam = bam, total.reads = total.reads, total.sampled = length(gr),
                     # n.pro = n.pro, pct.pro = pct,
                     n.pro.1k = n.pro.1k, pct.pro.1k = n.pro.1k/length(gr) # ,
                     # n.pro.2k = n.pro.2k, pct.pro.2k = n.pro.2k/length(gr),
                     # n.pro.3k = n.pro.3k, pct.pro.3k = n.pro.3k/length(gr)
    )
    if (rm.subset)
      system(paste0("rm -rf ", bam.s))
    return(df)
  }, threads = threads) %>% do.call(rbind, .)
  df <- norm.fac.cal(df, scale.fac = scale.fac)
  if (mode == "bed") {
    colnames(df) <- sub("pro.1k", paste0("bed.", utilsFanc::so.formatter(bed.ext.1side)), colnames(df))
  }
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    write.table(df, out.file, sep = "\t", row.names = F, col.names = T, quote = F)
  }
  return(df)
}

norm.fac.cal <- function(df, scale.fac = 3, slot = "pct.pro.1k") {
  # meant to be called by bam.promoter.pct
  # first estimate size factors from sequencing depth:
  # choice for scale.fac = 3: pct.pro.1k is around 0.3
  df$depth.fac <- df$total.reads/median(df$total.reads)
  df$scale.fac <- scale.fac
  df$final.fac <- df$depth.fac * df$scale.fac * df[[slot]]
  return(df)
}


gr.annot.by.bed <- function(gr, bed, bed.ext.each.side = 0, anno.col.name, name.hits.as = "hit") {
  if (is.character(bed)) {
    bed <- utilsFanc::import.bed.fanc(bed = bed, return.gr = T)
  }
  bed <- bed + bed.ext.each.side
  o <- findOverlaps(gr, bed, ignore.strand = T)
  mcols(gr)[, anno.col.name] <- NA
  mcols(gr)[, anno.col.name][sort(unique(queryHits(o)))] <- name.hits.as
  return(gr)
}

bam.filter.by.qname <- function(bam, read.names, out.bam) {
  f <- function(x) {
    res <- x$qname %in% read.names
    return(res)
  }
  dir.create(dirname(out.bam), showWarnings = F, recursive = T)
  Rsamtools::filterBam(file = bam, destination = out.bam, 
                       filter = FilterRules(f))
  return(out.bam)
}

bam.mean.coverage <- function(bam, bed) {
  # for each region in bed, calculate the mean coverage of reads in bam. Simply divide total bases 
  # by region lengths.
  # bed could also be just a gr. 
  if (is.character(bed)) {
    gr <- utilsFanc::import.bed.fanc(bed, return.gr = T)
  } else{
    gr <- bed
  }
  # if (is.null(gr$forth)) {
  #   gr$forth <- utilsFanc::gr.get.loci(gr)
  # }
  # utilsFanc::check.dups(gr$forth, "gr$forth")
  bam.gr <- GenomicAlignments::readGAlignments(
    file = bam, param = ScanBamParam(which = gr), with.which_label = T)
  if (any(njunc(bam.gr) != 0)) {
    stop("any(njunc(bam.gr) != 0) this means the width on reference is not accurate (such as splicing)")
  }
  
  cov.df <- data.frame(
    width = width(bam.gr), region = bam.gr@elementMetadata$which_label %>% as.character()) %>% 
    dplyr::group_by(region) %>% dplyr::summarise(sum = sum(width)) %>% 
    dplyr::ungroup() %>% as.data.frame()
  cov.df <- cov.df %>% utilsFanc::loci.2.df(loci.col.name = "region")
  cov.df <- cov.df %>% dplyr::mutate(mean = sum/(end - start + 1)) %>% 
    dplyr::mutate(chr = NULL, start = NULL, end = NULL)
  return(cov.df)
}

