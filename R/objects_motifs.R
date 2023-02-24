motif.enrich.vierstra <- function(fg, np.q.th = 5, n.fg, seed = NULL, n.bg.each, 
                                  remove.promoter = F, genome, tssRegion = c(-2000, 1000),
                                  np.buffer = 250,
                                  scan.beds, TFs.use = NULL, scan.score.cutoff,
                                  out.dir, root.name = NULL, threads = 10,
                                  tabix = TABIX) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  if (is.character(fg)) {
    fg <- rtracklayer::import(fg)
  }
  if (grepl("rds|Rds", scan.beds[1])) {
    scan.beds <- readRDS(scan.beds)
  }
  if (!is.null(TFs.use)) {
    scan.beds <- scan.beds[TFs.use]
  }
  if (is.null(names(scan.beds))) {
    stop("is.null(names(scan.beds))")
  }
  if (ncol(mcols(fg)) == 6) {
    print("dealing with narrowPeak format. Recentering to summit +- np.buffer bp")
    fg <- fg[fg$qValue > np.q.th]
    set.seed(seed = seed)
    fg <- fg[sort(sample(1:length(fg), size = n.fg, replace = F))]
    
    fg$start.old <- start(fg)
    start(fg) <- fg$peak +  fg$start.old - np.buffer
    end(fg) <- fg$peak +  fg$start.old + np.buffer
    
    if (remove.promoter) {
      if (genome == "mm10") {
        genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
        annoDb <- "org.Mm.eg.db"
        TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
      } else {
        stop("only mm10 supported right now")
      }
      anno <- ChIPseeker::annotatePeak(fg, tssRegion = tssRegion,
                                       TxDb = TxDb, annoDb = annoDb)@anno
      fg <- anno[!grepl("Promoter", anno$annotation)]
    }
    
    mcols(fg) <- NULL
  }
  
  bg <- motif.get.bg(fg.gr = fg, n.bg.each = n.bg.each)
  fg$.type <- "fg"
  bg$.type <- "bg"
  r <- c(fg, bg)
  r.out <- paste0(out.dir, "/", root.name, "_regions.bed")
  utilsFanc::write.zip.fanc(r, r.out, bed.shift = T)
  regionPositions <- utilsFanc::safelapply(seq_along(scan.beds), function(i) {
    bed <- scan.beds[i]
    print(paste0("process # ", i, " bed file"))
    o <- tempfile()
    cmd <- paste0(tabix, " -B ", bed, " ", r.out, " > ", o)
    system(cmd)
    gr <- utilsFanc::import.bed.fanc(o) %>% .[, 1:6] %>% 
      dplyr::filter(fifth >= scan.score.cutoff) %>% 
     makeGRangesFromDataFrame(keep.extra.columns = F)
    return(gr)
  }, threads = threads) %>% GRangesList()
  names(regionPositions) <- names(scan.beds)
  allPositions <- unlist(regionPositions)
  peakSet <- r
  overlapRegions <- findOverlaps(peakSet, allPositions, ignore.strand = TRUE)
  
  regionMat <- Matrix::sparseMatrix(i = queryHits(overlapRegions), 
                                    j = match(names(allPositions), 
                                              names(regionPositions))[subjectHits(overlapRegions)], 
                                    x = rep(TRUE, length(overlapRegions)), 
                                    dims = c(length(peakSet), length(regionPositions)))
  colnames(regionMat) <- names(regionPositions)
  
  regionMat <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(matches = regionMat), 
                                                          rowRanges = peakSet)
  saveRDS(regionMat, paste0(out.dir, "/", root.name, "_mat.Rds"))
  
  enrich <- motif.compute.enrich(matches.se = regionMat)
  write.table(enrich, paste0(out.dir, "/", root.name, "_enrich.tsv"), sep = "\t", 
              quote = F, row.names = F, col.names = T)
  return(enrich)
}

motif.get.bg <- function(fg.gr, n.bg.each) {
  if (length(unique(width(fg.gr)) ) != 1) {
    stop("length(unique(width(fg)) ) != 1")
  }
  w <- width(fg.gr)[1]
  buffer <- floor(w/2)
  fg.range <- fg.gr %>% as.data.frame() %>% 
    group_by(seqnames) %>% 
    dplyr::summarise(start = min(start), end = max(end),
                     n.fg.chr = n()) %>% 
    ungroup()
  bg <- fg.range %>% split(., f= 1:nrow(.))  %>% 
    lapply(function(x) {
      n.bg.chr <- x$n.fg.chr * n.bg.each
      centers <- sample(x = ceiling(x$start+w/2):floor(x$end-w/2),
                        size = n.bg.chr, replace = F)
      gr <- data.frame(chr = x$seqnames, start = centers - buffer, end = centers + buffer) %>% 
        makeGRangesFromDataFrame(keep.extra.columns = T)
      return(gr)
    }) %>% Reduce(c, .)
  return(bg)
}

motif.compute.enrich <- function (matches.se, compare = NULL, background = NULL) {
  # basically copied from ArchR
  if (is.null(compare)) {
    compare <- which(rowRanges(matches.se)$.type == "fg")
  }
  if (is.null(background)) {
    background <- 1:nrow(matches.se)
  }
  matches <- assay(matches.se)
  matchCompare <- matches[compare, , drop = FALSE]
  matchBackground <- matches[background, , drop = FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)
  pOut <- data.frame(feature = colnames(matches), CompareFrequency = matchCompareTotal, 
                     nCompare = nrow(matchCompare), CompareProportion = matchCompareTotal/nrow(matchCompare), 
                     BackgroundFrequency = matchBackgroundTotal, nBackground = nrow(matchBackground), 
                     BackgroundProporition = matchBackgroundTotal/nrow(matchBackground))
  pOut$Enrichment <- pOut$CompareProportion/pOut$BackgroundProporition
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x) {
    p <- -phyper(pOut$CompareFrequency[x] - 1, pOut$BackgroundFrequency[x], 
                 pOut$nBackground[x] - pOut$BackgroundFrequency[x], 
                 pOut$nCompare[x], lower.tail = FALSE, log.p = TRUE)
    return(p/log(10))
  }) %>% unlist %>% round(4)
  pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)), 
                          0)
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]
  pOut
}

n.vs.q <- function(np, cutoffs, titrate.col = "V9",
                   outfile = NULL) {
  if (is.character(np)) {
    np <- read.table(np, header = F, sep = "\t", quote =)
  }
  df <- lapply(cutoffs, function(cutoff) {
    n <- sum(np[, titrate.col] >= cutoff)
    df <- data.frame(cutoff = cutoff, n = n)
    return(df)
  }) %>% Reduce(rbind, .)
  if (!is.null(outfile)) {
    dir.create(dirname(outfile), showWarnings = F, recursive = T)
  }
  return(df)
}

motif.read.vierstra <- function(in.gr, region.col = "region", genome, nkc.only = F, cutoff) {
  if (nkc.only) {
    warning("looking at nkc only!!")
    motif.bed <- paste0("~/motifs/lodge/", genome, ".archetype_motifs.v1.0.nkc.bed.gz")
  } else {
    motif.bed <- paste0("~/motifs/lodge/", genome, ".archetype_motifs.v1.0.bed.gz")
  }
  if (is.null(mcols(in.gr)[[region.col]])) {
    mcols(in.gr)[[region.col]] <- paste0("region_", 1:length(in.gr))
  }
  mcols(in.gr) <- mcols(in.gr)[, region.col, drop = F]
  gr <- utilsFanc::import.bed.tabix(bed = motif.bed, gr = in.gr)
  gr <- plyranges::join_overlap_left(x = gr, y = in.gr)
  gr <- gr[gr$fifth >= cutoff]
  
  count.mat <- mcols(gr)[, c("forth", "region")] %>% as.data.frame() %>% 
    mutate(motif = motif.name.format(forth), presence = TRUE) %>% 
    reshape2::acast(formula = motif ~ region, value.var = "presence", fun.aggregate = sum)
  b.mat <- count.mat %>% as.logical.Array() %>% matrix(ncol = ncol(count.mat))
  dimnames(b.mat) <- dimnames(count.mat)
  res <- list(count.mat = count.mat, b.mat = b.mat)
  return(res)
}

motif.mat.filter <- function(motif.mat, rules.df, motifs.include = NULL) {
  # rules.df <- data.frame(rule = MAP1, include = "B6.Ly49h_MAP8:B6.Ly49m_MAP8", exclude = "B6.Ly49c_MAP8")
  if (!is.logical(motif.mat)) {
    stop("motif.mat must be logical")
  }
  if (!is.null(motifs.include)) {
    motifs.include <- motif.name.format(motifs.include)
    motif.mat <- motif.mat[rownames(motif.mat) %in% motifs.include, ]
  }
  rownames(motif.mat) <- motif.name.format(rownames(motif.mat))
  res <- rules.df %>% split(., f = 1:nrow(.)) %>% 
    lapply(function(rule) {
      fg <- rule$fg %>% strsplit(":") %>% unlist()
      bg <- rule$bg %>% strsplit(":") %>% unlist()
      utilsFanc::check.intersect(c(fg, bg), "region in rule", 
                                 colnames(motif.mat), "region in matif mat")
      
      bFg <- motif.mat[, fg, drop = F] %>% apply(1, all)
      bBg <- motif.mat[, bg, drop = F] %>% apply(1, any)
      bPass <- bFg & (!bBg)
      motifs.pass <- rownames(motif.mat)[bPass]
      return(motifs.pass)
    })
  names(res) <- rules.df$rule
  return(res)
}
# motif.diff <- function(in.gr, region.col = "region", genome, nkc.only = F, cutoff,
#                        rules.df) {
#   
# }

motif.name.format <- function(x) {
  res <- x %>% gsub("[^A-Za-z0-9]", "_", .)
  return(res)
}
