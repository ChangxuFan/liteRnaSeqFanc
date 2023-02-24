diffbind.from.target <- function(sample.info, file.out=NULL, target.pipe.dir,
                                 use.bed = T, zip.bed =F, bed2bam = T, force.bedTobam = F,
                                 genome, force.bam.index = F, max.threads = 6,
                                 bedToBam = "/bar/cfan/anaconda2/envs/jupyter/bin/bedToBam",
                                 samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  # sample.info has 2 columns: SampleID, Factor, Replicagte, alias
  # alias: a map between your convinient sample names to the sample names xiaoyun put for sequencing.
  if (use.bed == T) {
    reads.glob <- paste0("step3.3_*.open.bed")
  }
  else {
    reads.glob <- paste0("step2.1*.bam")
    warning("using the target bam file, which not deduped!!")
  }
  if (is.character(sample.info))
    sample.info <- read.table(sample.info, header = T, as.is = T)
  
  sample.info <- sample.info %>%
    mutate(bamReads = paste0(target.pipe.dir, "/*", alias, "*/", reads.glob) %>% Sys.glob() %>%  normalizePath(mustWork = T),
           Peaks = paste0(target.pipe.dir, "/*", alias, "*/step3.4*.narrowPeak") %>% Sys.glob() %>% normalizePath(mustWork = T),
           PeakCaller = "narrow",
           alias = NULL)
  
  sample.info$bamReads <- utilsFanc::safelapply(sample.info$bamReads, function(x) {
    if (use.bed == F) {
      if (!file.exists(paste0(x, ".bai")) || force.bam.index == T)  {
        cmd <- paste0("/bar/cfan/anaconda2/envs/jupyter/bin/samtools index ", x)
        print(cmd); system(cmd)
      }
    } else {
      if (bed2bam == T) {
        bam <- x %>% sub("bed$", "bam", .)
        if (!file.exists(bam) || force.bedTobam == T) {
          cmd <- paste0(bedToBam, " -i ", x, " -g ", "~/genomes/", genome, "/", genome, ".chrom.sizes",
                        " -mapq 42 > ", bam)
          print(cmd); system(cmd)
        }
        if (!file.exists(paste0(bam, ".bai")) || force.bam.index == T) {
          cmd <- paste0(samtools, " index ",bam )
          print(cmd); system(cmd)
        }
        return(bam)
      } else {
        if (!file.exists(paste0(x, ".gz")) && zip.bed == T) {
          cmd <- paste0("~/scripts/bed_browser_ez.sh ", x)
          print(cmd); system(cmd)
          return(paste0(x, ".gz"))
        }
      }
    }
    return(x)
  }, threads = min(max.threads, nrow(sample.info))) %>% unlist()
  
  if (use.bed == T && zip.bed == T) {
    sample.info$bamReads <- paste0(sample.info$bamReads, ".gz")
  }
  if (!is.null(file.out))
    write.table(sample.info, file.out, sep = "\t", row.names = F, col.names = T, quote = F)
  return(sample.info)
}

diffbind.write.report <- function(dbo, q.threshold=NULL, p.threshold=NULL, log2FC = NULL, 
                                  bSimple = F, bLabel.genes = T, genome,
                                  outdir = NULL, rootname = "") {
  if (is.null(q.threshold) && is.null(p.threshold))
    stop("at least one of q.threshold or p.threshold should be supplied")
  if (!is.null(q.threshold)) {
    th <- q.threshold
    bUsePval <- F
    flag <- "q"
  }
  if (!is.null(p.threshold)) {
    th <- p.threshold
    bUsePval <- T
    flag <- "p"
  }
  fold <- 0
  fold.flag <- ""
  if (!is.null(log2FC)) {
    fold <- log2FC
    fold.flag <- paste0("_f", log2FC)
  }
  contrast <- lapply(seq_along(dbo$contrasts), function(i) {
    options(scipen = 4)
    name1 <- dbo$contrasts[[i]]$name1
    name2 <- dbo$contrasts[[i]]$name2
    df <- DiffBind::dba.report(dbo, contrast = i, th = th, bUsePval = bUsePval, bNormalized = T,
                               bCalledDetail = T, bCalled = T, bCounts = T, fold = fold) %>% as.data.frame()
    colnames(df)[duplicated(colnames(df))] <- paste0(colnames(df)[duplicated(colnames(df))], "_called")
    
    df <- utilsFanc::add.column.fanc(df, data.frame(id = paste0("peak_", 1:nrow(df))#, notUsed = "miao"#, strand="."
                                                    ), pos = 4)
    df <- utilsFanc::add.column.fanc(df, data.frame(log10p.value = log10(df$p.value), log10FDR = log10(df$FDR)),
                                     after = "FDR")
    
    if (bLabel.genes) {
      gr <- df %>% makeGRangesFromDataFrame()
      gr <- utilsFanc::gr.fast.annot(gr, genome, use.klraps = F)
      df <- utilsFanc::add.column.fanc(df1 = df, df2 = data.frame(gene_closest = gr$SYMBOL),
                                       pos = 7)
    }
    
    
    if (!is.null(outdir)) {
      system(paste0("mkdir -p ", outdir))
      n.peaks <- lapply(c("all", "1up", "2up"), function(type) {
        if (type == "all") {
          file.name <- paste0(outdir, "/", rootname, "_",name1 , "_V_", name2, "_", flag, round(-log10(th), 2), fold.flag, ".bed")
        } else if (type == "1up") {
          file.name <- paste0(outdir, "/", rootname, "_",name1 , "_VV_", name2, "_", flag, round(-log10(th), 2), fold.flag, ".bed")
          df <- df %>% filter(df$Fold > 0)
        } else {
          file.name <- paste0(outdir, "/", rootname, "_",name2 , "_VV_", name1,"_", flag, round(-log10(th), 2), fold.flag, ".bed")
          df <- df %>% filter(df$Fold < 0)
        }
        n.peaks <- nrow(df)
        if (bSimple == T) {
          df <- df[, 1:4]
        }
        write.table(df, file.name,
                    row.names = F, col.names = T, quote = F, sep = "\t")
        system(paste0("/bar/cfan/scripts/bed_browser_v2.sh -s 1 ", file.name))
        if (nrow(df) > 100) {
          df <- df[sample(1:nrow(df), size = 100, replace = F) %>% sort(),] 
          file.name <- sub(".bed", "_peakwatch.bed", file.name)
          write.table(df, file.name,
                      row.names = F, col.names = T, quote = F, sep = "\t")
          system(paste0("/bar/cfan/scripts/bed_browser_v2.sh -s 1 ", file.name))
        }
        return(n.peaks)
      }) %>% unlist()
      stat.df <- data.frame(type = c("all", "up", "down"), 
                            n.peaks = n.peaks)
      write.table(stat.df, paste0(outdir, "/", rootname, "_stat.tsv"), sep = "\t", 
                  row.names = F, col.names = T, quote = F)
    }
    return(df)
  })
  return(contrast)
}


diffbind.fanc <- function(sampleSheet = NULL, dbo = NULL, cmd, 
                          attributes = DBA_FACTOR,
                          minOverlap = 0,
                          scoreCol, peak.filter = 0, q.th = 0.01,
                          # count params
                          file.type, summits = F, mapQCth = 30,
                          # DAR params:
                          block = NULL,
                          bParallel=T, ReportInit="DBA", 
                          bCorPlot = F, bLog.corr = T,
                          volcano.contrast = 1, volcano.log2FC = 1, 
                          volcano.flip = F, volcano.genes.label = NULL, genome,
                          rootname, outdir, ...) {
  # possible values for cmd: 
  ##s (sampleSheet, to create dbo from sampleSheet)
  ##c: count; a: contrast and analyse, r: report
  system(paste0("mkdir -p ", outdir))
  
  
  if ("s" %in% cmd) {
    if (is.character(sampleSheet))
      sampleSheet <- read.table(sampleSheet, as.is = T, sep = "\t", header = T, quote = "")
    # browser()
    dbo <- dba(sampleSheet = sampleSheet, minOverlap = minOverlap, scoreCol = scoreCol, filter = peak.filter,
               bLowerScoreBetter = F, bRemoveRandom = T, bRemoveM = T, bSummarizedExperiment = F, 
               config = data.frame(RunParallel = bParallel, DataType = DBA_DATA_GRANGES, 
                                   ReportInit= ReportInit, AnalysisMethod = DBA_DESEQ2,
                                   bCorPlot = bCorPlot, 
                                   bUsePval = F, minQCth = 30, th=q.th,
                                   fragments = F)
    )
  }
  
  if ("c" %in% cmd) {
    if (file.type == "bam") {
      bUseSummarizeOverlaps = T
      readFormat <- DBA_READS_BAM
      dbo$config$fragmentSize <- NULL
    } else {
      if (file.type == "bed") {
        bUseSummarizeOverlaps = F
        readFormat <- DBA_READS_BED
        dbo$config$fragmentSize <- 0
      } else {
        stop ("file.type has to be bam or bed")
      }
    }
    # browser()
    count.params <- list(DBA = dbo, minOverlap = minOverlap, score = DBA_SCORE_RPKM,  bLog = F, 
                   bRemoveDuplicates=F,
                   bScaleControl = F, mapQCth = mapQCth, filter = 0, filterFun=max,
                   #  bSubControl = F, minCount=0, not in my version
                   bUseSummarizeOverlaps = bUseSummarizeOverlaps,
                   bParallel=bParallel, readFormat = readFormat)
    if (file.type == "bed") {
      count.params <- c(params, list(summits = summits))
    }
    dbo <- do.call(dba.count, count.params)
    saveRDS(dbo, paste0(outdir, "/dbo_", rootname, ".Rds"))
  }
  if ("pca" %in% cmd)    {
    png(file=paste0(outdir, "/", rootname, "_pca.png"),
        height = 5, width = 5, units = "in", res = 200, pointsize = 5)
    try(print(dba.plotPCA(dbo, attributes = attributes, score = DBA_SCORE_RPKM, 
                          bLog = T, label = DBA_REPLICATE, cor = F)))
    dev.off()
  }
  
  if ("corr" %in% cmd ) {
    png(file=paste0(outdir, "/", rootname, "_corr.png"),
        height = 800, width = 800, units = "px", res = 150, pointsize = 10)
    try(print(dba.plotHeatmap(DBA = dbo, correlations=TRUE, score = DBA_SCORE_RPKM,
                              bLog = bLog.corr, distMethod="pearson")))
    dev.off()
    
  }
  
  if ("a" %in% cmd) {
    if (is.null(block)) {
      dbo <- dba.contrast(dbo, categories = DBA_FACTOR, minMembers = 2)
    } else {
      dbo <- dba.contrast(dbo, categories = DBA_FACTOR, minMembers = 2, block = block)
    }
    
    dbo <- dba.analyze(dbo, method = DBA_DESEQ2, bSubControl = F, bFullLibrarySize = F,
                       bTagwise = F, filter = 0, filterFun = max, bReduceObjects = F,
                       bParallel = bParallel)
    
    if (!is.null(block)) {
      dbo$config$AnalysisMethod <- DBA_DESEQ2_BLOCK
    }
    saveRDS(dbo, paste0(outdir, "/dbo_", rootname, ".Rds"))
  }
  
  if (!is.null(block)) {
    dbo$config$AnalysisMethod <- DBA_DESEQ2_BLOCK
  }
  
  if ("v" %in% cmd) {
    p <- dba.plotVolcano.fanc(dbo, contrast = volcano.contrast, th = q.th, 
                              genes.label = volcano.genes.label, genome = genome,
                         bUsePval = F, fold = volcano.log2FC, bFlip = volcano.flip)
    ggsave(paste0(outdir, "/", rootname, "_volcano.pdf"), p, width = 6, height = 5,  units = "in",
             dpi = 150, limitsize = F)
    scFanc::wrap.plots.fanc(list(p), plot.out = paste0(outdir, "/", rootname, "_volcano.html"), 
                            sub.width = 6, sub.height = 6, tooltip = "loci")
  }
  
  if ("m" %in% cmd) {
    png(file=paste0(outdir, "/", rootname, "_ma.png"),
        height = 800, width = 800, units = "px", res = 150, pointsize = 10)
    try(print(dba.plotMA(DBA = dbo, fold = volcano.log2FC, contrast = volcano.contrast, 
                         th = q.th, bUsePval = F, bFlip = volcano.flip)))
    dev.off()
    png(file=paste0(outdir, "/", rootname, "_xy.png"),
        height = 800, width = 800, units = "px", res = 150, pointsize = 10)
    try(print(dba.plotMA(DBA = dbo, fold = volcano.log2FC, contrast = volcano.contrast, 
                         th = q.th, bUsePval = F, bFlip = volcano.flip, bXY = T)))
    dev.off()
  }
  
  if ("r" %in% cmd) {
    diffbind.write.report(dbo = dbo, q.threshold = q.th, outdir = paste0(outdir, "/report/"),
                          log2FC = volcano.log2FC, rootname = rootname, genome = genome)
  }
  return(dbo)
}


# diffbind.fanc.bk <- function(samples=NULL, dbo = NULL, 
#                              rootname, outdir, ...) {
#   # ...: additional commands to pass on to diffbind.
#   stop("function untested!!")
#   if (is.null(dbo)) {
#     if (!is.null(samples)) {
#       if (!is.data.frame(samples))
#         samples <- read.table(samples, as.is = T, header = T)
#     } else
#       stop("either samples (df or file location) or an already counted diffbind object needs to be offered")
#     dba.all.pre <- dba(sampleSheet = samples, minOverlap = 0, 
#                        bLowerScoreBetter = F, bRemoveRandom = T, bRemoveM = T, bSummarizedExperiment = F, 
#                        config = data.frame(DataType = DBA_DATA_FRAME, AnalysisMethod = DBA_DESEQ2,
#                                            bUsePval = F, minQCth = 30, th=0.01))
#     dba.all <- dba.count(dba.all.pre, score = DBA_SCORE_READS, minOverlap = 0, bLog = F, fragmentSize = 0, summits = F, bRemoveDuplicates=F,
#                          bScaleControl = F, mapQCth = 30, filter = 0, bParallel=T)
#     
#   } else {
#     dba.all <- dbo
#   }
#   
#   png(file=paste0(outdir, "/", rootname, "_pca.png"),
#       height = 3000, width = 3000, units = "px")
#   dba.plotPCA(dba.all, attributes = DBA_FACTOR, score = DBA_SCORE_RPKM, bLog = T, label = DBA_REPLICATE, cor = F)
#   dev.off()
#   
#   dba.all.contrast <- dba.contrast(dba.all, categories = DBA_FACTOR)
#   dba.all.analyse <- dba.analyze(dba.all.contrast, method = DBA_DESEQ2, bSubControl = F, bFullLibrarySize = F,
#                                  bTagwise = F, filter = 0, filterFun = max, bReduceObjects = F,
#                                  bParallel = T)
#   
#   dba.all.report <- diffbind.write.report(dba.all.analyse, q.threshold = 0.01, outdir = "diffbind/all_target/", rootname = "all_target")
#   
#   return(list(dba = dba.all.analyse, report=dba.all.report))
# }

pv.DBAplotVolcano.fanc <- function (pv, contrast, method = "edgeR", th = 0.05, bUsePval = F, 
                                    fold = 0, facname = "", bLabels = FALSE, maxLabels = 50,
                                    genes.label = NULL, genome, 
                                    dotSize = 1, bSignificant = T, bFlip = FALSE, xrange, yrange) 
{
  if (missing(contrast)) {
    contrast <- 1:length(pv$contrasts)
  }
  else {
    if (contrast > length(pv$contrasts)) {
      stop("Specified contrast number is greater than number of contrasts")
      return(NULL)
    }
  }
  # added by FANC: 
  p.list <- list()
  for (con in 1:length(contrast)) {
    conrec <- pv$contrasts[[contrast[con]]]
    name1 <- conrec$name1
    name2 <- conrec$name2
    if (bFlip) {
      name1 <- conrec$name2
      name2 <- conrec$name1
    }
    for (meth in method) {
      res <- pv.DBAreport(pv, contrast = contrast[con], 
                          method = meth, bUsePval = T, th = 100, bNormalized = TRUE, 
                          bFlip = bFlip, precision = 0)
      if (!is.null(res)) {
        if (bUsePval) {
          vals <- res$"p-value"
          res$vals <- res$`p-value`
          idx <- vals <= th
          tstr <- "p"
          res = mutate(res, Legend = ifelse(res$"p-value" <= 
                                              th, sprintf(" p-val<=%1.2f", th), sprintf(" p-val >%1.2f", 
                                                                                        th)))
        }
        else {
          vals <- res$FDR
          res$vals <- res$FDR
          idx <- vals <= th
          tstr <- "FDR"
          res = mutate(res, Legend = ifelse(res$FDR < 
                                              th, sprintf(" FDR<=%1.2f", th), sprintf(" FDR >%1.2f", 
                                                                                      th)))
        }
        res$Legend[idx & abs(res$Fold) < fold] <- sprintf("abs(Fold)<%1.2f", 
                                                          2^fold)
        idx <- idx & abs(res$Fold) >= fold
        sigSites <- res[idx, ]
        rownames(sigSites) <- 1:sum(idx)
        res <- cbind(0, res)
        colnames(res)[1] <- "SiteNum"
        res[idx, 1] <- 1:sum(idx)
        plotTitle <- sprintf("%s Contrast: %s vs. %s [%s %s<=%1.3f", 
                             facname, name1, name2, sum(idx), tstr, th)
        if (fold > 0) {
          plotTitle <- sprintf("%s & abs(Fold)>=%1.2f]", 
                               plotTitle, 2^fold)
        }
        else {
          plotTitle <- sprintf("%s]", plotTitle)
        }
        xLabel <- sprintf("log2(%s/%s)", name1, 
                          name2)
        yLabel <- sprintf("-log10(%s)", tstr)

        res$Legend[res$Legend == " FDR<=0.05"] <- "DARs"
        res$Legend[res$Legend != "DARs"] <- "nonDARs"
        res$loci <- paste0(res$Chr, ":", res$Start, "-", res$End)
        p <- ggplot(res, aes(x = Fold, y = -log10(vals), text = loci)) + 
          geom_point(aes(col = Legend), size = dotSize) + 
          scale_color_manual(values = c(crukMagenta, crukGrey)) + 
          labs(title = plotTitle, x = xLabel, y = yLabel)

        if (!is.null(genes.label)) {
          res.wGene <- diffbind.res.add.gene(res.df = res, genome = genome)
          res.wGene <- res.wGene %>% dplyr::filter(gene_closest %in% genes.label,
                                            Legend == "DARs")
          if (nrow(res.wGene) < 1) {
            warning("none of the genes you wanted to label is in the dataset")
          } else {
            p <- p + ggrepel::geom_text_repel(data = res.wGene, inherit.aes = F,
                                 aes(x = Fold, y = -log10(vals), label = gene_closest),
                                 min.segment.length = 0, force = 1, point.padding = 0,
                                 box.padding = 2
                                 )
            # ggplot(data.frame(x = 1:3, y = 2:4), aes(x = x, y = y)) +
            #   geom_point() +
            #   ggrepel::geom_text_repel(aes(label = y), min.segment.length = 0, 
            #                            box.padding = 1.5, point.padding = 0, force = 1)
            
          }
        } else if (bLabels) {
          maxLabels <- min(sum(idx), maxLabels)
          if (maxLabels > 0) {
            xx <- which(idx)[1:maxLabels]
            p <- p + geom_text_repel(data = sigSites[1:maxLabels, 
                                                     ], aes(x = Fold, y = -log10(vals[xx]), 
                                                            label = rownames(sigSites)[1:maxLabels]))
          }
        }
      }
    }
  }
  all.font.size <- 7.5
  p <- p  +
    theme(plot.title = element_text(size=10 ), legend.title = element_text(size = all.font.size ),
          legend.text = element_text(size = all.font.size ), axis.title = element_text(size = all.font.size )) +
    theme_bw()
    
  return(p)
}


dba.plotVolcano.fanc <- function(DBA, contrast = 1, method = DBA$config$AnalysisMethod, 
                                 th = DBA$config$th, bUsePval = DBA$config$bUsePval, fold = 0, 
                                 factor = "", bFlip = FALSE, bLabels = FALSE, maxLabels = 50, 
                                 genes.label = NULL, genome,
                                 dotSize = 1) {
  DBA <- pv.check(DBA, bCheckEmpty = TRUE)
  
  res <- pv.DBAplotVolcano.fanc(DBA, contrast = contrast, method = method, 
                                th = th, bUsePval = bUsePval, fold = fold, facname = factor, 
                                dotSize = dotSize, bFlip = bFlip, bLabels = bLabels, 
                                genes.label = genes.label, genome = genome,
                                maxLabels = maxLabels)
  return(res)
}

# add gene names to volcano plots
diffbind.res.add.gene <- function(res.df, genome) {
  res.gr <- makeGRangesFromDataFrame(res.df, keep.extra.columns = T)
  if (genome == "mm10") {
    genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
    annoDb <- "org.Mm.eg.db"
    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  } else if (genome == "hg38") {
    genome.name <- "BSgenome.Hsapiens.UCSC.hg38"
    annoDb <- "org.Hs.eg.db"
    TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  } else {
    stop("only mm10 and hg38 supported")
  }
  anno <- ChIPseeker::annotatePeak(res.gr, TxDb = TxDb, annoDb = annoDb)@anno
  if (length(anno) != length(res.gr)) {
    stop("length(anno) != length(res.gr)")
  }
  if (!identical(res.df$SiteNum, res.gr$SiteNum)) {
    stop("!identical(anno$SiteNum, res.gr$SiteNum)")
  }
  res.df$gene_closest <- anno$SYMBOL
  return(res.df)
}


diffbind.fast.report <- function(gr, out.dir, root.name) {
  utilsFanc::write.zip.fanc(df = gr[gr$Fold > 0,], 
                            out.file = paste0(out.dir, "/", root.name, "_pos.bed"),
                            bed.shift = T, zip = T)
  utilsFanc::write.zip.fanc(df = gr[gr$Fold < 0,], 
                            out.file = paste0(out.dir, "/", root.name, "_neg.bed"),
                            bed.shift = T, zip = T)
  return()
}

dba.plotMA.fanc <- function (DBA, peaks.gr = NULL, highlight.idx = NULL, highlight.gr = NULL,
                             plot.out = NULL, height = 800, width = 800, 
                             contrast = 1, method = DBA$config$AnalysisMethod, 
                             th = DBA$config$th, bUsePval = DBA$config$bUsePval, fold = 0, 
                             bNormalized = TRUE, factor = "", bFlip = FALSE, bXY = FALSE, 
                             dotSize = 0.45, bSignificant = TRUE, bSmooth = TRUE, xrange, 
                             yrange, ...) 
{
  DBA <- pv.check(DBA, bCheckEmpty = TRUE)
  if (!is.null(plot.out)) {
    png(file= plot.out,
        height = height, width = width, units = "px", res = 150, pointsize = 10)
  }
  try(print(res <- pv.DBAplotMA.fanc(DBA, peaks.gr = peaks.gr, highlight.idx = highlight.idx, highlight.gr = highlight.gr,
                                     contrast = contrast, method = method, 
                                     bMA = !bXY, bXY = bXY, th = th, bUsePval = bUsePval, 
                                     fold = fold, facname = factor, bNormalized = bNormalized, 
                                     cex = dotSize, bSignificant = bSignificant, bSmooth = bSmooth, 
                                     bFlip = bFlip, xrange = xrange, yrange = yrange, ...)))
  
  if (!is.null(plot.out)) {
    dev.off()
  }
  invisible(res)
}

dba.plotXY.fanc <- function (DBA, peaks.gr = NULL, highlight.idx = NULL, highlight.gr = NULL,
                             plot.out = NULL, color.density = F, quantile.limit = NULL,
                             pt.size = 0.05,  transformation = function(x) return(2^x),
                             contrast = 1, method = DBA$config$AnalysisMethod, 
                             th = DBA$config$th, bUsePval = DBA$config$bUsePval, fold = 0, 
                             bNormalized = TRUE, factor = "", bFlip = FALSE, bXY = T, 
                             dotSize = 0.45, bSignificant = TRUE, bSmooth = TRUE, xrange, 
                             yrange, return.data = F,
                             ...) 
{
  DBA <- pv.check(DBA, bCheckEmpty = TRUE)

  res <- pv.DBAplotMA.fanc(DBA, peaks.gr = peaks.gr, highlight.idx = highlight.idx, highlight.gr = highlight.gr,
                           contrast = contrast, method = method, 
                           bMA = !bXY, bXY = bXY, th = th, bUsePval = bUsePval, 
                           fold = fold, facname = factor, bNormalized = bNormalized, 
                           cex = dotSize, bSignificant = bSignificant, bSmooth = bSmooth, 
                           bFlip = bFlip, xrange = xrange, yrange = yrange, return.xy.df = T, ...)
  if (return.data == F) {
    df <- res$res
    x.name <- names(df)[6]
    y.name <- names(df)[5]
    df$hl <- 0
    df$hl[res$highlight.idx] <- 1

    scFanc::xy.plot(df = df, x = x.name, y = y.name, transformation = transformation, 
                    quantile.limit = quantile.limit, highlight.var = "hl", highlight.values = 1, outfile = plot.out,
                    color.density = color.density, pt.size = pt.size)
    return(NULL)
  } else {
    return(res)
  }
  
}

pv.DBAplotMA.fanc <- function (pv, peaks.gr = NULL, highlight.idx = NULL, highlight.gr = NULL, 
                               contrast, method = "edgeR", bMA = T, bXY = F, th = 0.05, 
                               bUsePval = F, fold = 0, facname = "", bNormalized = T, cex = 0.15, 
                               bSignificant = T, bSmooth = T, bFlip = FALSE, xrange, yrange, 
                               return.xy.df = F,
                               ...) 
{
  if (missing(contrast)) {
    contrast <- 1:length(pv$contrasts)
  }
  else {
    if (contrast > length(pv$contrasts)) {
      stop("Specified contrast number is greater than number of contrasts")
      return(NULL)
    }
  }
  plotfun <- plot
  if (bSmooth) {
    plotfun <- smoothScatter
  }
  numSites <- nrow(pv$binding)
  for (con in 1:length(contrast)) {
    conrec <- pv$contrasts[[contrast[con]]]
    name1 <- conrec$name1
    name2 <- conrec$name2
    if (bFlip) {
      name1 <- conrec$name2
      name2 <- conrec$name1
    }
    for (meth in method) {
      ### written by FANC
      if (is.null(peaks.gr)) {
        res <- pv.DBAreport(pv, contrast = contrast[con], 
                            method = meth, bUsePval = T, th = 100, bNormalized = bNormalized, 
                            bFlip = bFlip, precision = 0)
      } else {
        res <- peaks.gr %>% as.data.frame() %>% dplyr::rename(Chr = seqnames, Start = start, End = end) %>% 
          dplyr::mutate(strand = NULL)
      }

      if (!is.null(res)) {
        ### written by FANC
        if (!is.null(highlight.idx)) {
          if (is.logical(highlight.idx))
            idx <- highlight.idx
          else
            idx <- (1:nrow(res)) %in% highlight.idx
        } else if (!is.null(highlight.gr)) {
          gr <- makeGRangesFromDataFrame(res)
          o.fanc <- findOverlaps(query = gr, subject = highlight.gr)
          id.n <- queryHits(o.fanc) %>% unique()
          idx <- (1:nrow(res)) %in% id.n
        } else {
          if (bUsePval) {
            idx <- res$"p-value" <= th
            # tstr <- "p"
          }
          else {
            idx <- res$FDR <= th
            # tstr <- "FDR"
          }
          idx <- idx & (abs(res$Fold) >= fold)
        }
        ### written by FANC
        if (bUsePval)
          tstr <- "p"
        else
          tstr <- "FDR"
        ###
        if (bMA) {
          if (missing(xrange)) {
            xmin <- floor(min(res$Conc))
            xmax <- ceiling(max(res$Conc))
          }
          else {
            if (length(xrange) != 2) {
              stop("xrange must be vector of two numbers")
            }
            xmin <- xrange[1]
            xmax <- xrange[2]
          }
          if (missing(yrange)) {
            ymin <- floor(min(res$Fold))
            ymax <- ceiling(max(res$Fold))
          }
          else {
            if (length(yrange) != 2) {
              stop("yrange must be vector of two numbers")
            }
            ymin <- yrange[1]
            ymax <- yrange[2]
          }
          if (bSmooth | !bSignificant) {
            plotfun(res$Conc, res$Fold, pch = 20, cex = cex, 
                    col = crukBlue, xaxp = c(xmin, xmax, xmax - 
                                               xmin), xlim = c(xmin, xmax), xlab = "log concentration", 
                    yaxp = c(ymin, ymax, (ymax - ymin)), ylim = c(ymin, 
                                                                  ymax), ylab = sprintf("log fold change: %s - %s", 
                                                                                        name1, name2), main = sprintf("%s Binding Affinity: %s vs. %s (%s %s < %1.3f)", 
                                                                                                                      facname, name1, name2, sum(idx), tstr, 
                                                                                                                      th), ...)
          }
          else {
            plotfun(res$Conc[!idx], res$Fold[!idx], pch = 20, 
                    cex = cex, col = crukBlue, xaxp = c(xmin, 
                                                        xmax, xmax - xmin), xlim = c(xmin, xmax), 
                    xlab = "log concentration", yaxp = c(ymin, 
                                                         ymax, (ymax - ymin)), ylim = c(ymin, 
                                                                                        ymax), ylab = sprintf("log fold change: %s - %s", 
                                                                                                              name1, name2), main = sprintf("%s Binding Affinity: %s vs. %s (%s %s < %1.3f)", 
                                                                                                                                            facname, name1, name2, sum(idx), tstr, 
                                                                                                                                            th), ...)
          }
          if (bSignificant) {
            points(res$Conc[idx], res$Fold[idx], pch = 20, 
                   cex = cex, col = crukMagenta)
          }
          abline(h = 0, col = "dodgerblue")
        }
        if (bXY) {
          if (missing(xrange)) {
            xmin <- floor(min(res[, 5]))
            xmax <- ceiling(max(res[, 5]))
          }
          else {
            if (length(xrange) != 2) {
              stop("xrange must be vector of two numbers")
            }
            xmin <- xrange[1]
            xmax <- xrange[2]
          }
          if (missing(yrange)) {
            ymin <- floor(min(res[, 6]))
            ymax <- ceiling(max(res[, 6]))
          }
          else {
            if (length(yrange) != 2) {
              stop("yrange must be vector of two numbers")
            }
            ymin <- yrange[1]
            ymax <- yrange[2]
          }
          xymin <- min(xmin, ymin)
          xymin <- max(xymin, 0)
          xymax <- max(xmax, ymax)
          if (return.xy.df == T) {
            return(list(res = res, highlight.idx = idx))
          } else {
            plotfun(res[!idx, 6], res[!idx, 5], pch = 20, 
                    cex = cex, col = crukBlue, xaxp = c(xymin, 
                                                        xymax, xymax - xymin), xlim = c(xymin, 
                                                                                        xymax), xlab = sprintf("log concentration :%s", 
                                                                                                               name2), yaxp = c(xymin, xymax, (xymax - 
                                                                                                                                                 xymin)), ylim = c(xymin, xymax), ylab = sprintf("log concentration :%s", 
                                                                                                                                                                                                 name1), main = sprintf("%s Binding Affinity: %s vs. %s (%s %s < %1.3f)", 
                                                                                                                                                                                                                        facname, name1, name2, sum(idx), tstr, 
                                                                                                                                                                                                                        th), ...)
            points(res[idx, 6], res[idx, 5], pch = 20, 
                   cex = cex, col = crukMagenta)
            abline(0, 1, col = "dodgerblue")
          }

        }
      }
    }
  }
}

scale.fanc <- function(mat, margin, max = NULL, min = NULL) {
  mat <- as.matrix(mat)
  if (margin == 1) {
    mat <- t(mat)
  }
  
  mat <- scale(mat, center = T, scale = T)
  if (!is.null(max))
    mat[mat > max] <- max
  if (!is.null(min))
    mat[mat < min] <- min
  return(mat)
  
}

diffbind.deeptools <- function(dbo, bw.vec, method, q.th, log2FC,
                               genome, # for blacklist purpose
                               sort.using.samples = NULL,
                               do.compute.matrix = T, do.plot.heatmap = T,
                               threads = 1, work.dir, root.name = NULL,
                               debug = F) {
  dbo$config$AnalysisMethod <- method
  
  if (is.null(root.name)) {
    root.name <- paste0("q", -log10(q.th),"_f", log2FC )
  }
  peak.dir <- paste0(work.dir, "/peaks/", root.name, "/")
  # system(paste0("rm -rf ", peak.dir))
  blacklist <- paste0("~/genomes/", genome, "/blacklist/", genome, ".blacklist.bed")
  if (!file.exists(blacklist)) {
    stop("!file.exists(blacklist)")
  }
  
  trash <- diffbind.write.report(dbo = dbo, q.threshold = q.th, 
                                 outdir = peak.dir,
                                 log2FC = log2FC, bSimple = T)
  regions.vec <- Sys.glob(paste0(peak.dir, "/*_VV_*.bed"))
  lapply(regions.vec, function(bed) {
    system(paste0("sed -i '1d' ", bed))
  })
  
  if (debug == T) {
    regions.vec <- regions.vec %>% .[grepl("peakwatch",.)]
  } else {
    regions.vec <- regions.vec %>% .[!grepl("peakwatch",.)]
  }
  
  deeptools.refpoint(bw.vec = bw.vec, regions.vec = regions.vec, 
                     blacklist = blacklist, sort.using.samples = sort.using.samples, 
                     compute.matrix = do.compute.matrix, 
                     plot.heatmap = do.plot.heatmap, threads = threads, 
                     out.dir = work.dir, root.name = root.name)
  return()

}

diffbind.homer <- function(dbo, report.bed, n.peaks.fg = 1000, n.peaks.bg = 5000,
                           trends = c("up", "down", "de"), n.sets.bg = 3,
                           fg.seed = 10, bg.seeds = c(42, 84, 126),
                           out.dir, threads.each = 4, npar = 1) {
  
}

diffbind.make.a2bl <- function(dbo, contrast.id, dar.bed, root.name, force = F,
                               res.slot = "res_diffbind", summary.slot = "summary_diffbind") {
  # no longer an issue: somehow the blocked versions come with colnames!
  # stop("wrong: columns of dds are often rearranged!")
  # caveat is that non-blocked versions are not supported~
  dds <- dbo$contrasts[[contrast.id]]$DESeq2$block$DEdata
  if (is.null(dds)) {
    if (force) {
      dds <- dbo$contrasts[[contrast.id]]$DESeq2$DEdata
      warning("force mode. There might be mismatch of samples!!")
    } else {
      stop("is.null(dds)")
    }
    
  }
  rownames(dds) <- dba.peakset(dbo, bRetrieve = T) %>% utilsFanc::gr.get.loci()
  if (force) {
    colnames(dds) <- dbo$samples$SampleID
  }
  res <- dba.report(dbo, DataType = DBA_DATA_GRANGES, th = 1)
  res$gene <- utilsFanc::gr.get.loci(res)
  res <- mcols(res) %>% as.data.frame()
  res <- res %>% dplyr::rename(log2FoldChange = Fold, pvalue = p.value, padj = FDR)
  dar <- read.table(dar.bed, header = T)
  dar$gene <- paste0(dar[, 1], ":", dar[, 2], "-", dar[, 3])
  summ <- list(de.genes = dar$gene, up.genes = dar$gene[dar$Fold > 0], down.genes = dar$gene[dar$Fold < 0])
  summ.n <- lapply(summ, length)
  names(summ.n) <- c("n.de", "n.up", "n.down")
  summ <- c(summ, summ.n)
  
  a2b <- list(root.name = root.name, dds = dds)
  a2b[[res.slot]] <- res
  a2b[[summary.slot]] <- summ
  a2bl <- list()
  a2bl[[root.name]] <- a2b
  return(a2bl)
}

diffbind.get.raw.mat <- function(dbo) {
  stop("wrong: columns of dds are often rearranged!")
  n.covered <- sum(dbo$contrasts[[1]]$group1, dbo$contrasts[[1]]$group2)
  if (n.covered != nrow(dbo$samples)) {
    stop("n.covered != nrow(dbo$samples)")
  }
  mat <- dbo$contrasts[[1]]$DESeq2$counts
  loci <- dba.peakset(dbo, bRetrieve = T, bRemoveRandom = T, bRemoveM = T) %>% 
    utilsFanc::gr.get.loci()
  rownames(mat) <- loci
  colnames(mat) <- dbo$samples$SampleID
  return(mat)
}



diffbind.match.bg <- function(dbo, peakset.list = NULL, contrast = 1, method, q.th, log2FC,
                              n = 3, mask = NULL) {
  if (is.null(peakset.list)) {
    peakset <- dba.report(DBA = dbo, contrast = contrast, method = method, th = q.th, 
                          fold = log2FC)
    peakset.list <- list(up = peakset %>% plyranges::filter(Fold > 0),
                         down = peakset %>% plyranges::filter(Fold < 0))
    rm(peakset)
  }
  gr <- dba.peakset(dbo, bRetrieve = T, bRemoveM = F, bRemoveRandom = F)
  df <- mcols(gr) %>% as.data.frame()
  rownames(df) <- gr %>% utilsFanc::gr.get.loci()
  
  if (!is.null(mask)) {
    use.mask <- dbo$masks[[mask]]
    if (is.null(use.mask)) {
      stop("is.null(use.mask)")
    }
    df <- df[, use.mask]
    
  }
  mat <- df %>% as.matrix()
  
  bg.list <- lapply(peakset.list, function (peakset) {
    fg <- peakset %>% utilsFanc::gr.get.loci()
    fg.mat <- mat[fg, ]
    search.space.mat <- mat[! rownames(mat) %in% fg, ]
    
    df$locus <- rownames(df)
    df$group <- "bg"
    df$group[df$locus %in% fg] <- "fg"
    t <- ArchR::.matchBiasCellGroups(input = df, groups = df$group, useGroups = "fg", bgdGroups = "bg",
                                bias = c("A4WT", "A6WT"), k = 100, n = 10000, bufferRatio = 0.8, logFile = NULL)
    for (i in 1:n) {
      
      t <- Biobase::matchpt(x = fg.mat, y = search.space.mat)
    }
    # bg <- scFanc::bg.gen.2(mat = df, fg.vec = fg, n.bg.each = n, no.replace = T)
  })
}


