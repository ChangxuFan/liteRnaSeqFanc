deeptools.refpoint <- function (bw.vec, regions.vec, blacklist = NULL, 
                                upstream = 2000, downstream = 2000, norm.to = NULL, binSize = 10,
                                missingDataAsZero = T, missingDataColor = 1, 
                                color.map = "Reds", sort.using.samples = 1, compute.matrix = T, 
                                plot.heatmap = T, threads, out.dir, root.name = NULL, 
                                other.params.computeMatrix = "",
                                other.params.plotHeatmap = "",
                                computeMatrix = "/bar/cfan/anaconda2/envs/jupyter/bin/computeMatrix", 
                                plotHeatmap = "/bar/cfan/anaconda2/envs/jupyter/bin/plotHeatmap") {
  # copied from scFanc::deeptools.refpoint. 
  if (!is.null(names(bw.vec))) {
    sample.labels <- names(bw.vec)
  } else {
    sample.labels <- basename(bw.vec) %>% 
      sub("\\.bw$|\\.bigwig$|\\.bigWig$|\\.BigWig$", "", .)
  }
  if (!is.null(names(regions.vec))) {
    region.labels <- names(regions.vec)
  } else {
    region.labels <- basename(regions.vec) %>% sub(".bed|.narrowPeak", "", .)
  }
  if (length(color.map) != 1) {
    if (length(color.map) != length(bw.vec)) {
      stop("length(color.map) != length(bw.vec)")
    }
  }
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  dir.create(out.dir, showWarnings = F, recursive = T)
  prefix <- paste0(out.dir, "/", root.name)
  mat <- paste0(prefix, ".mat.gz")
  mat_hm <- paste0(prefix, ".tsv.gz")
  sorted_bed <- paste0(prefix, ".sorted.bed")
  heatmap <- paste0(prefix, ".hm.pdf")
  if (compute.matrix == T) {
    cmd <- paste0(computeMatrix, " reference-point -S ", 
                  paste0(bw.vec, collapse = " "), " -R ", paste0(regions.vec, collapse = " "),
                  " -o ", mat, 
                  ifelse(missingDataAsZero, " --missingDataAsZero ", ""), 
                  " -b ", downstream, " -a ", upstream, " --referencePoint center ", 
                  " --binSize ", binSize,
                  " --samplesLabel ", paste0(sample.labels, collapse = " "), 
                  " -p ", threads, 
                  other.params.computeMatrix)
    if (!is.null(blacklist)) {
      cmd <- paste0(cmd, " -bl ", blacklist)
    }
    print(cmd)
    system(cmd)
  }
  
  if (!is.null(norm.to)) {
    mat <- deeptools.mat.norm(mat.file = mat, norm.to = norm.to)
  }
  if (plot.heatmap == T) {
    cmd <- paste0(plotHeatmap, " -m ", mat, " -out ", heatmap, 
                  " --outFileNameMatrix ", mat_hm, " --outFileSortedRegions ", 
                  sorted_bed, " --colorMap ", color.map, " --refPointLabel center",
                  " --missingDataColor ", missingDataColor,
                  " --regionsLabel ", paste0(region.labels, collapse = " "),
                  paste0(other.params.plotHeatmap, collapse = " "))
    if (!is.null(sort.using.samples)) {
      cmd <- paste0(cmd, " --sortUsingSamples ", sort.using.samples)
    }
    print(cmd)
    system(cmd)
  }
  return()
}

deeptools.scaleRegion <- function (bw.vec, regions.vec, blacklist = NULL, 
                                upstream = 2000, downstream = 2000, norm.to = NULL, binSize = 10,
                                missingDataAsZero = T, missingDataColor = 1, 
                                color.map = "Reds", sort.using.samples = 1, compute.matrix = T, 
                                plot.heatmap = T, threads, out.dir, root.name = NULL, 
                                other.params.computeMatrix = "",
                                other.params.plotHeatmap = "",
                                computeMatrix = "/bar/cfan/anaconda2/envs/jupyter/bin/computeMatrix", 
                                plotHeatmap = "/bar/cfan/anaconda2/envs/jupyter/bin/plotHeatmap") {
  stop("just copied from deeptools.refpoint. Haven't finished modifying yet")
  if (!is.null(names(bw.vec))) {
    sample.labels <- names(bw.vec)
  } else {
    sample.labels <- basename(bw.vec) %>% 
      sub("\\.bw$|\\.bigwig$|\\.bigWig$|\\.BigWig$", "", .)
  }
  if (!is.null(names(regions.vec))) {
    region.labels <- names(regions.vec)
  } else {
    region.labels <- basename(regions.vec) %>% sub(".bed|.narrowPeak", "", .)
  }
  if (length(color.map) != 1) {
    if (length(color.map) != length(bw.vec)) {
      stop("length(color.map) != length(bw.vec)")
    }
  }
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  dir.create(out.dir, showWarnings = F, recursive = T)
  prefix <- paste0(out.dir, "/", root.name)
  mat <- paste0(prefix, ".mat.gz")
  mat_hm <- paste0(prefix, ".tsv.gz")
  sorted_bed <- paste0(prefix, ".sorted.bed")
  heatmap <- paste0(prefix, ".hm.pdf")
  if (compute.matrix == T) {
    cmd <- paste0(computeMatrix, " scale-regions -S ", 
                  paste0(bw.vec, collapse = " "), " -R ", paste0(regions.vec, collapse = " "),
                  " -o ", mat, 
                  ifelse(missingDataAsZero, " --missingDataAsZero ", ""), 
                  " -b ", downstream, " -a ", upstream, " --referencePoint center ", 
                  " --binSize ", binSize,
                  " --samplesLabel ", paste0(sample.labels, collapse = " "), 
                  " -p ", threads, 
                  other.params.computeMatrix)
    if (!is.null(blacklist)) {
      cmd <- paste0(cmd, " -bl ", blacklist)
    }
    print(cmd)
    system(cmd)
  }
  
  if (!is.null(norm.to)) {
    mat <- deeptools.mat.norm(mat.file = mat, norm.to = norm.to)
  }
  if (plot.heatmap == T) {
    cmd <- paste0(plotHeatmap, " -m ", mat, " -out ", heatmap, 
                  " --outFileNameMatrix ", mat_hm, " --outFileSortedRegions ", 
                  sorted_bed, " --colorMap ", color.map, " --refPointLabel center",
                  " --missingDataColor ", missingDataColor,
                  " --regionsLabel ", paste0(region.labels, collapse = " "),
                  paste0(other.params.plotHeatmap, collapse = " "))
    if (!is.null(sort.using.samples)) {
      cmd <- paste0(cmd, " --sortUsingSamples ", sort.using.samples)
    }
    print(cmd)
    system(cmd)
  }
  return()
}




deeptools.pileup.qc <- function(bw.vec, all.in.one = T, genome, plot.type = "TSS", 
                                # if plot.type == "motif"
                                motif.name,
                                n.regions = 10000, n.trials = 1, 
                                threads.master = NULL, threads.each = 1,
                                compute.matrix = T, 
                                out.dir, root.name = NULL, ...) {
  dir.create(out.dir, showWarnings = F, recursive = T)
  if (is.null(threads.master))
    threads.master <- length(bw.vec)
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  if (plot.type %in% c("TSS", "motif")) {
    if (plot.type == "TSS") {
      df <- read.table(paste0("~/genomes/", genome, "/", genome, "_TSS.bed"), sep = "\t")
      out.suffix <- plot.type
    } else {
      df <- read.table(paste0("~/motifs/sth/archetypes/", genome, "/", motif.name, ".bed"), sep = "\t")
      out.suffix <- paste0(plot.type, "_", motif.name)
    }
    regions.vec <- utilsFanc::safelapply(1:n.trials, function(i) {
      seed <- i * 10
      region.file <- paste0(out.dir, "/", root.name, "_", out.suffix, "_seed_", seed, ".bed")
      set.seed(seed = seed)
      df <- df %>% .[sample(1:nrow(.), n.regions, F),]
      write.table(df, region.file, quote = F, sep = "\t", row.names = F, col.names = F)
      return(region.file)
    }, threads = n.trials) %>% unlist()
    blacklist <- paste0("~/genomes/", genome, "/blacklist/", genome, ".blacklist.bed")
    
    if (all.in.one == T) {
      bw <- list(bw.vec)
      deeptools.threads <- threads.master
    } else {
      bw <- bw.vec
      deeptools.threads <- threads.each
      root.name <- names(bw) %||% (basename(bw) %>% sub("\\.bw$|\\.bigwig$|\\.bigWig$|\\.BigWig$", "", .))
    }
    utilsFanc::safelapply(seq_along(bw), function(i) {
      if (is.list(bw)) {
        x <- bw[[i]]
      } else {
        x <- bw[i]
      }
      root.name <- root.name[i]
      deeptools.refpoint(bw.vec = x, regions.vec = regions.vec, blacklist = blacklist, 
                         upstream = 2000, downstream = 2000, 
                         color.map = "Reds", sort.using.samples = NULL, 
                         compute.matrix = compute.matrix, 
                         plot.heatmap = T,
                         threads = deeptools.threads, 
                         out.dir = out.dir, root.name = root.name, ...)
    }, threads = threads.master)
    
  } else {
    stop("only TSS has been developed")
  }
}

deeptools.mat.read <- function(mat.file) {
  # mat.file <- "/bar/cfan/4dn/nk/chip/histone_occu/sth/test_TF/k4me3.mat.gz"
  title <- readLines(con = mat.file, n = 1)
  df <- readr::read_table2(file = mat.file, comment = "@", col_names = FALSE) %>% 
    as.data.frame()
  mat <- df[, -(1:6)] %>% as.matrix()
  colnames(mat) <- NULL
  gr <- df[, 1:6] %>% `colnames<-`(c("chr", "start", "end", "forth", "fifth", "strand")) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  res <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = mat), rowRanges = gr,
                                                    metadata = list(title = title))
  return(res)
}

deeptools.mat.write <- function(se, assay, out.file) {
  if (grepl(".gz$", out.file)) {
    out.file <- sub(".gz$", "", out.file)
  }
  anno <- rowRanges(se) %>% `names<-`(NULL) %>% as.data.frame()
  if (ncol(anno) > 7) {
    anno <- anno[, 1:7]
  }
  if (ncol(anno) < 7) {
    for (i in (ncol(anno) + 1): 7) {
      ncol[, paste0("ADD", i)] <- "."
    }
  }
  cols <- c("seqnames", "start", "end", colnames(anno)[6:7], "strand")
  anno <- anno[, cols]
  if (is.null(assays(se)[[assay]])) {
    stop("is.null(assays(se)[[assay]])")
  }
  values <- assays(se)[[assay]] %>% as.data.frame()
  df <- cbind(anno, values)
  title <- se@metadata$title
  if (is.null(title)) {
    stop("se@metadata$title is empty")
  }
  dir.create(dirname(out.file), showWarnings = F, recursive = T)
  write(title, out.file)
  write.table(df, out.file, sep = "\t", row.names = F, col.names = F, quote = F,
              append = T)
  system(paste0("gzip -f ", out.file))
  out.file <- paste0(out.file, ".gz")
  return(out.file)
}

deeptools.mat.norm <- function(mat.file, out.file = NULL, norm.to = 1000) {
  if (is.null(out.file)) {
    out.file <- tools::file_path_sans_ext(mat.file) %>% 
      paste0("_norm", utilsFanc::so.formatter(norm.to), ".gz")
  }
  se <- deeptools.mat.read(mat.file)
  mat <- assays(se)$raw
  norm.mat <- diag(norm.to/rowSums(mat)) %*% mat
  assays(se)$norm <- norm.mat
  out.file <- deeptools.mat.write(se = se, assay = "norm", out.file = out.file)
  return(out.file)
}


deeptools.mat.renorm.plot <- function(mat.file, renorm.to = 1000, norm.after.aggr = T,
                                      plot.out) {
  # compared to deeptools.mat.read, this one handles the type of mat.file where there are multiple samples
  title <- readLines(con = mat.file, n = 1)
  df <- readr::read_table2(file = mat.file, comment = "@", 
                           col_names = FALSE, col_types = readr::cols(
                             X2 = readr::col_character(), X3 = readr::col_character())) %>% 
    as.data.frame()
  
  mat <- df[, -(1:6)]
  
  json <- jsonlite::fromJSON(sub("^@", "", title))
  
  df <- lapply(1:length(json$group_labels), function(i) {
    lapply(1:length(json$sample_labels), function(j) {
      i.start <- json$group_boundaries[i] + 1
      i.end <- json$group_boundaries[i+1]
      j.start <- json$sample_boundaries[j] + 1
      j.end <- json$sample_boundaries[j+1]
      print(j.start)
      print(j.end)
      mat <- mat[i.start:i.end,j.start:j.end]
      
      if (norm.after.aggr) {
        aggr <- colSums(mat)
        if (is.null(renorm.to)) {
          renorm.to <- sum(aggr)
        }
        aggr <- renorm.to * (aggr/sum(aggr))
      } else {
        stop("only norm.after.aggr has been developed")
      }
      df <- data.frame(group = json$group_labels[i], sample = json$sample_labels[j],
                       values = aggr, position = 1:length(aggr))
      return(df)
    }) %>% do.call(rbind, .) %>% return()
  }) %>% do.call(rbind, .)
  
  p <- ggplot(df, aes(x = position, y = values)) +
    geom_line(aes(color = sample, group = sample)) + theme_bw() +
    theme(aspect.ratio = 1)
  scFanc::wrap.plots.fanc(list(p), plot.out = plot.out, sub.width = 7)
  invisible(p)
}