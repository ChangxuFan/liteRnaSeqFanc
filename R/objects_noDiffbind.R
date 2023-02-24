peak.merge <- function(peak.files, score.col = 9, th = 8,
                       threads = 1) {
  gr <- utilsFanc::safelapply(peak.files, function(peak.file) {
    peaks <- read.table(peak.file, sep = "\t")
    peaks$V2 <- peaks$V2 + 1
    bPass <- peaks[, score.col] >= th
    peaks <- peaks[bPass, 1:3]
    colnames(peaks) <- c("chr", "start", "end")
    peaks <- makeGRangesFromDataFrame(peaks)
    return(peaks)
  }, threads = threads) %>% do.call(c, .)
  gr <- reduce(gr, ignore.strand = T) %>% sort()
  return(gr)
}

diffbind.2.raw.a2bl <- function(dbo = NULL, mat = NULL,
                                sample.info, bam.col = "bamReads", sample.col = "SampleID", bSingleEnd, 
                                features = NULL, feature.file.col,
                                feature.file.filter.col, feature.file.filter.th,
                                cluster.ident, clusters = NULL,
                                sample.order = NULL,
                                filter.nz = F, filter.size = NULL, filter.samples = NULL,
                                filter.fun = "max", sequential.filter,
                                independentFiltering = T,
                                quantile.norm = F, deseq2.norm.method = "ratio", deseq2.locfunc = NULL,
                                coldata.columns = NULL, coldata.df = NULL,
                                pca.ntop = 10000, pca.groupings = NULL,
                                design.formula, contrast,
                                work.dir,threads = 1, plot.dir = NULL,
                                single.sample = F, debugger.demo = F) {
  # features: you can use 1 vector of loci
  # or a list of loci vectors, each corresponding to a cluster in cluster.ident
  # or you can leave it null and specify feature.file related arguments.
  # this is basically doing diffbind's job: feature.file.col is a column of 
  # sample.info that points to the peak file; feature.file.filter.col is the 
  # column of the peak file that has peak score in it. feature.file.filter.th
  # is the score threshold.
  if (!is.null(dbo)) {
    stop("currently the diffbind route is wrong.")
  }
  if (!is.null(dbo))
    sample.info <- dbo$samples
  if (is.character(sample.info)) {
    sample.info <- read.table(sample.info, header = T)
  }
  if (is.null(plot.dir))
    plot.dir <- paste0(work.dir, "/plots")
  if (is.null(clusters)) {
    clusters <- sample.info[, cluster.ident] %>% unique()
  }
  
  if (!debugger.demo) {
    if (is.null(pca.groupings))
      pca.groupings <- coldata.columns
  }
  
  if (!is.null(sample.info$Replicate)) {
    # back compatibility with DiffBind nomenclature
    if (!grepl("rep", sample.info$Replicate[1]))
      sample.info$Replicate <- paste0("rep", sample.info$Replicate)
  }
  
  if (is.null(mat)) {
    if (!is.null(dbo))
      mat <- diffbind.get.raw.mat(dbo)
    else {
      if (is.null(features)) {
        stop("untested")
        features <- peak.merge(peak.files = sample.info[, feature.file.col],
                               score.col = feature.file.filter.col,
                               th = feature.file.filter.th, threads = threads)
        
      }
      if (is.list(features)) {
        if (!identical(sort(clusters), sort(names(features)))) {
          stop("length(features)!=nrow(sample.info)")
        }
        features.all <- unlist(features) %>% unique()
      } else {
        features.all <- features
      }
      
      mat <- bam.count(sample.info = sample.info, bam.col = bam.col, sample.col = sample.col, 
                       features.gr = features.all, bSingleEnd = bSingleEnd, return.mat = T,
                       sort = T, threads = threads)
      
    }
  }
  
  a2b.list <- sample.info %>% filter(!!as.name(cluster.ident) %in% clusters) %>% 
    split(., f = factor(.[, cluster.ident], levels = unique(.[, cluster.ident]))) %>% 
    utilsFanc::safelapply(function(df) {
      cluster <- df[, cluster.ident][1]
      a2b <- list()
      a2b$root.name <- df[, cluster.ident][1]
      a2b$bulk.mat <- mat[, df[[sample.col]]]
      
      if (!is.null(features)) {
        if (is.list(features)) {
          feature.use <- features[[cluster]]
        } else {
          feature.use <- features
        } 
        a2b$bulk.mat <- a2b$bulk.mat %>% .[rownames(.) %in% feature.use,]
      }
      
      colnames(a2b$bulk.mat) <- colnames(a2b$bulk.mat) %>% sub(paste0(cluster, "_"), "", .)
      a2b$coldata <- df %>% dplyr::mutate(Sample = colnames(a2b$bulk.mat))
      rownames(a2b$coldata) <- colnames(a2b$bulk.mat)
      
      a2b <- s2b.deseq(s2b.obj = a2b, quantile.norm = quantile.norm,
                       norm.method = deseq2.norm.method, locfunc = deseq2.locfunc,
                       filter.nz = filter.nz, filter.size = filter.size, filter.samples = filter.samples,
                       filter.fun = filter.fun, sequential.filter = sequential.filter,
                       independentFiltering = independentFiltering,
                       pca.ntop = pca.ntop, pca.groupings = pca.groupings,
                       design = design.formula, contrast = contrast,
                       sample.order = sample.order, try.hm = F,
                       force.hm = F, force = F, 
                       plot.dir = plot.dir,
                       single.sample = single.sample)
      
      return(a2b)
    }, threads = threads)
  saveRDS(a2b.list, paste0(work.dir, "/bulk.list.Rds"))
  return(a2b.list)
}
