trust.count <- function(sample.info, work.dir, levels = c("gene", "subgroup", "allele" ),
                        single.cell.mode = F,
                        use.consensus_count = T, only.use.groups = NULL,
                        skip.deseq2 = F, TCR.mode = F, na.pct.threshold = 10,
                        # deseq2 arguments (mandatory)
                        design.formula, contrast, pca.groupings,
                        # deseq2 arguments (optional)
                        quantile.norm = F, deseq2.norm.method = "ratio", deseq2.locfunc = NULL,
                        filter.nz = F, filter.size = NULL, filter.samples = NULL,
                        filter.fun = "max", sequential.filter,
                        independentFiltering = F,
                        pca.ntop = 10000, 
                        sample.order = NULL, single.sample = F) {
  # # by default, we separate by heavy vs light chains
  # if (is.null(plot.dir))
  #   plot.dir <- paste0(work.dir, "/plots")
  # 
  # if (is.character(sample.info)) {
  #   sample.info <- read.table(sample.info, header = T)
  # }
  # fields <- c("sample", "airr")
  # utilsFanc::check.intersect(x = fields, "required fields", 
  #                            colnames(sample.info), "colnames(sample.info)")
  # airr <- sample.info[, fields] %>% unique() %>% 
  #   split(., f = 1:nrow(.)) %>% lapply(function(s) {
  #     airr <- read.table(s$airr, header = T, sep = "\t") 
  #     airr$sample <- s$sample
  #     return(airr)
  #   }) %>% do.call(rbind, .)
  # 
  # fields <- c("sequence_id", "sample", 
  #             "v_call", "d_call", "j_call", "c_call", 
  #             "consensus_count")
  # utilsFanc::check.intersect(
  #   fields, "mandatory columns",
  #   colnames(airr), "colnames(airr)")
  # 
  # airr <- airr[, fields]
  # airr$group <- NA
  # call.cat <- paste0(airr$v_call, airr$d_call, airr$j_call, airr$c_call)
  # if (TCR.mode) {
  #   print("Counting T cell receptors")
  #   stop("TCR mode not developed yet")
  #   
  # } else {
  #   print("Counting B cell receptors")
  #   airr$group[grepl("IGL|IGK", call.cat)] <- "IGLK"
  #   airr$group[grepl("IGH", call.cat)] <- "IGH"
  # }
  # 
  # 
  # total.contigs <- nrow(airr)
  # na.contigs <- sum(is.na(airr$group))
  # na.pct <- round(100*na.contigs/total.contigs, digits = 2)
  # na.thresh <- 1
  # if (na.pct > na.thresh) {
  #   stop(paste0("Too many NA's in airr$group. ", na.pct, "% of the contigs are NA. ", 
  #               "This is higher than the threshold (", na.thresh, "%) set by Changxu Fan", 
  #               "Maybe you are looking at TCRs somehow?"))
  # }
  # airr <- airr[!is.na(airr$group),]
  
  airr <- trust.airr.read(sample.info = sample.info, TCR.mode = TCR.mode, 
                          na.pct.threshold = na.pct.threshold)
  
  if (!is.null(only.use.groups)) {
    airr <- airr %>% filter(group %in% only.use.groups)
    if (nrow(airr) < 1) 
      stop("nrow(airr) < 1")
  }
  
  s2b.list <- lapply(levels, function(level) {
    s2b.list <- airr %>% split(., f = factor(.$group, levels = unique(.$group))) %>% 
      lapply(function(df) {
        
        group <- df$group[1]
        
        if (level == "gene") {
          segments <- c("v", "j", "c")
        } else if (level ==  "subgroup") {
          segments <- c("v")
        } else if (level == "allele") {
          segments <- c("v", "j")
        } else {
          stop(paste0("level ", level, " is not recognized. "))
        }
        if (group %in% c("IGH", "TRB")) {
          segments <- c(segments, "d") %>% utilsFanc::sort.by(c("v", "d", "j", "c"))
        }
        
        s2b.list <- lapply(segments, function(segment) {
          call <- paste0(segment, "_call")
          if (level == "allele") {
            abbr <- "al"
          } else if (level == "gene") {
            abbr <- "gn"
            df[, call] <- sub("\\*.+$", "", df[, call])
          } else if (level == "subgroup") {
            abbr <- "sg"
            df[, call] <- gsub("\\-.+$", "", df[, call]) %>% gsub("\\*.+$", "", .)
          } else {
            stop(paste0("level ", level, " is not recognized. "))
          }
          
          # if (single.cell.mode) {
          #   # utilsFanc::check.intersect("cell_id", "cell_id", colnames(df), "colnames(df)")
          #   # df$cell_id <- paste0(df$sample, "#", df$cell_id, "-1")
          #   # utilsFanc::check.dups(df$cell_id, "df$cell_id")
          #   # df <- df[df[, call] != "",]
          #   # s2b <- list(
          #   #   root.name = paste0(abbr, "_", group, "_", segment),
          #   #   sc.mat = df[, c(call, "cell_id", "consensus_count")] %>% 
          #   #     reshape2::acast(as.formula(paste0(call, " ~ cell_id")), value.var = "consensus_count"),
          #   #   airr = df
          #   # )
          #   # s2b$sc.mat[is.na(s2b$sc.mat)]  <- 0
          #   # return(s2b)
          #   
          # } else {
          # }
          df <- df %>% dplyr::group_by(!!as.name(call), sample)
          
          if (single.cell.mode) {
            utilsFanc::check.intersect("cell_id", "cell_id", colnames(df), "colnames(df)")
            df$cell_id <- paste0(df$sample, "#", df$cell_id, "-1")
            utilsFanc::check.dups(df$cell_id, "df$cell_id")
            df <- suppressMessages(dplyr::summarise(df, n = n()))
            counting.flag <- "sc"
            
          } else {
            if (use.consensus_count) {
              df <- suppressMessages(dplyr::summarise(df, n = sum(consensus_count)))
              counting.flag <- "nR"
            } else {
              df <- suppressMessages(dplyr::summarise(df, n = n()))
              counting.flag <- "nC"
            }
          }
          df <- df %>% dplyr::ungroup() %>% as.data.frame()
          df <- df[!df[, call] %in% c(""),]
          
          df <- utilsFanc::change.name.fanc(df = df, cols.from = call, cols.to = "gene")
          s2b <- list()
          s2b$root.name <- paste0(counting.flag, "_", abbr,
                                  "_", group, "_", segment)
          s2b$bulk.mat <- df %>% reshape2::acast(gene ~ sample, value.var = "n")
          s2b$bulk.mat[is.na(s2b$bulk.mat)] <- 0
          s2b$coldata <- sample.info[, !colnames(sample.info) %in% c("airr")] %>% 
            dplyr::filter(sample %in% colnames(s2b$bulk.mat))
          rownames(s2b$coldata) <- s2b$coldata$sample
          s2b$bulk.mat <- s2b$bulk.mat[, s2b$coldata$sample]
          if (!skip.deseq2) {
            plot.dir <- paste0(work.dir, "/deseq2_assess/")
            s2b <- s2b.deseq(
              s2b.obj = s2b, quantile.norm = quantile.norm,
              norm.method = deseq2.norm.method, locfunc = deseq2.locfunc,
              filter.nz = filter.nz, filter.size = filter.size, filter.samples = filter.samples,
              filter.fun = filter.fun, sequential.filter = sequential.filter,
              independentFiltering = independentFiltering,
              pca.ntop = pca.ntop, pca.groupings = pca.groupings,
              design = design.formula, contrast = contrast,
              sample.order = sample.order, try.hm = F,
              force.hm = F, force = F, 
              plot.dir = plot.dir,
              single.sample = single.sample
            )
          } else {
            genes <- rownames(s2b$bulk.mat )
            s2b$bulkNorm <- s2b$bulk.mat %>% as.data.frame() %>%
              lapply(function(x) return(x/sum(x))) %>% as.data.frame()
            s2b$bulkNorm <- cbind(data.frame(gene = genes), s2b$bulkNorm)
            colnames(s2b$bulkNorm) <- sub("^X", "", colnames(s2b$bulkNorm))
            s2b$bulkNorm.mat <- s2b$bulkNorm[, -1] %>% as.matrix()
            rownames(s2b$bulkNorm.mat) <- genes
          }
          
          return(s2b)
        })
        return(s2b.list)
      }) %>% Reduce(c, .)
    return(s2b.list)
  }) %>% Reduce(c, .)
  
  names(s2b.list) <- lapply(s2b.list, function(s2b) return(s2b$root.name)) %>% unlist()
  dir.create(work.dir, recursive = T, showWarnings = F)
  saveRDS(s2b.list, paste0(work.dir, "/bulk.list.Rds"))
  return(s2b.list)
}


trust.airr.read <- function(sample.info, TCR.mode = F, na.pct.threshold = 10, file.out = NULL) {
  if (is.character(sample.info)) {
    sample.info <- read.table(sample.info, header = T)
  }
  fields <- c("sample", "airr")
  utilsFanc::check.intersect(x = fields, "required fields", 
                             colnames(sample.info), "colnames(sample.info)")
  airr <- sample.info[, fields] %>% unique() %>% 
    split(., f = 1:nrow(.)) %>% lapply(function(s) {
      airr <- read.table(s$airr, header = T, sep = "\t") 
      airr$sample <- s$sample
      return(airr)
    }) %>% do.call(rbind, .)
  fields <- c("sequence_id", "sample", 
              "v_call", "d_call", "j_call", "c_call", 
              "consensus_count")
  utilsFanc::check.intersect(
    fields, "mandatory columns",
    colnames(airr), "colnames(airr)")
  
  # airr <- airr[, fields]
  airr$group <- NA
  call.cat <- paste0(airr$v_call, airr$d_call, airr$j_call, airr$c_call)
  
  if (TCR.mode) {
    print("Counting T cell receptors")
    airr$group[grepl("TRB", call.cat)] <- "TRB"
    airr$group[grepl("TRA", call.cat)] <- "TRA"
  } else {
    print("Counting B cell receptors")
    airr$group[grepl("IGL|IGK", call.cat)] <- "IGLK"
    airr$group[grepl("IGH", call.cat)] <- "IGH"
  }
  
  total.contigs <- nrow(airr)
  na.contigs <- sum(is.na(airr$group))
  na.pct <- round(100*na.contigs/total.contigs, digits = 2)
  na.thresh <- na.pct.threshold
  if (na.pct > na.thresh) {
    stop(paste0("Too many NA's in airr$group. ", na.pct, "% of the contigs are NA. ", 
                "This is higher than the threshold (", na.thresh, "%) set by Changxu Fan", 
                "Maybe you are looking at TCRs somehow?"))
  }
  airr <- airr[!is.na(airr$group),]
  if (!is.null(file.out)) {
    dir.create(dirname(file.out), showWarnings = F, recursive = T)
    write.table(airr, file.out, sep = "\t", quote = F, col.names = T, row.names = F)
  }
  invisible(airr)
}

trust.junction.sum <- function(airr.cat, use.consensus_count = F, out.dir, root.name = NULL) {
  # airr.cat: catted airr files. generated using trust.airr.read()
  airr <- airr.cat
  rm(airr.cat)
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  
  fields <- c("junction_aa", "sample", "consensus_count")
  utilsFanc::check.intersect(fields, " required columns",
                             colnames(airr), "colnames(airr.cat)")
  airr <- airr[, fields]
  
  if (use.consensus_count) {
    airr <- airr[rep(seq_along(airr$consensus_count), airr$consensus_count), ]
  }
  
  airr <- airr %>% dplyr::filter(junction_aa != "")
  sum.df <- airr %>% dplyr::group_by(sample) %>% 
    dplyr::summarise(n_distinct = length(unique(junction_aa)),
                     n_total = n()) %>% 
    dplyr::ungroup() %>% as.data.frame() %>% 
    dplyr::mutate(ratio = round(n_distinct/n_total, digits = 3))
  
  dir.create(out.dir, showWarnings = F, recursive = T)
  write.table(sum.df, paste0(out.dir, "/", root.name, "_junctionAA_summary.tsv"),
              quote = F, col.names = T, row.names = F, sep = "\t")
  
  airr$length <- nchar(airr$junction_aa)
  
  p <- ggplot(airr, aes(x = length, y = sample)) +
    ggridges::geom_density_ridges(alpha = 0.3, show.legend = F, 
                                  aes(color = sample, fill = sample)) +
    theme(aspect.ratio = 1) +
    scale_x_continuous(breaks = scales::breaks_pretty())
  
  scFanc::wrap.plots.fanc(list(p), plot.out = paste0(
    out.dir, "/", root.name, "_junctionAA_length_distro.png"))
  
  return(sum.df)
}

trust.dist.aa.1.core <- function(aa.vec, length.mercy = 2) {
  # the simpliest way to measure how similar 2 amino acid sequences are:
  # how many amino acids are shared between the 2
  # length.mercy: ignore (gives an NA) pairs of aa sequences with dramatic differences in length. 
  if (any(aa.vec == "")) {
    warning("some elements in aa.vec are empty")
  }
  aa.vec <- aa.vec[aa.vec != ""]
  aa <- strsplit(aa.vec, "")
  aal <- sapply(aa, function(x) return(length(unique(x))))
  df <- data.frame(seq = rep(1:length(aa), sapply(aa, length)), 
                   aa = unlist(aa))
  df <- unique(df)
  df$value <- 1
  mat <- reshape2::acast(df, formula = seq ~ aa, value.var = "value")
  mat[is.na(mat)] <- 0
  
  n.shared <- mat %*% t(mat)
  # t <- unique(diag(n.shared) - sapply(aa, function(x) length(unique(x))))
  # t is zero
  
  length.sum <- outer(aal, aal, "+")
  dist <- 1 - n.shared/length.sum
  dist <- dist %>% scales::rescale(to=c(0,1))
  dist <- round(dist, digits = 4)
  # mask out those with high length differences.
  
  length.diff <- outer(aal, aal, "-")
  dist[abs(length.diff) > length.mercy] <- NA
  return(dist)
}

trust.dist.aa.1 <- function(airr.cat, length.mercy = 2, n.subsample = 1000,
                            out.dir, root.name = NULL) {
  airr <- airr.cat
  rm(airr.cat)
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  fields <- c("junction_aa", "sample")
  utilsFanc::check.intersect(fields, " required columns",
                             colnames(airr), "colnames(airr.cat)")
  airr <- airr[, fields]
  airr <- airr[airr$junction_aa != "",]
  diff.df <- airr %>% split(., f = factor(.$sample, levels = unique(.$sample))) %>% 
    lapply(function(airr) {
      print(paste("processing sample: ", airr$sample[1]))
      dist <- trust.dist.aa.1.core(aa.vec = airr$junction_aa, length.mercy = length.mercy)
      diag(dist) <- NA
      dist <- as.vector(dist) 
      dist <- dist[!is.na(dist)]
      dist <- dist[sample(1:length(dist), n.subsample, replace = T)]
      res <- data.frame(sample = airr$sample[1], dist = dist)
      return(res)
    }) %>% do.call(rbind, .)
  p <- ggplot(diff.df, aes(x = dist, y = sample, color = sample, fill = sample)) +
    ggridges::geom_density_ridges(alpha = 0.3) +
    theme(aspect.ratio = 1)
  dir.create(out.dir, showWarnings = F, recursive = T)
  scFanc::wrap.plots.fanc(plot.list = list(p), 
                          plot.out = paste0(out.dir, "/", root.name, "_junctionAA_dist1_distro.png"))
  invisible(p)
}

trust.kmer.umap <- function(sample.info, TCR.mode = F, k = 3, residues = "AA", seed = 42, 
                            junction.length.min = NULL, junction.length.max = NULL,
                            out.dir, root.name = NULL) {
  if (is.null(root.name)) root.name <- basename(out.dir)
  file.name <- paste0(root.name, "_", residues,k, "mer")
  dir.create(out.dir, showWarnings = F, recursive = T)
  
  utilsFanc::t.stat("Reading in airr files")
  airr <- trust.airr.read(sample.info = sample.info, TCR.mode = TCR.mode)
  airr <- airr[airr$junction_aa != "",]
  if (!is.null(junction.length.min))
    airr <- airr[nchar(airr$junction_aa) >= junction.length.min,]
  if (!is.null(junction.length.max))
    airr <- airr[nchar(airr$junction_aa) <= junction.length.max,]
  
  rownames(airr) <- paste0(airr$sample, ".", airr$sequence_id)
  airr <- airr[, c("junction_aa", "sample")] %>% unique()
  airr$name <- rownames(airr)
  aa <- ape::as.AAbin(strsplit(airr$junction_aa, ""))
  
  utilsFanc::t.stat("Generating Kmer matrix")
  kmat <- kmer::kcount(aa, k = k, residues = residues, compress = T)
  
  kmat <- kmat[, colSums(kmat) > 0]
  rownames(kmat) <- airr$name
  saveRDS(kmat, paste0(out.dir, "/", file.name, "_kmat.Rds"))
  
  utilsFanc::t.stat("Calculating UMAP")
  um <- uwot::umap(kmat)
  colnames(um) <- c("UMAP1", "UMAP2")
  um <- as.data.frame(um)
  um$name <- airr$name
  um$sample <- airr$sample
  saveRDS(um, paste0(out.dir, "/", root.name, "_", residues,k, "mer.Rds"))
  
  utilsFanc::t.stat("Generating Seurat object")
  library(Seurat)
  so <- scFanc::fake.so.gen.2(meta.data = airr)
  rownames(um) <- um$name
  so <- seurat.add.embed(so = so, embed.df = um, embedding.name = "UMAP")
  saveRDS(so, paste0(out.dir, "/", file.name, "_so.Rds"))
  
  utilsFanc::t.stat("Plotting")
  pl <- list()
  pl$all <- xy.plot(um, x = "UMAP1", y = "UMAP2", add.abline = F) + ggtitle("all")
  xlim <- c(min(um$UMAP1), max(um$UMAP1))
  ylim <- c(min(um$UMAP2), max(um$UMAP2))
  pl.each <- um %>% split(., f = factor(.$sample, levels = unique(.$sample))) %>% 
    lapply(function(um) {
      p <- xy.plot(um, x = "UMAP1", y = "UMAP2", x.limit = xlim, y.limit = ylim, add.abline = F) +
        ggtitle(um$sample[1])
      return(p)
    })
  pl <- c(pl, pl.each)
  
  p <- scFanc::wrap.plots.fanc(pl, plot.out = paste0(out.dir, "/", root.name, "_", residues,k, "mer.png"))
  invisible(p)
  
}

trust.umap.dist.assess <- function(kmat, umap.df, n = 1000, dist.method = "euclidean",
                                   seed = 42, seed.scramble = 0, 
                                   out.dir, root.name = NULL) {
  if (is.null(root.name)) root.name <- basename(root.name)
  um <- umap.df
  rm(umap.df)
  set.seed(seed = seed)
  cells <- sample(1:nrow(um), size = n, replace = F) %>% sort()
  um <- um[cells, c("UMAP1", "UMAP2")] %>% as.matrix()
  kmat <- kmat[cells, ]
  
  dist.u <- dist(um, method = dist.method)
  dist.k <- dist(kmat, method = dist.method)
  
  set.seed(seed = seed)
  points <- sample(1:(n*n), size = n, replace = F) %>% sort()
  
  set.seed(seed = seed.scramble)
  points.scr <- sample(1:(n*n), size = n, replace = F) %>% sort()
  
  ddf <- data.frame(
    k = as.vector(as.matrix(dist.k))[points],
    u = as.vector(as.matrix(dist.u))[points],
    scr = as.vector(as.matrix(dist.u))[points.scr]
  )
  ddf <- ddf %>% dplyr::filter(k != 0, u != 0, scr != 0)
  pl <- lapply(c("u", "scr"), function(x) {
    p <- xy.plot(ddf, x = "k", y = x, add.abline = F, add.corr = T)
  })
  
  p <- scFanc::wrap.plots.fanc(pl, plot.out = paste0(
    out.dir, "/", root.name, "_n", n, "_s", seed, "_scr", seed.scramble, ".png"))
  invisible(p)
}

trust.cluster.assess <- function(so, group.by = "seurat_clusters", groups = NULL,
                                 split.by = NULL, splits = NULL,
                                 subset.n = 20, 
                                 out.dir, root.name = NULL ) {
  if (is.null(root.name)) root.name <- basename(out.dir)
  
  x <- scFanc::get.cell.list(obj = so, is.ao = F, style = "Seurat",
                             n.cells.each = subset.n, 
                             split.by = split.by, splits = splits,
                             group.by = group.by, groups = groups)
  seq.df <- lapply(names(x), function(y) {
    x <- x[[y]]
    seqs <- so@meta.data[x, "junction_aa"]
    df <- data.frame(group = y, name = x, seq = seqs)
    return(df) 
  }) %>% do.call(rbind, .)
  
  # length distro
  seq.df$len <- nchar(seq.df$seq)
  p.len <- ggplot(seq.df, aes(x = group, y = len)) +
    geom_violin() +
    geom_jitter()
  
  scFanc::wrap.plots.fanc(list(p.len), plot.out = paste0(out.dir, "/", root.name, "_length_violin.png"))
  
  fa <- paste0(out.dir, "/", root.name, ".fa")
  aligned.fa <- sub(".fa$", "_aligned.fa", fa)
  seqinr::write.fasta(as.list(seq.df$seq), as.string = T, file.out = fa,
                      names = paste0(seq.df$group, "..", seq.df$name))
  abaFanc2::mafft.fanc(fa, aligned.fa = aligned.fa)
  return()
}

trust4.clone.sharing <- function(sample.info, out.dir, root.name = NULL, seed = 42) {
  if (is.null(root.name)) root.name <- basename(out.dir)
  airr <- trust.airr.read(sample.info = sample.info)
  airr <- airr[, c("sample", "junction_aa")] %>% 
    dplyr::filter(junction_aa != "") %>% 
    unique()
  
  samples <- unique(airr$sample)
  airr.list <- airr %>% split(., f = factor(.$sample, levels = samples))
  
  df <- lapply(samples, function(x) {
    aax <- airr.list[[x]]$junction_aa
    nx <- length(aax)
    df <- lapply(samples, function(y) {
      aay <- airr.list[[y]]$junction_aa
      ny <- length(aay)
      n <- min(nx, ny)
      if (nx > n) {
        set.seed(seed = seed)
        aax <- aax[sort(sample(1:nx, size = n, replace = F))]
      } else {
        set.seed(seed = seed)
        aay <- aay[sort(sample(1:ny, size = n, replace = F))]
      }
      stats <- data.frame(x = x, y = y)
      stats$nshare <- length(intersect(aax, aay))
      stats$frac <- round(stats$nshare/n, digits = 2)
      return(stats)
    }) %>% do.call(rbind, .)
    return(df)
  }) %>% do.call(rbind, .)
  mat <- reshape2::acast(df, formula = x ~ y, value.var = "frac")
  out.file <- paste0(out.dir, "/", root.name, "_cloneShare.tsv")
  dir.create(out.dir, showWarnings = F, recursive = T)
  wdf <- cbind(data.frame(sample = rownames(mat)), as.data.frame(mat))
  write.table(wdf, file = out.file, quote = F, sep = "\t", col.names = T, row.names = F)
  return(mat)
}

trust4.dist.distro <- function(sample.info, seq.type = "junction_aa", 
                               group.by = "sample", n.each = 1000, seed = 42,
                               k = 3, dist.method = "euclidean",
                               out.dir, root.name = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  dir.create(out.dir, showWarnings = F, recursive = T)
  if (seq.type == "junction_aa") {
    kmer.type <- "AA"
  } else if (seq.type == "junction") {
    kmer.type <- "DNA"
  } else {
    stop("seq.type must be 'junction_aa' or 'junction'")
  }
  
  if (is.character(sample.info)) {
    sample.info <- read.table(sample.info, header = T, sep = "\t")
  }
  if (!group.by %in% colnames(sample.info)) {
    stop(paste0("group.by not found in sample.info"))
  }
  
  airr <- trust.airr.read(sample.info = sample.info)
  airr <- dplyr::left_join(airr, sample.info, by = "sample")
  rownames(airr) <- paste0(airr$sample, ".", airr$sequence_id)
  airr <- airr[, c("sample", seq.type, group.by)] %>% unique() # "name" shouldn't 
  # be in here, otherwise unique() wouldn't work
  airr$name <- rownames(airr)
  airr <- airr[airr[, seq.type] != "",]
  groups <- airr[, group.by] %>% unique()
  if (!is.null(n.each)) {
    n.dist <- 5 * n.each
    airr <- airr %>% split(., f = factor(.[, group.by], levels = groups)) %>% 
      lapply(function(airr) {
        if (nrow(airr) > n.each) {
          set.seed(seed = seed)
          airr <- airr[sort(sample(1:nrow(airr), size = n.each, replace = F)),]
        }
        return(airr)
      }) %>% do.call(rbind, .)
  }
  
  utilsFanc::t.stat("Generating Kmer matrix", " based on ", nrow(airr), " sequences")
  if (seq.type == "junction_aa") {
    seqbin <- ape::as.AAbin(strsplit(airr[, seq.type], ""))
  } else if (seq.type == "junction") {
    seqbin <- ape::as.DNAbin(strsplit(airr[, seq.type], ""))
  }
  
  kmat <- kmer::kcount(seqbin, k = k, residues = kmer.type, compress = T)
  kmat <- kmat[, colSums(kmat) > 0]
  rownames(kmat) <- airr$name
  # saveRDS(kmat, "~/test/sth/kmat.Rds")
  # saveRDS(seqbin, "~/test/sth/seqbin.Rds")
  if (dist.method == "euclidean") {
    dmat <- distances::distances(kmat)
  }
  seeds.mat <- matrix(3 * (1:(length(groups)^2)), length(groups))
  colnames(seeds.mat) <- rownames(seeds.mat) <- groups
  utilsFanc::t.stat(paste0("Calculating distances. "))
  dist.df <- lapply(groups, function(x) {
    lapply(groups, function(y) {
      s <- lapply(1:2, function(j) {
        i <- ifelse(j == 1, x, y)
        z <- which(airr[, group.by] == i)
        if (!is.null(n.each)) {
          set.seed(seed = seed + seeds.mat[x, y] + j)
          z <- z[sample(1:length(z), size = n.dist, replace = T)]
          # note you shouldn't sort!
        }
        return(z)
      }) 
      names(s) <- c("x", "y")

      if (dist.method == "euclidean") {
        dist <- sapply(1:length(s$x), function(i) {
          return(dmat[s$x[i], s$y[i]])
        })
        dist <- data.frame(xid = airr$name[s$x], yid = airr$name[s$y], dist = dist)
      }
      # dmat <- dist(kmat[unlist(s), ], method = dist.method) %>% as.matrix()
      
      
      if (any(is.na(dist$dist))) {
        stop("any(is.na(dist$dist))")
      }
      
      dist$x <- x
      dist$y <- y
      return(dist)
    }) %>% do.call(rbind, .) %>% return()
  }) %>% do.call(rbind, .)
  data <- list(airr = airr, kmat = kmat, dmat = dmat, dist.df = dist.df,
               args = list(
                 seq.type = seq.type, 
                 group.by = group.by, n.each = n.each, seed = seed,
                 k = k, dist.method = dist.method,
                 out.dir = out.dir, root.name = root.name))
  saveRDS(data, paste0(out.dir, "/", root.name, "_data.Rds"))
  p <- ggplot(data$dist.df, aes(x = dist)) +
    geom_density() +
    facet_grid(x ~ y) +
    theme(aspect.ratio = 1) +
    ggsave(paste0(out.dir, "/", root.name, ".png"),
           width = 2 * length(groups), height = 2 * length(groups), dpi = 300)
  invisible(p)
}