featureCounts.2.mat <- function(in.tsv, rm.Gm = F) {
  df <- read.table(in.tsv, header = T, sep = "\t")
  rownames(df) <- df$gene
  df$gene <- NULL
  mat <- as.matrix(df)
  if (rm.Gm) {
    mat <- mat[!grepl("^Gm\\d+$", rownames(mat)),]
  }
  return(mat)
}


expr.bar <- function(mat, genes, bSort = T,
                     plot.out = NULL,
                     sub.width = 10, sub.height = 4, n.col = 1,
                     n.show = NULL) {
  pl <- lapply(genes, function(gene) {
    df <- data.frame(sample = colnames(mat), expr = mat[gene,])
    if (bSort)
      df <- df %>% arrange(desc(expr))
    if (is.null(n.show)) {
      n.show <- nrow(df)
    }
    df <- df[1:n.show, ]
    if (bSort)
      df <- df %>% arrange(expr) 
    df <- df %>% dplyr::mutate(sample = factor(sample, levels = sample))
    p <- ggplot(df, aes(x = sample, y = expr)) +
      geom_bar(stat = "identity", alpha = 0.5, fill = "purple4") +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      theme_classic() +
      ggtitle(gene)
    return(p)
  })
  p <- scFanc::wrap.plots.fanc(plot.list = pl, n.col = n.col,
                               sub.width = sub.width, sub.height = sub.height,
                               plot.out = plot.out)
  invisible(p)
}

expr.w.meta <- function(mat, genes, sample.info, x.var, 
                        plot.out = NULL,
                        sub.width = 10, sub.height = 4, n.col = NULL) {
  if (is.character(sample.info)) {
    sample.info <- read.table(sample.info, header = T)
  }
  utilsFanc::check.intersect(x = sample.info$sample, x.name = "sample.info$sample",
                             y = colnames(mat), y.name = "colnames(mat)")
  utilsFanc::check.intersect(genes, "genes", rownames(mat), "rownames(mat)")
  mat <- mat[, sample.info$sample]
  pl <- lapply(genes, function(gene) {
    df <- sample.info
    df$expr <- mat[gene, df$sample]
    library(ggpubr)
    p <- ggbarplot(data = df, x = x.var, y = "expr", add = c("jitter", "mean")) +
      ggtitle(gene)
    return(p)
  })
  
  p <- scFanc::wrap.plots.fanc(plot.list = pl, plot.out = plot.out,n.col = n.col, 
                               sub.width = sub.width, sub.height = sub.height)
  invisible(p)
}

##########
# functions first written for ultralow rna
detection.rate <- function(exp.df, ticks, ticks.log2p1 = T,
                           y.lim = c(0, 20000), x.log2p1 = T,
                           use.bar = F,
                           plot.out = NULL) {
  if (!is.data.frame(exp.df)) {
    exp.df <- as.data.frame(exp.df)
  }
  if ("gene" %in% colnames(exp.df)) {
    rownames(exp.df) <- exp.df$gene
    exp.df$gene <- NULL
  }
  if (is.null(rownames(exp.df))) {
    stop("is.null(rownames(exp.df))")
  }
  if (ticks.log2p1) {
    ticks <- 2^(ticks) - 1
  }
  df <- lapply(ticks, function(tick) {
    n <- sapply(exp.df, function(x) return(sum(x >= tick)))
    return(n)
  }) %>% Reduce(rbind, .) %>% as.data.frame()
  df$cutoff <- ticks
  df.melted <- reshape2::melt(df, id.vars = "cutoff", 
                              variable.name = "sample", 
                              value.name = "n.genes")
  if (x.log2p1) {
    df.melted$cutoff <- log2(df.melted$cutoff + 1)
  }
  p <- ggplot(df.melted, aes(x = cutoff, y = n.genes))
  if (use.bar) {
    p <- p + geom_bar(aes(fill = sample), position = "dodge", stat = "identity")
  } else {
    p <- p + geom_point() +
      geom_line(aes(group = sample, color = sample))
  }
  p <- p + 
    ylim(y.lim)
  if (!is.null(plot.out)) {
    scFanc::wrap.plots.fanc(list(p), plot.out = plot.out, sub.height = 5, sub.width = 8)
  }
  
  return(p)
}

simple.scatter <- function(exp.df, comps.df, plot.out = NULL, transform = c("linear", "log2")) {
  exp.df <- as.data.frame(exp.df)
  pl.scatter <- lapply(transform, function(tr) {
    pl <- comps.df %>% split(., f= 1:nrow(.)) %>% 
      lapply(function(comp) {
        if (tr == "linear") {
          tra <- NULL
        } else if (tr == "log2"){
          tra <- function(x) log2(x + 1)
        }
        exp.df[, comp$x] <- tra(exp.df[, comp$x])
        exp.df[, comp$y] <- tra(exp.df[, comp$y])
        bool.df <- exp.df > 0
        single.rate <- round(1- sum((bool.df[, comp$x] * bool.df[, comp$y]) > 0)/
                               sum((bool.df[, comp$x] + bool.df[, comp$y]) > 0), digits = 2)
        nz <- exp.df[(bool.df[, comp$x] + bool.df[, comp$y]) > 0,]
        corr <- round(cor(nz[, comp$x], nz[, comp$y]), digits = 2)
        p <- scFanc::xy.plot(df = exp.df, x = comp$x, y = comp$y, 
                             title = tr) +
          ggtitle(corr, subtitle = paste0("single detect: ", single.rate))
      })
    return(pl)
  }) %>% Reduce(c, .)
  if (!is.null(plot.out))
    scFanc::wrap.plots.fanc(pl.scatter, plot.out = plot.out)
  invisible(pl.scatter)
}

xy.plot.grid <- function(mat, plot.out = NULL, transformation = NULL, ...) {
  pl <- utilsFanc::safelapply(colnames(mat), function(x) {
    lapply(colnames(mat), function(y) {
      p <- scFanc::xy.plot(df = as.data.frame(mat), x = x, y = y, 
                           transformation = transformation,
                           ...)
      return(p)
    }) %>% return()
  }, threads = 6) %>% Reduce(c, .)
  p <- scFanc::wrap.plots.fanc(plot.list = pl, plot.out = plot.out, n.col = ncol(mat))
  invisible(p)
}

quadrant.mat <- function(exp.df, comp.df, cutoff, bLog2p1 = F,
                         out.file = NULL, plot = T) {
  exp.df <- exp.df %>% as.data.frame()
  if (! "gene" %in% colnames(exp.df)) {
    exp.df$gene <- rownames(exp.df)
  }
  mats <- comp.df %>% split(., f = 1:nrow(.)) %>% 
    lapply(function(comp) {
      x <- comp$x
      y <- comp$y
      df <- exp.df[, c("gene", x, y)]
      if (bLog2p1) {
        df[, x] <- log2(df[, x] + 1)
        df[, y] <- log2(df[, y] + 1)
      }
      
      df$x <- "l"
      df$x[df[, x] > cutoff] <- "h"
      df$y <- "L"
      df$y[df[, y] > cutoff] <- "H"
      df$cat <- paste0(df$x, "_", df$y)
      sum.df <- df %>% group_by(cat) %>% summarise(n = n()) %>% 
        ungroup() %>% as.data.frame()
      rownames(sum.df) <- sum.df$cat
      mat <- outer(c("l", "h"), c("L", "H"), function(x, y) {
        cats <- paste0(x, "_", y)
        
        values <- sum.df[cats, "n"]
        return(values)
      })
      
      colnames(mat) <- paste0(y, c("l", "h"))
      rownames(mat) <- paste0(x, c("l", "h"))
      mat <- (mat/sum(mat)) %>% round(digits = 3)
      return(mat)
    })
  names(mats) <- paste0(comp.df$y, ":", comp.df$x)
  if (!is.null(out.file)) {
    system(paste0("mkdir -p ", dirname(out.file)))
    for (i in 1:length(mats)) {
      if (i == 1) {
        app <- F
      } else {
        app <- T
      }
      write.table(mats[i], out.file, sep = "\t", quote = F, col.names = NA, append = app)
    }
    if (plot) {
      plot.file <- paste0(tools::file_path_sans_ext(out.file), "_hm.png")
      pl <- lapply(mats, function(mat) {
        df <- reshape2::melt(mat)
        df$Var1 <- factor(df$Var1, levels = rownames(mat))
        df$Var2 <- factor(df$Var2, levels = colnames(mat))
        p <- ggplot(df, aes(x = Var1, y = Var2)) +
          geom_tile(aes(fill = value)) +
          scale_fill_gradient(low = "white", high = "red") +
          geom_text(aes(label = value)) +
          xlab("") + ylab("") +
          ggtitle(label = paste0(cutoff, " ", ifelse(bLog2p1, "log2p1", "")))
        return(p)
      })
      scFanc::wrap.plots.fanc(plot.list = pl, plot.out = plot.file)
    }
  }
  return(mats)
}

mean.diff.plot <- function(mats, comp.list, breaks = 0:10, use.diff = F,
                           plot.out = NULL, width = 15) {
  # either a single matrix or a list of matrices. eg. you might have 1 matrix for mouse samples
  # and another for human. When more than one samples are provided, comps should be a list.
  # eg: mats = list(mouse = mat1, human = mat2), 
  # comps = list(mouse = c("rep1:rep2:rep3", "rep4:rep5:rep6"), huamn = c("hrep1:hrep2"))
  flat <- comp.list %>% unlist()
  dup <- flat[duplicated(flat)]
  if (length(dup) > 0) {
    stop(paste0("duplicated comp: ", paste0(dup, collapse = ", ")))
  }
  
  if (!is.list(mats)) {
    mats <- list(miao = mats)
  }
  if (!is.list(comp.list)) {
    comp.list <- list(miao = comp.list)
  }
  
  utilsFanc::check.intersect(names(comp.list), "names(comp.list)", names(mats), "names(mats)")
  
  df <- lapply(names(comp.list), function(name) {
    comps <- comp.list[[name]]
    mat <- mats[[name]]
    df <- lapply(comps, function(comp) {
      samples <- strsplit(comp, ":") %>% unlist()
      utilsFanc::check.intersect(samples, "samples", colnames(mat), name)
      mat <- mat[, samples]
      mat <- mat[rowSums(mat > 0) > 0,]
      Mean <- rowMeans(mat)
      if (ncol(mat) == 2 && use.diff) {
        print("2 samples, using diff as instructed instead of sd")
        d <- abs(mat[, 1] - mat[, 2])/2
      } else {
        d <- rowSds(mat)
      }
      res <- data.frame(gene = rownames(mat), mean = Mean, diff = d, comp = comp)
      return(res)
    }) %>% do.call(rbind, .)
    return(df)
  }) %>% do.call(rbind, .)
  df$mean.log2 <- log2(df$mean + 1)
  df$log2.expr <- cut(df$mean.log2, breaks = breaks)
  df$diff.mean.r <- df$diff/df$mean 
  # df$mean won't be zero. it was filtered by mat <- mat[rowSums(mat > 0) > 0,]
  
  df.n <- df %>% group_by(log2.expr, comp) %>% summarise(n = n()) %>% 
    ungroup() %>% group_by(comp) %>% mutate(n = n/sum(n)) %>% 
    as.data.frame()
  
  p <- ggplot(df, aes(x = log2.expr, y = diff.mean.r, color = comp)) + 
    geom_boxplot(outlier.size = 0.01) +
    geom_point(data = df.n, aes(x = log2.expr, y = n, color = comp),
                position = position_dodge(width = 0.8))
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    ggsave(plot.out, p, width = width, height = 3, 
           dpi = 100)
  }
  invisible(p)
}


mat.plot.density <- function(mat, bLog2p1 = F, group.df = NULL, plot.out = NULL,
                             text.size = 8) {
  # group.df: must have columns: sample, group. other columns ignored.
  if (bLog2p1) {
    mat <- log2(mat + 1)
  }
  df <- mat %>% reshape2::melt(varnames = c("gene", "sample"))
  if (is.null(group.df)) {
    df$sample <- factor(df$sample, levels = colnames(mat))
    p <- ggplot(df, aes(x = value, color = sample))
  } else {
    utilsFanc::check.intersect(c("sample", "group"), "required columns",
      colnames(group.df), "colnames(group.df)")
    dups <- group.df$sample %>% .[duplicated(.)]
    if (length(dups) > 0) {
      stop("length(dups) > 0")
    }
    df <- left_join(df, group.df)
    df$sample <- factor(df$sample, levels = colnames(mat))
    p <- ggplot(df, aes(x = value, color = sample, linetype = group))
  }
  p <- p + geom_density() + 
    theme(aspect.ratio = 1) +
    theme(text = element_text(size = text.size)) +
    theme(legend.key.size = unit(0.1, "in"))
    # guides(color = guide_legend(override.aes = list(size = 1)))
  if (!is.null(plot.out)) {
    scFanc::wrap.plots.fanc(list(p), plot.out = plot.out, sub.width = 6)
  }
  invisible(p)
}

mat.subset.to.mim <- function(mats) {
  if (!is.list(mats)) {
    mats <- list(mats)
    flag <- T
  } else {
    flag <- F
  }
  min.size <- min(sapply(mats, function(mat) return(min(colSums(imm$raw.mat)))))
  sub.mats <- lapply(mats, function(mat) {
    genes <- rownames(mat)
    sub.mat <- apply(mat, 2, function(x) {
      return(rbinom(length(x), x, min.size/sum(x)))
    })
    rownames(sub.mat) <- genes
    return(sub.mat)
  })
  if (flag) {
    sub.mats <- sub.mats[[1]]
  }
  return(sub.mats)
}

# END: functions written for ultralow rna
###########