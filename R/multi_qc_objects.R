#
m.qc.general <- function(qc.dir.list) {
  lapply(qc.dir.list, function(qc.dir) {
    samples <- system(paste0("ls ", qc.dir), intern = T)
    dfs <- lapply(samples, function(s) {
      j <- jsonlite::fromJSON(paste0(qc.dir,"/",s, "/qc.json"))
      df <- j$general %>% as.data.frame()
      return(df)
    }) %>% rbind.union.fanc() %>% return()

  }) %>% rbind.union.fanc() %>% return()
}

m.qc.samstat <- function(qc.dir.list) {
  options(digits = 3)
  fields.samstat <- c("total_reads", "pct_mapped_reads", "pct_properly_paired_reads",
                      "with_itself", "pct_singletons")

  my.samstat <- c("sample", "replicate", "n.reads", "pct.dup", "pct.mapped", "pct.properly.paired",
                  "pct.self.mapped", "pct.singleton", "pct.mito" )

  df <- lapply(qc.dir.list, function(qc.dir) {
    samples <- system(paste0("ls ", qc.dir), intern = T)

    dfs <- lapply(samples, function(s) {
      j <- jsonlite::fromJSON(paste0(qc.dir,"/",s, "/qc.json"))
      reps <- names(j$align$samstat)



      raw.stat <- lapply(reps, function (r) {
        stat <- j$align$samstat[[r]] %>% unlist()
        stat <- stat[fields.samstat]
        # modify with_itself to use percentage instead of number of reads.
        stat["with_itself"] <- stat["with_itself"]/stat["total_reads"]

        # insert pct.duplicate after total_reads
        dup <- j$align$dup[[r]]$pct_duplicate_reads
        names(dup) <- "pct.dup"
        stat <- append(stat, dup, after = which(names(stat) == "total_reads"))
        # insert frac.mito at the end.
        mito <- j$align$frac_mito[[r]]$frac_mito_reads
        names(mito) <- "pct.mito"
        stat <- append(stat, mito, after = length(stat))
        stat <- round(stat, digits = 3)
        return(stat)
      })

      filtered.stat <- lapply(reps, function(r) {
        stat <- j$align$nodup_samstat[[r]] %>% unlist()
        stat <- stat[fields.samstat]
        # modify with_itself to use percentage instead of number of reads.
        stat["with_itself"] <- stat["with_itself"]/stat["total_reads"]
        # insert pct.duplicate after total_reads
        dup <- 0
        names(dup) <- "pct.dup"
        stat <- append(stat, dup, after = which(names(stat) == "total_reads"))
        # insert frac.mito at the end.
        mito <- 0
        names(mito) <- "pct.mito"
        stat <- append(stat, mito, after = length(stat))
        stat <- round(stat, digits = 3)
        return(stat)
      })

      cat.stat <- mapply(raw.stat,  filtered.stat, reps, FUN = function(x, y, r) {
        cat.stat <- paste0(y, "/", x)
        libname <- c(s, r)
        #print(libname)
        names(libname) <- c("sample", "replicate")
        cat.stat <- c(libname,cat.stat)
        return(cat.stat)
      }, SIMPLIFY = F)

      df <- Reduce(rbind, cat.stat) %>% as.data.frame()
      colnames(df) <- my.samstat

      return(df)
    }) %>% rbind.union.fanc() %>% return()

  }) %>% rbind.union.fanc()
  rownames(df) <- NULL
  return(df)
}

m.qc.idr <- function(qc.dir.list, rotate=F) {

  df <- lapply(qc.dir.list, function(qc.dir) {
    idr.dir <- sub("qc", "idr", qc.dir)
    #print(idr.dir)
    samples <- system(paste0("ls ", qc.dir), intern = T)

    dfs <- lapply(samples, function(s) {
      j <- jsonlite::fromJSON(paste0(qc.dir,"/",s, "/qc.json"))
      df <- j$replication$reproducibility$idr %>% as.data.frame()
      df <- add.column.fanc(df, data.frame(sample = s), pos=1)

      n.peaks <- j$replication$num_peaks %>% as.data.frame()
      colnames(n.peaks) <- paste0("rep", 1:ncol(n.peaks))
      df <- add.column.fanc(df, n.peaks, after = "sample")

      contrasts <- system(paste0("ls -d ", idr.dir,"/",s,"*"), intern = T) %>%
        sub(".+(rep\\d_vs_rep\\d)", "\\1",.)

      contrasts.df <- lapply(contrasts, function(x) {
        cmd <- paste0("zcat ",  idr.dir,"/",s,"_", x ,"/", "*bfilt.narrowPeak.gz | wc -l ")
        #print(cmd)
        data.frame(V=system(cmd, intern = T)) %>%
          return()
      }) %>% Reduce(cbind, .)
      colnames(contrasts.df) <- contrasts

      df <- add.column.fanc(df, contrasts.df, after = "Np")

      return(df)
    }) %>% rbind.union.fanc() %>% return()

  }) %>% rbind.union.fanc()
  rownames(df) <- df$sample
  df$sample <- NULL
  if (rotate==T)
    df <-rotate.df.fanc(df)
  return(df)
}

m.qc.frag.length <- function(qc.dir.list) {
  options(digits = 3)


  df <- lapply(qc.dir.list, function(qc.dir) {
    samples <- system(paste0("ls ", qc.dir), intern = T)

    dfs <- lapply(samples, function(s) {
      j <- jsonlite::fromJSON(paste0(qc.dir,"/",s, "/qc.json"))
      reps <- names(j$align$frag_len_stat)
      df <- j$align$frag_len_stat %>% Reduce(rbind, .)
      rownames(df) <- NULL
      df <- cbind(data.frame(sample = s, replicate = reps), df)

      return(df)
    }) %>% rbind.union.fanc() %>% return()

  }) %>% rbind.union.fanc()
  rownames(df) <- NULL
  return(df)
}

m.qc.get.file.each.library <- function(dir.list, file.pattern, sample.level = F) {
  options(stringsAsFactors = F)
  lapply(dir.list, function (dir) {
    libs <- system(paste0("ls ", dir), intern = T)
    files <- Sys.glob(paste0(dir, "/", libs, "/", file.pattern))
    if (sample.level == F) {
      samples <- grab.sample(libs, is.short.name = T)
      reps <- grab.rep(libs, is.short.name = T)
      df <- data.frame(sample = samples, rep = reps, file = files)
    } else {
      samples <- libs
      df <- data.frame(sample = samples, file = files)
    }

    return(df)
  }) %>% Reduce(rbind, .) %>% return()
}

m.qc.plot.grid <- function(df, margin.width = 10) {
  # 3 columns: sample; rep; file
  df$sample.n <- factor(df$sample, levels = unique(df$sample)) %>% as.numeric()
  df$rep.n <- sub("rep(\\d)+", "\\1", df$rep) %>% as.numeric()

  samples <- df$sample %>% unique()
  n.samples <- max(df$sample.n)
  n.reps <- max(df$rep.n)
  reps <- paste0("rep", 1:n.reps)

  par(las=1)
  par(mar=c(3,margin.width,1,1))
  plot.new()
  plot(x=c(0,(n.reps)), y = c(0,(n.samples)), type = "n", axes=FALSE, frame.plot=F,
       xlab = "", ylab = "")

  Axis(side=1, labels=reps, at= (1:n.reps) - 0.5)
  Axis(side=2, labels=rev(samples) %>%lapply(smart.break.fanc, margin.width), at = (1:n.samples) - 0.5)
  # add labels:
  # for (i in 1:n.samples) {
  #   text(x = n.reps, y = n.samples - i + 0.5, samples[i])
  # }
  #
  # for (i in 1:n.reps) {
  #   text(x=i-0.5, y= n.samples + 0.1, reps[i])
  # }
  # plot all the images
  for (i in 1:nrow(df)) {
    # df[i,"file"] %>% print()
    img <- readPNG(df[i,"file"])

    rasterImage(img, df[i, "rep.n"] -1, n.samples - df[i, "sample.n"],
                df[i, "rep.n"], n.samples - df[i, "sample.n"] + 1)
  }

}

m.qc.plot.sample <- function(df, n.plots.per.row) {
  n.plots <- nrow(df)
  y.lim <- ceiling(n.plots/n.plots.per.row)
  df$sample.num <- 1:nrow(df)
  df$col.num <- (df$sample.num - floor(df$sample.num/n.plots.per.row) * n.plots.per.row) %>%
    if_else(.==0, n.plots.per.row, .)
  df$row.num <- ceiling(df$sample.num/n.plots.per.row)
  # return(df)
  par(mar=rep(0,4))
  plot.new()
  plot(x=c(0,n.plots.per.row), y = c(0,y.lim), type = "n", frame.plot=F, axes=FALSE,
       xlab = "", ylab = "")

  for (i in 1:nrow(df)) {
    img <- readPNG(df[i, "file"])
    rasterImage(img, df[i, "col.num"] -1, y.lim - df[i, "row.num"],
                df[i, "col.num"], y.lim - df[i, "row.num"] + 1)
  }

  text(x=df$col.num -1  , y = y.lim - df$row.num+0.9,labels =df$sample,  adj = c(0,0), cex = 1)

}

m.qc.json.lib.feature.table <- function(qc.dir.list, json.slice.vector,
                                        feature.subset.pattern = NULL,
                                        sample.subset.pattern = NULL,
                                        replicate.subset.pattern = NULL) {
  df <- lapply(qc.dir.list, function(qc.dir) {
    samples <- system(paste0("ls ", qc.dir), intern = T)

    dfs <- lapply(samples, function(s) {
      j <- jsonlite::fromJSON(paste0(qc.dir,"/",s, "/qc.json"))

      for (i in json.slice.vector) {
        j <- j[[i]]
      }
      # print(j)
      reps <- names(j)
      df <- Reduce(rbind, lapply(j, unlist))
      rownames(df) <- NULL
      df <- cbind(data.frame(sample = s, replicate = reps), df)

      return(df)
    }) %>% rbind.union.fanc() %>% return()

  }) %>% rbind.union.fanc()
  if (!is.null(feature.subset.pattern))
    df <- df[,grep(feature.subset.pattern, colnames(df))]
  if (!is.null(sample.subset.pattern))
    df <- filter(df, grepl(sample.subset.pattern, sample))
  if (!is.null(replicate.subset.pattern))
    df <- filter(df, grepl(replicate.subset.pattern, replicate))

  rownames(df) <- NULL
  colnames(df) <- lapply(colnames(df), smart.break.fanc, 12) %>% unlist()
  return(df)
}

