nanopore.filter.bam.by.decision <- function(
  bam, out.dir = NULL, decision.df, decisions.use = NULL, threads = 1) {
  # decision.df: mandatory columns: read_id; decision
  if (is.null(out.dir)) {
    out.dir <- dirname(bam)
  }
  cols.required <- c("read_id", "decision")
  utilsFanc::check.intersect(cols.required, "required columns", 
                             colnames(decision.df), "colnames(decision.df)")
  if (!is.null(decisions.use)) {
    decision.df <- decision.df %>% dplyr::filter(decision %in% decisions.use)
  }
  if (nrow(decision.df) < 1) {
    stop("nrow(decision.df) < 1")
  }
  root <- basename(bam) %>% tools::file_path_sans_ext()
  out.bams <- decision.df %>% split(., f = .$decision) %>% 
    utilsFanc::safelapply(function(df) {
      decision <- df$decision[1]
      print(paste0("processing: ", decision))
      out.bam <- paste0(out.dir, "/", root, "_", decision, ".bam")
      bam.filter.by.qname(bam = bam, out.bam = out.bam, read.names = df$read_id)
      return(out.bam)
    })
  return(out.bams)
}

nanopore.on.target.rate <- function(seqsum, as.final, bam, as.bed, plot.df.tsv = NULL,
                                    plot.out) {
  # seqsum: the sequencing summary file
  # as.final: you generate this from adaptive sampling summary csv. It's the final decision
  # of each read.
  # as.bed: the target region that you use during adaptive sampling.
  if (is.null(plot.df.tsv)) {
    if (is.character(seqsum)) {
      seqsum <- read.table(seqsum, header = T, sep = "\t")
    }
    if (is.character(as.final)) {
      as.final <- read.table(as.final, sep = "\t", header = T)
    }
    nr <- nrow(seqsum)
    seqsum <- seqsum[, "read_id", drop = F]
    as.final <- as.final[, c("read_id","decision"), drop = F]
    seqsum <- dplyr::left_join(seqsum, as.final, by = "read_id")
    if (nrow(seqsum) != nr) {
      stop("check point 1: nrow(seqsum) != nr")
    }
    seqsum$decision <- as.character(seqsum$decision)
    seqsum$decision[is.na(seqsum$decision)] <- "missing"
    as.gr <- utilsFanc::import.bed.fanc(as.bed, return.gr = T)
    # bam.gr <- GenomicAlignments::readGAlignments(bam, use.names = T)
    o <- findOverlaps(bam.gr, as.gr)
    bam.df <- data.frame(read_id = names(bam.gr), on.target = "off")
    bam.df$on.target[unique(queryHits(o))] <- "on"
    
    bam.df <- bam.df %>% dplyr::arrange(desc(on.target))
    bam.df <- bam.df[!duplicated(bam.df$read_id), ]
    # if one of the alignments of a read is on target, we consider this read on-target.
    
    seqsum <- dplyr::left_join(seqsum, bam.df, by = "read_id")
    seqsum$on.target[is.na(seqsum$on.target)] <- "no_aln"
    if (nrow(seqsum) != nr) {
      stop("check point 2: nrow(seqsum) != nr")
    }
    root <- tools::file_path_sans_ext(plot.out)
    dir.create(dirname(root), showWarnings = F, recursive = T)
    write.table(seqsum, paste0(root, "_", "seqsum", ".tsv"), sep = "\t",
                row.names = F, col.names = T, quote = F)
  } else{
    seqsum <- read.table(plot.df.tsv, header = T, sep = "\t")
  }
  seqsum <- seqsum %>% dplyr::group_by(decision, on.target) %>% 
    dplyr::summarise(n = n()) %>% dplyr::ungroup() %>% 
    dplyr::group_by(decision) %>% 
    dplyr::summarise(on.target = on.target, n = n, frac = n/sum(n)) %>% 
    dplyr::ungroup() %>% as.data.frame()
  
  root <- tools::file_path_sans_ext(plot.out)
  dir.create(dirname(root), showWarnings = F,recursive = T)
  write.table(seqsum, file = paste0(root, "_frac.tsv"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  seqsum$decision <- factor(
    seqsum$decision, 
    levels = c("missing","no_decision", "unblock", "unblock_hit_outside_bed", "stop_receiving"))
  seqsum$on.target <- factor(seqsum$on.target, levels = c("on", "off", "no_aln"))
  # seqsum$frac.label <- paste0(round(seqsum$frac, digits = 3), "\n", seqsum$n)
  p.frac <- ggplot(seqsum, aes(x = decision, fill = on.target, y = frac)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.75)) +
    geom_text(aes(label = n), position = position_dodge(width = 0.75)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          text = element_text(size = 12)) +
    theme(aspect.ratio = 1) +
    ggsave(plot.out, width = 6, height = 6, dpi = 300)
  invisible(p.frac)
}

nanopore.pore.activity.graph <- function(
  seqsum, pore.col = "pore_id", pore.ids = NULL, read.ids, split.by.pore.id = F,
  x = "start_time", y = "sequence_length_template",
  plot.out = NULL) {
  if (is.character(seqsum)) {
    seqsum <- read.table(seqsum, sep = "\t", header = T)
  }
  if (pore.col == "pore_id") {
    if (! "pore_id" %in% colnames(seqsum)) {
      seqsum$pore_id <- paste0(seqsum$channel, "_", seqsum$mux)
    }
  }
  
  
  if (is.null(pore.ids)) {
    pore.ids <- seqsum %>% dplyr::filter(read_id %in% read.ids) %>% 
      dplyr::pull(!!as.name(pore.col))
  }
  seqsum <- seqsum[, c(x, y, pore.col)] %>% dplyr::filter(!!as.name(pore.col) %in% pore.ids)
  if (split.by.pore.id) {
    split.by <- seqsum[, pore.col] %>% factor(., levels = unique(.))
  } else{
    split.by <- 1
  }
  pl <- seqsum %>% split(f = split.by) %>% 
    lapply(function(seqsum) {
      p <- ggplot(seqsum, aes_string(x = x, y = y, group = pore.col, color = pore.col)) +
        geom_line() +
        geom_point()
      return(p)
    })
  p <- scFanc::wrap.plots.fanc(plot.list = pl, plot.out = plot.out)
  invisible(p)
}

nanopore.death <- function(seqsum, read.type.df, pore.id.col = "pore_id",
                           n.time.blocks = 100, blocks.plot = seq(5, 95, 10),
                           quantile = 0.25, out.dir, root.name = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  if (pore.id.col == "pore_id" && !"pore_id" %in% colnames(seqsum)) {
    seqsum$pore_id <- paste0(seqsum$channel, "_", seqsum$mux)
  }
  read.type.df <- read.type.df[, c("read_id", "read_type")] %>% na.omit()
  seqsum <- seqsum[, c("read_id", "start_time", pore.id.col, "sequence_length_template")]
  utilsFanc::check.dups(read.type.df$read_id)
  utilsFanc::check.dups(seqsum$read_id)
  
  seqsum <- dplyr::inner_join(seqsum, read.type.df, by = "read_id")
  seqsum$time.block <- cut(seqsum$start_time, breaks = n.time.blocks)
  seqsum$time.block.int <- seqsum$time.block %>% as.integer()
  bk <- seqsum
  
  seqsum <- seqsum %>% 
    dplyr::group_by(!!as.name(pore.id.col), read_type, time.block, time.block.int) %>% 
    dplyr::summarise(n.base = sum(sequence_length_template)) %>% 
    dplyr::ungroup() %>% as.data.frame()
  pores <- seqsum[, pore.id.col] %>% unique()
  time.blocks <- seqsum$time.block %>% unique() %>% sort()
  read.types <- seqsum$read_type %>% unique()
  time.frame <- lapply(read.types, function(read.type) {
    time.frame <- data.frame(time.block = rep(time.blocks, each = length(pores)),
                             pore = rep(pores, length(time.blocks)),
                             read_type = read.type)
    colnames(time.frame) <- c("time.block", pore.id.col, "read_type")
    time.frame$time.block.int <- as.integer(time.frame$time.block)
    return(time.frame)
  }) %>% do.call(rbind, .)
  
  seqsum <- dplyr::left_join(
    time.frame, seqsum, by = c("time.block", "time.block.int", pore.id.col, "read_type"))
  seqsum$n.base[is.na(seqsum$n.base)] <- 0
  
  seqsum <- seqsum %>% dplyr::group_by(!!as.name(pore.id.col)) %>% 
    dplyr::arrange(time.block) %>% 
    dplyr::mutate(time.of.death = as.integer(time.block[max(which(n.base != 0))])) %>% 
    dplyr::ungroup()
  
  seqsum <- seqsum %>% dplyr::filter(time.block.int <= time.of.death)
  seqsum <- seqsum %>% dplyr::group_by(!!as.name(pore.id.col), read_type) %>% 
    dplyr::arrange(time.block) %>% 
    dplyr::mutate(n.base.cum = cumsum(n.base)) %>% 
    dplyr::ungroup()
  
  # seqsum <- seqsum %>% filter(time.block.int == 1)

  seqsum <- seqsum %>% dplyr::group_by(read_type, time.block) %>% 
    dplyr::mutate(top.quantile = quantile(n.base.cum, 1-quantile), 
                  bottom.quantile = quantile(n.base.cum, quantile),
                  section = ifelse(n.base.cum >= top.quantile[1], "top", 
                                   ifelse(n.base.cum <= bottom.quantile[1], "bottom", 
                                          NA))) %>% 
    dplyr::ungroup()
  
  seqsum <- as.data.frame(seqsum)
  dir.create(out.dir, showWarnings = F, recursive = T)
  tsv <- paste0(out.dir, "/", root.name, "_seqsum.tsv")
  write.table(seqsum, tsv, sep = "\t", col.names = T, row.names = F, quote = F)
  seqsum <- seqsum[!is.na(seqsum$section),]
  seqsum <- seqsum[seqsum$time.block.int %in% blocks.plot,] %>% 
    dplyr::arrange(time.block.int) %>% 
    dplyr::mutate(time.block = as.character(time.block)) %>% 
    dplyr::mutate(time.block = factor(time.block, levels = sort(unique(time.block))))

  pl <- seqsum %>% split(., f = .$read_type) %>% 
    lapply(function(df) {
      df$time.block.fac <- as.character(df$time.block.int) %>% 
        factor(., gtools::mixedsort(unique(.)))
      section.sum <- df %>% dplyr::group_by(time.block.fac, section) %>% 
        dplyr::summarise(mean.base.cum = utilsFanc::so.formatter(round(mean(n.base.cum), digits = 0))) %>% 
        dplyr::ungroup() %>% as.data.frame()
      section.sum$time.of.death <- 2
      # browser()
      p <- ggplot(df, aes(x = time.block.fac, y = time.of.death)) +
        geom_point(aes(color = section), position = position_jitterdodge(),
                   size = 0.1, alpha = 0.5) +
        geom_violin(aes(x = time.block.fac, color = section)) +
        geom_text(data = section.sum, 
                  aes(label = mean.base.cum, color = section),
                  position = position_dodge(width = 1), size = 2) +
        ggtitle(df$read_type[1]) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      return(p)
    })
  plot.out <- paste0(out.dir, "/", root.name, ".png")
  p <- scFanc::wrap.plots.fanc(pl, plot.out = plot.out, sub.width = length(blocks.plot))
  invisible(p)
}