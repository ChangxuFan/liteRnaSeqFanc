great.plot.bar <- function(rgreat.res.table, rank.by = "Binom_Adjp_BH", x.lab = "-log10(FDR)", minus.log10 = T,
                           top.n = 10, desc = F,
                           text.size = 0.36 * 6, string.wrap.width = 40, title = "",
                           plot.out, width = 2, height = 1, color = "darkorchid4") {
  df <- rgreat.res.table
  df$x <- df[, rank.by]
  df$name <- df$name %>% stringr::str_wrap(width = string.wrap.width) %>% 
    gsub("\n", "    \n", .)
  
  if (!is.null(top.n)) {
    if (desc) {
      df <- df %>% arrange(desc(x))
    } else {
      df <- df %>% arrange(x)
    }
    
    df <- df[1:min(top.n, nrow(df)),]
  }
  
  df$name <- factor(df$name, levels = rev(df$name))
  if (minus.log10)
    df$x <- -log10(df$x)
  
  p <- ggplot(df, aes(x = x, y = name)) +
    geom_bar(stat = "identity", fill = color, color = color, width = 0.5) +
    xlab(x.lab) +
    ylab("GO gene set") +
    ggtitle(title)
  
  p <- p %>% utilsFanc::theme.fc.1(italic.x = F, remove.axis.titles = F) +
    theme(plot.margin = unit(c(0.01, 0.1, 0, 0.01), "in"),
          axis.text.y = element_text(lineheight = 0.8))
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = F)
    ggsave(plot.out, p, device = cairo_pdf, width = width, height = height)
    write.table(df, paste0(tools::file_path_sans_ext(plot.out), ".tsv"), 
                col.names = T, row.names = F, sep = "\t", quote = F)
  }
  invisible(p)
}

great.plot.heatmap <- function(rgreat.table.list, pathways = NULL,
                               binom.FDR.cutoff = 0.05, 
                               column.color = "Binom_Fold_Enrichment",
                               cluster.cols = F, hm.column.order = NULL,
                               hm.colors = c("white", "deeppink4"), hm.values = NULL,
                               scale.rows = F, cluster.rows = F, 
                               add.FDR = T, string.wrap.width = 40,
                               plot.out, width = 2, height = 2, ...) {
  columns <- c("ID", "name", 
               "Binom_Genome_Fraction", "Binom_Expected", 
               "Binom_Observed_Region_Hits", "Binom_Fold_Enrichment", 
               "Binom_Region_Set_Coverage", "Binom_Raw_PValue", "Binom_Adjp_BH", 
               "Hyper_Total_Genes", "Hyper_Expected", "Hyper_Observed_Gene_Hits", 
               "Hyper_Fold_Enrichment", "Hyper_Gene_Set_Coverage", "Hyper_Term_Gene_Coverage", 
               "Hyper_Raw_PValue", "Hyper_Adjp_BH")
  utilsFanc::check.intersect(column.color, "column.color", columns, "available columns")
  
  if (is.null(names(rgreat.table.list))) {
    stop("is.null(names(rgreat.table.list))")
  }
  
  utilsFanc::check.dups(names(rgreat.table.list), "names(rgreat.table.list)")
  
  res <- lapply(names(rgreat.table.list), function(table.name) {
    res <- rgreat.table.list[[table.name]]
    utilsFanc::check.intersect(columns, "required columns", colnames(res), "colnames(res)")
    res <- res %>% dplyr::rename(pathway = name) %>% 
      dplyr::mutate(name = table.name)
  }) %>% do.call(rbind, .)
  
  utilsFanc::check.intersect(pathways, "pathways", res$pathway, "res$pathway")
  res <- res %>% dplyr::filter(pathway %in% pathways)
  res$pathway <- res$pathway %>% stringr::str_to_title() %>% 
    stringr::str_wrap(width = string.wrap.width) %>% 
    gsub("\n", "\n  ", .)
  
  utilsFanc::check.dups(paste0(res$name, res$pathway), "paste0(res$name, res$pathway)")
  
  mat <- reshape2::acast(res, pathway ~ name, value.var = column.color)
  
  if (scale.rows) {
    mat <- t(mat) %>% scale() %>% t()
    hm.colors <- NULL
    hm.values <- NULL
  }
  
  if (any(is.na(mat))) {
    stop("any(is.na(mat))")
  }
  
  if (!is.null(hm.column.order)) {
    utilsFanc::check.intersect(hm.column.order, "hm.column.order", colnames(mat), "colnames(mat)")
    mat <- mat[, hm.column.order]
  }
  
  if (add.FDR) {
    symbol.mat <- reshape2::acast(res, pathway ~ name, value.var = "Binom_Adjp_BH")
    if (!is.null(hm.column.order)) {
      symbol.mat <- symbol.mat[, hm.column.order]
    }
    symbol.mat[symbol.mat < binom.FDR.cutoff] <- "*"
    tmp <- dimnames(symbol.mat)
    symbol.mat <- stringr::str_extract(symbol.mat, "\\*+") %>%
      matrix(nrow = length(tmp[[1]]))
    dimnames(symbol.mat) <- tmp
    
    symbol.mat[is.na(symbol.mat)] <- ""
    
  } else {
    symbol.mat <- NULL
  }
  
  if (is.null(hm.values) && !scale.rows) {
    hm.values <- c(0, max(mat))
  }
  
  plot.mat.rank.row(mat = mat, symbol.mat = symbol.mat,
                    no.col.cluster = !cluster.cols, show_column_names = T,
                    hm.colors = hm.colors, hm.values = hm.values,
                    show_row_names = T, 
                    plot.out = plot.out, width = width, height = height,
                    row_name_fontSize = 5, cluster_rows = cluster.rows,
                    ...)
  
  write.table(res, paste0(tools::file_path_sans_ext(plot.out), ".tsv"), sep = "\t",
              row.names = F, col.names = T, quote = F)
  
  invisible(res)
  
}


list.expand <- function(dense.list) {
  if (is.null(names(dense.list))) {
    stop("is.null(names(dense.list))")
  }
  
  utilsFanc::check.dups(names(dense.list), "names(dense.list)")
  
  res <- list(c())
  for (name in rev(names(dense.list))) {
    n <- dense.list[[name]] %>% length()
    res <- lapply(1:n, function(i) {
      res <- lapply(res, function(res.sub) {
        to.add <- dense.list[[name]][[i]]
        names(to.add) <- name
        res.sub <- c(to.add, res.sub)
      }) 
      return(res)
    }) %>% Reduce(c, .)
  }
  return(res)
}
