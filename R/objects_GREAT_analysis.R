great.plot.bar <- function(rgreat.res.table, padj.col = "Binom_Adjp_BH", top.n = 10,
                           text.size = 0.36 * 6, string.wrap.width = 40, title = "",
                           plot.out, width = 2, height = 1, color = "darkorchid4") {
  df <- rgreat.res.table
  df$x <- df[, padj.col]
  df$name <- df$name %>% stringr::str_wrap(width = string.wrap.width)
  
  if (!is.null(top.n)) {
    df <- df %>% arrange(x)
    df <- df[1:min(top.n, nrow(df)),]
  }
  
  df$name <- factor(df$name, levels = rev(df$name))
  
  df$x <- -log10(df$x)
  
  p <- ggplot(df, aes(x = x, y = name)) +
    geom_bar(stat = "identity", fill = color, color = color, width = 0.5) +
    xlab("-log10(FDR)") +
    ylab("GO gene set") +
    ggtitle(title)
  
  p <- p %>% utilsFanc::theme.fc.1(italic.x = F, remove.axis.titles = F) +
    theme(plot.margin = unit(c(0.01, 0, 0, 0.01), "in"))
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = F)
    ggsave(plot.out, p, device = cairo_pdf, width = width, height = height)
  }
  invisible(p)
}