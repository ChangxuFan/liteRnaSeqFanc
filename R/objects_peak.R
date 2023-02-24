narrow.peak.watch <- function(np, key, n, min=NULL, max=NULL, out.file = NULL, cols.add = paste0("V", 7:10)) {
  if (is.character(np)) {
    np <- read.table(np, as.is = T, header = F)
  }
  if (!is.null(min)) {
    np <- np %>% .[.[,key] > min,]
  }
  if (!is.null(max)) {
    np <- np %>% .[.[,key] < min,]
  }
  np <- np[sample(1:nrow(np), n , replace = F),]
  np <- np[order(np[, key]), ]
  info <- c(key, cols.add) %>% unique()
  # np <- np %>% mutate(info = paste0(!!as.name(info), collapse = "|")) %>% 
  #   select(V1, V2, V3, info)
  np$info <- apply(X = np[,info],MARGIN = 1, FUN = function(x) {
    # print("miao")
    paste0(x[info], collapse = "|")
  })
  np <- np[, c("V1", "V2", "V3", "info")]
  if (!is.null(out.file)) {
    write.table(np, out.file, quote = F, row.names = F, col.names = F, sep = "\t")
    system(paste0("/bar/cfan/scripts/bed_browser_v2.sh ", out.file ))
  }
  return(np)
}