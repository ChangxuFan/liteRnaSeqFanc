cutrun.json.gen <- function(master.dir, names = NULL) {
  if (is.null(names))
    names <- list.dirs(path = master.dir, full.names = F, recursive = F)
  df <- lapply(names, function(name) {
    cat(paste0("sample: ", name, "\n"))
    files <- list()
    files$macs2.narrow.bw <- paste0(name, "/", "peakcalling/macs2.narrow/", name, ".cpm.norm.bw")
    files$macs2.narrow.peak <- paste0(name, "/", "peakcalling/macs2.narrow/", name, "_peaks.narrowPeak")
    files$macs2.narrow.peak.rmbl <- paste0(name, "/", "peakcalling/macs2.narrow/blacklist_filtered/",
                                           name, "_peaks.narrowPeak")
    
    files$macs2.broad.bw <- paste0(name, "/", "peakcalling/macs2.broad/", name, ".cpm.norm.bw")
    files$macs2.broad.peak <- paste0(name, "/", "peakcalling/macs2.broad/", name, "_peaks.broadPeak")
    
    files$seacr.bw <- paste0(name, "/", "peakcalling/seacr/", name, ".cpm.norm.bw")
    files$seacr.peak <- paste0(name, "/", "peakcalling/seacr/", name, "_treat.stringent.sort.bed")
    
    df <- data.frame(name = paste0(name, "::",names(files)), url = unlist(files))
    zip.idx <- which(grepl("Peak$|bed$", df$url))
    
    for (i in zip.idx) {
      system(paste0("~/scripts/bed_browser_v2.sh ", master.dir, "/", df$url[i])) 
      df$url[i] <- paste0(df$url[i], ".gz")
    }
    return(df)
  }) %>% Reduce(rbind, .)
  
  
  out.file <- paste0(master.dir, "/tracks.json")
  json <- utilsFanc::jsongen(df = df, outfile = out.file)
  utilsFanc::bash2ftp(out.file) %>% cat(sep = "\n")
  return(out.file)
}



