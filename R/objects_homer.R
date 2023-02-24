homer.motif.enrich <- function(fg.bed, genome, outdir, bg.bed = NULL, simplify.bed, size, thread, denovo = F,
                               find.motif = NULL, find.motif.name = NULL, other.params = "",
                               find.motifs.path = "/bar/cfan/software/homer/bin/findMotifsGenome.pl",
                               homer.path = "/bar/cfan/software/homer/bin/",
                               run = T, stdout.file = NULL) {
  
  system(paste0("mkdir -p ", outdir))
  PATH.bk <- Sys.getenv("PATH")
  if (!grepl("software/homer/bin", PATH.bk))
    Sys.setenv(PATH = paste0(PATH.bk, ":", homer.path))
  if (!is.character(fg.bed)) {
    fg.bed <- as.data.frame(fg.bed %>% `names<-`(NULL))
    if (simplify.bed == T)
      fg.bed <- fg.bed[, 1:3] %>% cbind(strand = "+")
    fg.bed.file <- tempfile()
    write.table(fg.bed, fg.bed.file, quote = F, sep = "\t", col.names = F, row.names = F)
  } else {
    fg.bed.file <- fg.bed
  }
  
  cmd <- paste0(find.motifs.path, " ", fg.bed.file, " ","/bar/cfan/genomes/",genome, "/", genome, ".fa", 
                " ", outdir, " -size ", size, " -p ", thread)
  if (denovo == F) {
    cmd <- paste0(cmd, " -nomotif ")
  }
  if (!is.null(bg.bed)) {
    if (!is.character(bg.bed)) {
      bg.bed <- as.data.frame(bg.bed %>% `names<-`(NULL))
      if (simplify.bed == T)
        bg.bed <- bg.bed[, 1:3] %>% cbind(strand = "+")
      bg.bed.file <- tempfile()
      write.table(bg.bed, bg.bed.file, quote = F, sep = "\t", col.names = F, row.names = F)
    } else {
      bg.bed.file <- bg.bed
    }
    
    cmd <- paste0(cmd, " -bg ", bg.bed.file)
  }
  
  cmd <- paste0(cmd, " ", other.params)
  if (!is.null(find.motif)) {
    if (is.null(find.motif.name))
      find.motif.name <- basename(find.motif) %>% sub(".motif", "", .)
    cmd <- paste0(cmd, " -find ", find.motif, " > ", outdir, "/", find.motif.name, "_hits.tsv")
    stdout.file <- NULL
  }
  try(utilsFanc::cmd.exec.fanc(cmd, intern = F, run = run, stdout.file = stdout.file))
  Sys.setenv(PATH = PATH.bk)
  
  return()
}

homer.scan.parse <- function(hit.tsv, hit.name = NULL, in.bed, out.hit.regions = NULL, out.hit.motifs = NULL) {
  # output bed files
  if (is.null(hit.name)) {
    hit.name <- sub("_hits.tsv", "", basename(hit.tsv))
  }
  if (is.null(out.hit.regions)) {
    out.hit.regions <- tools::file_path_sans_ext(basename(in.bed)) %>% 
      paste0(dirname(hit.tsv), "/",  ., "_", hit.name, "_regions.bed")
  }
  if (is.null(out.hit.motifs)) {
    out.hit.motifs <- tools::file_path_sans_ext(basename(in.bed)) %>% 
      paste0(dirname(hit.tsv), "/",  .,"_", hit.name, "_motifs.bed")
  }
  
  df.hits <- read.table(hit.tsv, header = T, sep = "\t", quote = "")
  df.in <- utilsFanc::import.bed.fanc(in.bed, no.shift = F)
  trash <- utilsFanc::write.zip.fanc(df.in[df.hits$PositionID %>% unique(),], out.hit.regions, bed.shift = T)
  df.in$PositionID <- 1:nrow(df.in)
  j <- left_join(df.hits, df.in)
  df.hit.motifs.pos <- j %>% filter(Strand == "+") %>% mutate(pos = start) %>% mutate(start = pos + Offset) %>% mutate(end = start + nchar(Sequence) - 1) %>% 
    select(chr, start, end, MotifScore, Motif.Name, Strand)
  df.hit.motifs.neg <- j %>% filter(Strand == "-") %>% mutate(pos = start) %>% mutate(end = pos + Offset) %>% mutate(start = end - nchar(Sequence) + 1) %>% 
    select(chr, start, end, MotifScore, Motif.Name, Strand)
  df.hit.motifs <- rbind(df.hit.motifs.pos, df.hit.motifs.neg)
  trash <- utilsFanc::write.zip.fanc(df.hit.motifs, out.hit.motifs, bed.shift = T)
  return()
}


