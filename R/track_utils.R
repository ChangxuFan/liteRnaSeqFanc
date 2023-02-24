bed.2.bw.fanc <- function(bed.vec, bw.vec = NULL, threads = 8, genome, chrom.sizes.suffix = ".convM",
                          scale.to = 10000000, force.count = F, sort = F, run = T,
                          bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools",
                          bedgraph2bigwig = "/opt/apps/kentUCSC/334/bedGraphToBigWig") {
  if (is.null(bw.vec)) 
    bw.vec <- bed.vec %>% sub(".bed$", ".bw", .)
  bdg.vec <- sub(".bw$", ".bedgraph", bw.vec)
  trash <- mclapply(seq_along(bed.vec), function(i) {
    bed <- bed.vec[i]
    bw <- bw.vec[i]
    bdg <- bdg.vec[i]
    chrom.sizes <- paste0("~/genomes/", genome, "/", genome, ".chrom.sizes", chrom.sizes.suffix)
    cmd <- paste0(bedtools, " genomecov -i ", bed, " -g ", chrom.sizes, " -bg ")
    if (!is.null(scale.to)) {
      if (file.exists(paste0(bed, ".wcl")) && force.count == F)
        wc.l <- readLines(paste0(bed, ".wcl"))[1] %>% as.numeric()
      else {
        wc.l <- paste0("wc -l ", bed) %>% system(intern = T) %>% 
          sub(" .+$", "", .) %>% as.numeric() %>% round(digits = 2)
        write(wc.l, paste0(bed, ".wcl"))
      }
      
      scale <- scale.to/wc.l
      cmd <- paste0(cmd, " -scale ", scale)
    }
    
    cmd <- paste0(cmd, " > ", bdg)
    utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run)
    if (!file.exists(bdg))
      stop(paste0(bdg, " failed to generate"))
    if (sort == T) {
      system(paste0("sort -k1,1 -k2,2n ", bdg, " > ", bdg, ".sort"))
      bdg <- paste0(bdg, ".sort")
    } 
    cmd <- paste0(bedgraph2bigwig, " ", bdg, " ", chrom.sizes, " ", bw)
    utilsFanc::cmd.exec.fanc(cmd, intern = F, run = run)
    if (!file.exists(bw))
      stop(paste0(bw, " failed to generate"))
  }, mc.cores = threads, mc.cleanup = T)
  return(bw.vec)
}