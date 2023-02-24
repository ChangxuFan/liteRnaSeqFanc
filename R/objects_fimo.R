fimo.parse.ez <- function(fimo.dir, out.rgbPeak, genome,
                          bedToBigBed = "/bar/cfan/software/kent/bin/bedToBigBed") {
  # region names should be simple genomic coordinates in the form of 
  # chr:start-end. This is what you get from bedtools getfasta by default
  # however, the default fasta header from bedtools getfasta is written in 1 based, 
  # despite in chr:start-end format
  df <- read.table(paste0(fimo.dir, "/fimo.tsv"), header = T)
  df.rgb <- utilsFanc::loci.2.df(loci.vec = df$sequence_name, remove.loci.col = T)
  df.rgb <- df.rgb %>% mutate(end  = start + df$stop, start = start + df$start,  # should be start = start + 1 + df$start -1
                    name = df$matched_sequence, score = round(-log10(df$p.value), digits = 0), 
                    strand = df$strand, thickStart = start, thickEnd = end, rgb = '0,255,0')
  df.rgb$start <- df.rgb$start - 1
  
  df.rgb$rgb[df.rgb$strand == "-"] <- '255,0,0'
  dir.create(dirname(out.rgbPeak), showWarnings = F, recursive = T )
  df.rgb <- df.rgb %>% arrange(chr, start)
  temp.bed <- paste0(out.rgbPeak, ".bed")
  write.table(df.rgb, temp.bed, sep = "\t", quote = F, col.names = F, row.names = F)
  chrom.sizes <- paste0("~/genomes/", genome, "/", genome, ".chrom.sizes")
  cmd <- paste0(bedToBigBed, " ", temp.bed, " ", chrom.sizes, " ", out.rgbPeak)
  system(cmd)
  cat(utilsFanc::bash2ftp(filename = out.rgbPeak))
  return()
}