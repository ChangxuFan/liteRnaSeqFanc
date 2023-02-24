#
filter.by.region.dna <- function(bam, region.df,
 what = c("qname","rname", "flag", "pos", "qwidth", "mapq") ,
 tag = c("AS", "XS")){
  # region.df has to use "chr left right" tradition.

  which <- df.to.IrangeList.fanc(region.df)
  what <- c("qname","rname", "flag", "pos", "qwidth", "mapq")

  param <- ScanBamParam(which=which, what=what, tag = tag)

  reads <- scanBam(bam, param=param)
  df <- lapply(reads, bam.chunk.to.df) %>% Reduce(rbind,.) %>% unique() %>%
    mutate(right = pos + qwidth) %>%
    rename(chr = rname, left = pos) %>%
    `[`(, c("chr", "left", "right", "flag", "mapq", "qname", paste0("tag.", tag))) %>%
    arrange(chr, left)
  return(df)
}

# filter.by.region.dna(t.bam, t.klra1)
bam.chunk.to.df <- function(chunk) {
    null.names <- names(chunk$tag[sapply(chunk$tag, is.null)])
    for (i in null.names) {
      chunk$tag[[i]] <- rep(NA, length(chunk$qname))
    }
    chunk <- as.data.frame(chunk)
    return(chunk)
}

filter.by.flag <- function(flag.list, pos) {
  lapply(flag.list, function(x) {
    bits <- R.utils::intToBin(x) %>% strsplit("") %>% unlist() %>% rev()
    if (length(pos) != 1)
      stop("pos must be of exactly 1 element!!")
    if (pos > length(bits))
      return(F)
    else {
      if (bits[pos] == "1")
        return(T)
      else return(F)
    }

  }) %>% unlist() %>% return()
}

df.to.IrangeList.fanc <- function(x) {
  x$id <- 1:nrow(x)
  irange.list <- split(x, f = x$id) %>% lapply(function(y) {
    irange <- IRanges::IRanges(start = y$left, end = y$right, names = y$chr)
    return(irange)
  }) %>% IRanges::IRangesList()
  names(irange.list) <- x$chr
  return(irange.list)
}

# test with the garcia dataset:
# t.bam <- "~/4dn/nk/garcia/bowtie2_aligned/EL4_rep1_ATAC_GSM3573162_garcia_SRR8466980_fastp_filtered_dedup.bam"
# t.klra1 <- data.frame(chr="chr6", left=130382338, right=130382765)
#
# t.which <- df.to.IrangeList.fanc(t.klra1)
# t.what <- c("rname", "flag", "pos", "qwidth", "seq", "qname")
#
# t.param <- ScanBamParam(which=t.which, what=t.what)
# t <- scanBam(t.bam, param=t.param)

throwing.rate.dna <- function(unfiltered.bam, filtered.bam = NULL,
                              region.df) {
  pre <- filter.by.region.dna(unfiltered.bam, region.df = region.df) %>%
    filter(!filter.by.flag(flag, 9)) # filter out secondary alignments in this region.
  post <- filter.by.region.dna(filtered.bam, region.df = region.df)

  df <- data.frame(n.pre = nrow(pre), n.post = nrow(post), throw.rate = 1- nrow(post)/nrow(pre))
  return(list(pre = pre, post = post, stats = df))

}

# my.js <- c("Sp_NK_Ly49A_neg", "Sp_NK_Ly49A_pos", "Sp_NK_Ly49D_neg", "Sp_NK_Ly49D_pos") %>%
#   paste0("/bar/cfan/4dn/nk/fanc/lodge/encode/",., "/",.,".json")
# encode.parser(my.js, assay = "atac", unfiltered.bam = "/bar/cfan/4dn/nk/fanc/encode/unfiltered_bam/")

# throwing.rate.dna(t.bam, t.bam, t.klra1)
# t <- throwing.rate.dna(Sys.glob("~/4dn/nk/fanc/encode/unfiltered_bam/Sp_NK_Ly49A_pos_rep1/*.bam"),
#                   Sys.glob("~/4dn/nk/fanc/encode/bam/Sp_NK_Ly49A_pos_rep1/*.bam"), t.klra1)
# t.pro1 <- data.frame(chr="chr6", left=130386751, right =130387272)
# t <- throwing.rate.dna(Sys.glob("~/4dn/nk/fanc/encode/unfiltered_bam/Sp_NK_Ly49A_pos_rep1/*.bam"),
#                   Sys.glob("~/4dn/nk/fanc/encode/bam/Sp_NK_Ly49A_pos_rep1/*.bam"), t.pro1)
