ta2bigwig <- function(in.ta, genome, run=T, bedgraph2bigwig = "/opt/apps/kentUCSC/334/bedGraphToBigWig",
                      bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools") {
  out.bedgraph <- sub("tagAlign.gz", "bedgraph", in.ta)
  out.bw <- sub("tagAlign.gz", "bw", in.ta)
  cmd <- paste0("zcat ", in.ta, " | sort -k 1,1 | ",bedtools, " genomecov", " -bg",
                " -i stdin", " -g ", "~/genomes/",genome, "/",genome, ".chrom.sizes",
                " | sort -k1,1 -k2,2n > ", out.bedgraph)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run)
  
  cmd <- paste0(bedgraph2bigwig, " ", out.bedgraph, " ", "~/genomes/",genome, "/",genome, ".chrom.sizes", " ", out.bw)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run)
  
  if (!file.exists(out.bw))
    stop(paste0(out.bw, " failed to generate"))
  return(out.bw)
}

# ta2bigwig(in.ta = Sys.glob("~/4dn/nk/fanc/encode/ta/Sp_NK_Ly49A_neg_rep1/*.tagAlign.gz"), genome = "mm10", run = F)
