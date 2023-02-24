# functions used for atacseq analyses downstream of peak calling and tagAlign file generation.
peak.rmsk <- function(narrowPeak, genome = NULL, rmsk.bed = NULL, conv.only = T, summit.in.rmsk = F,
                      out.peak = NULL, out.rmsk = NULL,  bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools",
                      add.param.peak = "", add.param.rmsk = "") {
  if (conv.only == T) {
    simple <- "simple."
  } else {
    simple <- ""
  }
  if (!is.null(genome))
    rmsk.bed <- paste0("~/genomes/",genome, "/rmsk/", genome, ".rm.", simple,"bed")
  
  if (is.null(out.peak))
    out.peak <- utilsFanc::insert.name.before.ext(name = narrowPeak, insert = "rmsk", delim = ".")
  if (is.null(out.rmsk))
    out.rmsk <- paste0(narrowPeak, ".rmsk")
  
  cmd <- paste0(bedtools, " intersect ", add.param.peak, " -a ", narrowPeak, " -b ", rmsk.bed,
                " -wa > ", out.peak)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)
  
  cmd <- paste0(bedtools, " intersect ", add.param.rmsk, " -b ", narrowPeak, " -a ", rmsk.bed,
                " -wa > ", out.rmsk)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)
  
  cmd <- paste0("~/scripts/bed_browser_v2.sh ", out.peak)
  try(system(cmd))
  cmd <- paste0("~/scripts/bed_browser_v2.sh ", out.rmsk)
  try(system(cmd))
  
  if (!file.exists(out.peak))
    stop(paste0(out.peak, " failed to generate"))
  return(out.peak)
}