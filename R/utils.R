fastq.group.gen <- function(fastqs, pe.regex, return.name = F) {
  rm.pe.name <- fastqs %>% sub(pe.regex, "", .)
  if (return.name == T)
    groups <- basename(rm.pe.name)
  else
    groups <- rm.pe.name %>% as.factor() %>% as.numeric()
  return(groups)
}