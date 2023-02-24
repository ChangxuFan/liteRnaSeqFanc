# this script is called from each of the individual directories (such as ~/4dn/nk/chip/sth/assorted) to generate all alignments
library(rnaSeqFanc)

options(stringsAsFactors = F)

t <- dna.rna.seq(fastq.info = "sample.info.tsv", type = "dna", genome.index = "~/genomes/mm10/bowtie2/mm10", qc = T,
                 pe.rmdup = T, se.rmdup = F,
                 skip.qc = T, skip.align = F, skip.dedup = F, generate.bw = T, mc.cores.GSM = 4, mc.cores.SRR = 2,
                 filter.expression.se = FILTER.LIGHT, filter.expression.pe = FILTER.LIGHT, thread = 4)
