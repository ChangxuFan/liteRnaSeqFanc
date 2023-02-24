FILTER.LIGHT <- "mapping_quality > 0 and (not secondary_alignment)"
FILTER.LIGHT.1 <- "mapping_quality > 1 and (not secondary_alignment)"
FILTER.LIGHT.8 <- "mapping_quality > 8 and (not secondary_alignment)"
FILTER.LIGHT.30 <- "mapping_quality > 30 and (not secondary_alignment)"

FILTER.FULLMATCH <- "[AS] == 0 and (not secondary_alignment)"
ATAC.ADAPTER.1="CTGTCTCTTATACACATCT"
ATAC.ADAPTER.2="CTGTCTCTTATACACATCT"
CUTADAPT.PARAMS.TARGET <- "--quality-cutoff=15,10 --minimum-length=36"
TE.CLASS.FILTER <- c("SINE", "LINE", "LTR", "DNA")
BOWTIE2.DEFAULTS <- " --xeq --dovetail"

host.name <- Sys.info()["nodename"][1]
if (grepl("ris", host.name)) {
    # ris:
    SAMTOOLS <- "samtools"
    BWA <- "bwa"
    BOWTIE2 <- "bowtie2"
    BOWTIE2.BUILD <- "bowtie2-build"
    SAMBAMBA <- "~/software/sambamba/sambamba"
} else {
    SAMTOOLS = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools"
    BWA = "~/software/bwa/bwa_0.7.17-r1198/bwa"
    BWA.TARGET <- "~/software/bwa/bwa_0.7.16a-r1181/bwa-0.7.16/bwa"
    BOWTIE2 <- "/opt/apps/bowtie2/2.3.4.1/bowtie2"
    BOWTIE2.BUILD <- "/opt/apps/bowtie2/2.3.4.1/bowtie2-build"

    # BOWTIE2 <- "~/software/Bowtie2/bowtie2-2.4.2-linux-x86_64/bowtie2"
    # BOWTIE2.BUILD <- "~/software/Bowtie2/bowtie2-2.4.2-linux-x86_64/bowtie2-build"
    

    SAMBAMBA <- "~/software/sambamba/sambamba"

    TABIX <- "/bar/cfan/anaconda2/envs/jupyter/bin/tabix"
    BGZIP <- "/bar/cfan/anaconda2/envs/jupyter/bin/bgzip"

    BEDTOOLS.DIR <- "/bar/cfan/anaconda2/envs/jupyter/bin/"
    SEQTK <- "/bar/cfan/anaconda2/envs/jupyter/bin/seqtk"
}



