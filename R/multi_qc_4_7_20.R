# library(jsonlite)
# library(png)
#












# t$replication$reproducibility$idr %>% as.data.frame()
# t$replication$num_peaks %>% as.data.frame()
#
# m.qc.general(c("~/4dn/nk/fanc/encode/qc/", "~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/qc/"))
#
# m.qc.idr(c("~/4dn/nk/fanc/encode/qc/", "~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/qc/"))



# m.qc.samstat(c("~/4dn/nk/fanc/encode/qc/", "~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/qc/"))

# m.qc.frag.lenth(c("~/4dn/nk/fanc/encode/qc/", "~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/qc/"))








#
# m.qc.get.file.each.library(c("~/4dn/nk/fanc/encode/frag_length/",
#                              "~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/frag_length//"), "*png") %>%
#   m.qc.plot.grid()

#
# my.js <- c("Sp_NK_Ly49A_neg", "Sp_NK_Ly49A_pos", "Sp_NK_Ly49D_neg", "Sp_NK_Ly49D_pos") %>%
#   paste0("/bar/cfan/4dn/nk/fanc/lodge/encode/",., "/",.,".json")
# oshea.js <- c("BM_HSC", "BM_memory_CD8", "Sp_NK")  %>%
#   paste0("~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/",., "/",.,".json")
#
# encode.parser(my.js, assay = "atac", jsd.dir = "/bar/cfan/4dn/nk/fanc/encode/jsd/")
# encode.parser(oshea.js, assay = "atac", jsd.dir = "/bar/cfan/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/jsd/")
# encode.parser(my.js, assay = "atac", gc.dir = "/bar/cfan/4dn/nk/fanc/encode/gc/")
# encode.parser(oshea.js, assay = "atac", gc.dir = "/bar/cfan/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/gc/")
# m.qc.get.file.each.library(c("~/4dn/nk/fanc/encode/peak/",
#                              "~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/peak//"), "*png") %>%
#   m.qc.plot.grid()
#
#
#
#
#
# m.qc.get.file.each.library(c("~/4dn/nk/fanc/encode/jsd/",
#                              "~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/jsd//"), "*png", T) %>%
#   m.qc.plot.sample(3)
#
#
#
#
# t.qc <- c("~/4dn/nk/fanc/encode/qc/", "~/4dn/nk/oshea/lodge/encode_HSC_CD8_NK/qc/")
# m.qc.json.lib.feature.table(t.qc, c("peak_enrich", "frac_reads_in_peaks", "macs2"), replicate.subset.pattern = "^rep\\d$")
# m.qc.json.lib.feature.table(t.qc, c("peak_enrich", "frac_reads_in_peaks", "idr"), replicate.subset.pattern = "^rep\\d_vs_rep\\d$")
# m.qc.json.lib.feature.table(t.qc, c("peak_enrich", "frac_reads_in_annot"))
# tt <- m.qc.json.lib.feature.table(t.qc, c("align", "frag_len_stat"))
