vcf.ann.extract <- function(vcf, fields, use.LOF = T, use.epistasis = T, remove.up.downstream.variants = T,
                            use.zygocity = F, out.xlsx = NULL) {
  # See specification in https://pcingola.github.io/SnpEff/snpeff/inputoutput/
  # each variant might have multiple lines for multiple annotations.
  field.map <- c("Allele", "Effect", "Putative_Impact", "Gene_Name", "Gene_ID", 
                 "Feature_Type", "Feature_ID", "Transcript_Biotype", "Rank", "HGVSc", 
                 "HGVSp", "cDNA_Pos", "CDS_Pos", "Protein_Pos", "Distance",
                 "Msg")
  epistasis <- list()
  epistasis$Putative_Impact <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  epistasis$Effect <- c(
    "chromosome_number_variation",
    "exon_loss_variant",
    "frameshift_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "rare_amino_acid_variant",
    "missense_variant",
    "disruptive_inframe_insertion",
    "conservative_inframe_insertion",
    "disruptive_inframe_deletion",
    "conservative_inframe_deletion",
    "5_prime_UTR_truncation+exon_loss_variant",
    "3_prime_UTR_truncation+exon_loss",
    "splice_branch_variant",
    "splice_region_variant",
    "stop_retained_variant",
    "initiator_codon_variant",
    "synonymous_variant",
    "initiator_codon_variant+non_canonical_start_codon",
    "coding_sequence_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "5_prime_UTR_premature_start_codon_gain_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TF_binding_site_variant",
    "regulatory_region_variant",
    "miRNA",
    "custom",
    "sequence_feature",
    "conserved_intron_variant",
    "intron_variant",
    "intragenic_variant",
    "conserved_intergenic_variant",
    "intergenic_region",
    "non_coding_transcript_exon_variant",
    "nc_transcript_variant",
    "gene_variant",
    "chromosome"
  )
  
  if (is.numeric(fields)) {
    fields <- field.map[fields]
  }
  
  if (is.character(vcf)) {
    vcf <- vcfR::read.vcfR(vcf)
  }
  fix <- vcf@fix %>% as.data.frame()
  variants <- paste0(fix$CHROM, ":", fix$POS, "_", fix$REF, ">", fix$ALT)
  utilsFanc::check.dups(variants, "variants")
  
  ann <- fix$INFO %>% stringr::str_extract("ANN[^;]+") %>% sub("ANN.", "", .) %>% strsplit(",")
  
  ann.df <- lapply(1:length(ann), function(i) {
    ann <- ann[[i]] %>% strsplit("\\|") %>% lapply(function(x) {
      # if (length(x) > 16) 
      if (length(x) < 16) {
        x <- c(x, rep("", 16-length(x)))
      }
      return(x)
    }) %>% as.data.frame() %>% t()
    rownames(ann) <- NULL
    
    if (ncol(ann) != length(field.map)) {
      stop("ncol(ann) != length(field.map)")
    }
    colnames(ann) <- field.map
    ann <- ann[, fields, drop = F] %>% unique()
    df <- cbind(data.frame(variant = variants[i]), ann)
    return(df)
  }) %>% do.call(rbind, .)
  
  ann.df <- tidyr::separate_rows(ann.df, !!!syms(fields), sep = "&") %>% as.data.frame()
  
  if (remove.up.downstream.variants && "Effect" %in% fields) {
    cat(paste0("Removing upstream_gene_variant and downstream_gene_variant")); cat("\n")
    ann.df <- ann.df[!ann.df$Effect %in% c("upstream_gene_variant", "downstream_gene_variant"),]
  }
  
  if (any(fields %in% names(epistasis)) && use.epistasis) {
    if (!"Gene_Name" %in% fields) {
      stop("Gene_Name must be in fields to use epistasis")
    }
    epistasis.fields <- names(epistasis) %>% .[.%in% fields]
    for (efield in epistasis.fields) {
      utilsFanc::check.intersect(unique(ann.df[, efield]), "ann.df[, efield]", epistasis[[efield]], "epistasis[[efield]]")
      ann.df$level <- factor(ann.df[, efield], levels = epistasis[[efield]]) %>% as.numeric()
      ann.df <- ann.df %>% dplyr::group_by(variant) %>% dplyr::filter(level == min(level)) %>% 
        dplyr::ungroup() %>% as.data.frame()
      ann.df$level <- NULL
    }
    
  }
  
  if (use.LOF) {
    LOF.df <- data.frame(variant = variants,
                         LOF = stringr::str_extract(fix$INFO, "LOF[^;]+") %>% gsub("LOF..|\\)", "", .),
                         NMD = stringr::str_extract(fix$INFO, "NMD[^;]+") %>% gsub("NMD..|\\)", "", .))
    # LOF.df <- data.frame(variant = variants,
    #                      LOF = grepl("LOF[^;]+", fix$INFO),
    #                      NMD = grepl("NMD[^;]+", fix$INFO))
    ann.df <- dplyr::left_join(ann.df, LOF.df, by = "variant")
  }
  
  if (use.zygocity) {
    zy <- cbind(data.frame(variant = variants), sub(":.+", "", vcf@gt[, -1]))
  }
  ann.df <- dplyr::left_join(ann.df, zy, by = "variant")
  if (!is.null(out.xlsx)) {
    dir.create(dirname(out.xlsx), showWarnings = F, recursive = T)
    xlsx::write.xlsx2(x = ann.df, file = out.xlsx, sheetName = "sheet1", col.names = T, row.names = F)
  }
  return(ann.df)
}