GEO.microarray.2.de <- function(GSE, GPL, sample.map, collapse.method = "mean") {
  # example usage: ~/hmtp/scAR/reference/schuettpelz2014/sz1.1_geo2R_2024-03-22.R
  # The function returns as gene x sample matrix.
  utilsFanc::check.intersect(
    c("GSM", "sample"), "required columns",
    colnames(sample.map), "colnames(sample.map)")
  
  if (!collapse.method %in%  c("mean", "sum")) {
    stop("only mean and sum have been developed for collpase.method")
  }
  # if (collapse.before.log2) {
  #   warning("Collapsing probes before log2 normalization!")
  # }
  
  #>>>>>>>> from GEO2R:
  gset <- GEOquery::getGEO(GSE, GSEMatrix =TRUE, AnnotGPL=TRUE)
  if (length(gset) > 1) idx <- grep(GPL, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  GPL <- gset@experimentData@other$platform_id
  # make proper column names to match toptable 
  Biobase::fvarLabels(gset) <- make.names(Biobase::fvarLabels(gset))
  ex <- Biobase::exprs(gset)
  # > ex[1:3, 1:3]
  #          GSM1329461 GSM1329462 GSM1329463
  # 10338001 6646.02700 7237.18900 7468.84100
  # 10338002   35.70558   34.97189   36.50325
  # 10338003 2080.23600 2249.72900 2751.11800
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  # it seems that it will try to decide if the dataset is already log transformed.
  
  if (LogC) { 
    ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex)
  }
  
  gset <- gset[complete.cases(Biobase::exprs(gset)), ]
  #<<<<<<<< END: from GEO2R
  
  #>>>>>>>>>>>> FANC code starts here:
  mat <- Biobase::exprs(gset)
  
  utilsFanc::check.dups(sample.map$GSM, "sample.map$GSM")
  utilsFanc::check.dups(sample.map$sample, "sample.map$sample")
  utilsFanc::check.intersect(colnames(mat), "GSMs", sample.map$GSM, "sample.map$GSM")
  
  map <- sample.map$sample
  names(map) <- sample.map$GSM
  colnames(mat) <- map[colnames(mat)]
  featureData <- gset@featureData@data
  
  if (!"Gene.symbol" %in% colnames(featureData)) {
    if ("GENE_SYMBOL" %in% colnames(featureData)) {
      colnames(featureData)[colnames(featureData) == "GENE_SYMBOL"] <- "Gene.symbol"
    } else if ("gene_assignment" %in% colnames(featureData)) {
      colnames(featureData)[colnames(featureData) == "gene_assignment"] <- "Gene.symbol"
    } else {
      stop("Can't find gene symbols.")
    }
  }
  gene.map <- featureData[, c("ID", "Gene.symbol")]
  gene.map$ID <- as.character(gene.map$ID)
  if (!identical(as.character(gene.map$ID), rownames(mat))) {
    stop("!identical(as.character(gene.map$ID), rownames(mat))")
  }
  
  mat <- mat[gene.map$Gene.symbol != "",]
  gene.map <- gene.map %>% filter(Gene.symbol != "")
  if (!identical(gene.map$ID, rownames(mat))) {
    stop("!identical(as.character(gene.map$ID), rownames(mat))")
  }
  
  if (any(is.na(mat))) {
    stop("any(is.na(mat))")
  }
  
  if (collapse.method %in% c("mean", "sum")) {
    groupings <- gene.map$Gene.symbol
    names(groupings) <- gene.map$ID
    mat <- scFanc::aggregate.fanc(mat = mat, margin = 1, groupings = groupings, take.mean = (collapse.method == "mean"))
  }
  
  if (any(is.na(mat))) {
    stop("any(is.na(mat))")
  }
  
  ## some gene names are written as something like "Bfar///3110001I22Rik". We remove the Rik part
  # weird.names <- rownames(mat) %>% .[grepl("^..........Rik$|///", .)]
  # name.df <- data.frame(before = weird.names, after = genename.remove.Rik(weird.names))
  # View(name.df)
  rownames(mat) <- rownames(mat) %>% microarray.genename.pick(GPL = GPL)
  mat <- mat[!grepl("^..........Rik\\d*$", rownames(mat)),]
  mat <- mat[!duplicated(rownames(mat)),]
  mat <- mat[!is.na(rownames(mat)),]
  mat <- mat[rownames(mat) != "",]
  return(as.matrix(mat))
}

microarray.genename.pick <- function(x, Rik.pattern = "^..........Rik\\d*$", GPL = "GPL16570") {
  if (GPL %in% c("GPL16570")) {
    strsplit(x, " *///* *") %>% lapply(function(x) {
      return(x[2])
    }) %>% unlist()
  } else {
    strsplit(x, "///") %>% lapply(function(x) {
      bRik <- grepl(Rik.pattern, x)
      if (sum(!bRik) > 0) {
        return(x[!bRik][1])
      } else {
        return(x[bRik][1])
      }
    }) %>% unlist()
  }
  
}

microarray.quantile.norm <- function(de) {
  de <- lapply(de, function(s2b) {
    s2b$bulkNorm <- as.data.frame(qn.fanc(s2b$bulk.mat, T))
    s2b$bulkNorm$gene <- rownames(s2b$bulkNorm)
    return(s2b)
  })
}

microarray.read.geo2r <- function(geo2r.df, symbol.col = "automatic", collapse.probes = "best.p",
                                  remove.Rik = T, remove.Gm = T) {
  if (is.character(geo2r.df)) {
    geo2r.df <- read.table(geo2r.df, header = T, sep = "\t", quote = "")
  }
  
  required.columns <- c("ID", "adj.P.Val", "P.Value", "logFC")
  required.rename <- c("probe", "padj", "pvalue", "log2FoldChange")
  utilsFanc::check.intersect(required.columns, "required columns",
                             colnames(geo2r.df), "colnames(geo2r.df)")
  if (symbol.col == "automatic") {
    candidates <- c("GENE_SYMBOL", "Gene.symbol")
    if (!any(candidates %in% colnames(geo2r.df))) 
      stop("!any(candidates %in% colnames(geo2r.df))")
    symbol.col <- candidates[candidates %in% colnames(geo2r.df)][1]
  } else {
    if (!symbol.col %in% colnames(geo2r.df)) {
      stop("!symbol.col %in% colnames(geo2r.df)")
    }
  }
  geo2r.df <- utilsFanc::change.name.fanc(geo2r.df, cols.from = c(required.columns, symbol.col), 
                              cols.to = c(required.rename, "gene"))
  if (remove.Rik) {
    geo2r.df$gene <- microarray.genename.pick(x = geo2r.df$gene, GPL = "miao")
    geo2r.df <- geo2r.df %>% dplyr::filter(!grepl("^..........[rR]ik\\d*", gene))
  }
  if (remove.Gm == T) {
    geo2r.df <- geo2r.df %>% dplyr::filter(!grepl("^Gm\\d+$", gene))
  }
  geo2r.df <- geo2r.df %>% dplyr::filter(gene != "", !is.na(gene))
  
  geo2r.df$log2FoldChange <- geo2r.df$log2FoldChange / log(2) # dui shu huan di gong shi
  
  if (collapse.probes == "best.p") {
    geo2r.df <- geo2r.df %>% split(geo2r.df$gene) %>% 
      lapply(function(df) {
        return(df[which.min(df$pvalue)[1],])
      }) %>% do.call(rbind, .)
    geo2r.df <- geo2r.df %>% dplyr::filter(gene != "", !is.na(gene))
    
  } else {
    stop("only best.p has been developed.")
  }
  
  rownames(geo2r.df) <- NULL
  geo2r.df$Gene.title <- NULL
  
  return(geo2r.df)
}
