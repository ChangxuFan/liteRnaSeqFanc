flowjo.table.read <- function(file, colname.file = NULL) {
  # perhaps not a very useful function
  
  # colname.file: a file with 2 columns: 
  # c1: column names from flowjo (format usually: lymphcyte/single_cells/Tcells/...)
  # c2: column names you want to give
  if (!is.null(colname.file)) {
    df <- read.table(file, header = F, sep = "\t")
    names.ori <- unlist(df[1, 2:ncol(df)])
    names(names.ori) <- NULL
    rename.df <- data.frame(ori = names.ori)
    col.df <- read.table(colname.file, sep = "\t")
    colnames(col.df) <- c("ori", "new")
    rename.df <- rename.df %>% dplyr::left_join(col.df, by = "ori")
    if (!identical(names.ori, rename.df$ori)) {
      stop("!unlist(df[1, 2:ncol(df)])")
    }
    df[1, 2:ncol(df)] <- rename.df$new
    df[1,1] <- "Sample"
    colnames(df) <- df[1, ]
    df <- df[2:nrow(df), ]
    rm(col.df); rm(rename.df); rm(names.ori)
  } else {
    df <- read.table(file, header = F, sep = "\t", header = T)
    colnames(df) <- sub("^Sample.$", "Sample", colnames(df))
  }
  df <- df %>% dplyr::filter(! Sample %in% c("Mean", "SD"))
  df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df)], as.numeric)
  return(df)
}

flowjo.colnames <- function(file, write = F, out.file = NULL) {
  # flowjo by default gives weird column names (format usually: lymphcyte/single_cells/Tcells/...)
  # this function just extract these names to put into a column. That's it.
  x <- readLines(file, n = 1)
  x <- x %>% strsplit("\t") %>% unlist()
  x <- x[-1]
  if (write) {
    if (is.null(out.file)) {
      out.file <- tools::file_path_sans_ext(file) %>% paste0("_columns.tsv")
    }
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    write(x, out.file, sep = "\t")
  }
  invisible(x)
}

flowjo.table.read.m <- function(files, rename.dup = F, sampleName.regex.fun = NULL) {
  # files: for example the tables from stain 1/2/3
  # metadata: tsv files; the "sample" column must be present. could also be a df
  dfs <- lapply(files, function(file) {
    df <- read.table(file, header = F, sep = "\t", quote = "")
    df <- df[1:(nrow(df) - 2), ]
    if (all(is.na(df[, ncol(df)]))) {
      df <- df[, -ncol(df)]
    }
    
    gates <- df[1, -1] %>% unlist()
    df <- df[-1, ]
    gates <- strsplit(gates, "(subset)*/| \\| ")
    gates <- lapply(gates, function(x) {
      if (grepl("/Q\\d", x[length(x) - 1])) {
        x <- x[length(x)-1] %>% stringr::str_extract(pattern = ": .+[\\+\\-] , .+[\\+\\-]") %>% 
          gsub(" +", "", .) %>% gsub("[:,]", "", .)
        return(x)
      } else {
        return(x[length(x) - 1])
      }
    }) %>% unlist() %>% 
      sub("^/", "", .) %>% gsub(" +", "", .)
    colnames(df) <- c("sample", gates)
    df$sample <- df$sample %>% sub("Panel\\d+_", "", .)
    
    if (!is.null(sampleName.regex.fun)) {
      df$sample <- sampleName.regex.fun(df$sample)
    }
    
    df[,2:length(df)] <- df[, 2:length(df), drop = F] %>% lapply(as.numeric)
    return(df)
  })

  gates <- lapply(dfs, function(df) return(colnames(df)[-1])) %>% unlist()
  dups <- gates %>% .[duplicated(gates)]
  if (length(dups) > 0) {
    if (!rename.dup) {
      stop(paste0("duplicated gate names found: ", 
                  paste0(unique(dups), collapse = ", ")))
    }
  }
  df <- Reduce(function(x, y) left_join(x, y, by = "sample"), dfs)
  return(df)
}
