#
target.jsongen <- function(target.out.full.paths, samples=NULL, out.json=NULL, prefix, suffix, track_type) {
  # prefix refers to what file you are trying to make a json out of.
  ##something like "step3.2_rmbl_". suffix is a similar idea, such as "peaks.narrowPeak".
  # target.out.full.paths is a vector of paths.
  # samples is a LIST of vectors of sample names, each element of samples corespond to each element of target.out.full.paths.
  # the default behavior is to take all of the directories starting with "Processed_" under each dir in target.out.full.paths.
  # if samples is specified, it will filter for only the samples from the sample list.
  files <- lapply(seq_along(target.out.full.paths), function(i) {
    processed_dirs <- system(paste0("ls -d ", target.out.full.paths[i], "/Processed_*"), intern = T) %>% basename()
    #print(processed_dirs)
    if (!is.null(samples[[i]])) {
      dirs_samples <- sapply(processed_dirs, function(x) {
        any(sapply(samples[[i]], grepl, x))
      })
      processed_dirs <- processed_dirs[dirs_samples]
    }

    files <- paste0(target.out.full.paths[i], "/", processed_dirs, "/", prefix, "*", suffix) %>% Sys.glob()
  }) %>% unlist()

  jsongen <- data.frame(name = basename(files), url = utilsFanc::bash2ftp(files), type = track_type)
  json <- jsongen %>% jsonlite::toJSON() %>% jsonlite::prettify()
  if (!is.null(out.json))
    write(json, out.json)
  return(jsongen)
}


# target.jsongen("~/4dn/nk/fanc/sth/target_pipe/",
#                samples = list(c("Ly49A_neg_rep1", "Ly49A_neg_rep2", "Ly49A_neg_rep3",
#                                 "Ly49A_pos_rep1", "Ly49A_pos_rep2", "Ly49A_pos_rep3")),
#                prefix = "", suffix = "bigWig", track_type = "bigwig",
#                out.json = "~/4dn/nk/fanc/sth/target_pipe/Ly49A_bigwig.json")

qatacview.jsongen <- function(df, add.to = NULL, out.json=NULL, replace=F, help = T) {
  if (help) {
    stop("takes in a tsv file with 3 columns: collection; pipedir, samples. It looks something like this:
   collection                        pipedir        samples
   NK_fanc /bar/cfan/4dn/nk/fanc/sth/target_pipe/ Ly49A_neg_rep1
   NK_fanc /bar/cfan/4dn/nk/fanc/sth/target_pipe/ Ly49A_neg_rep2
   NK_fanc /bar/cfan/4dn/nk/fanc/sth/target_pipe/ Ly49A_neg_rep3
   NK_fanc /bar/cfan/4dn/nk/fanc/sth/target_pipe/ Ly49A_pos_rep1
   NK_fanc /bar/cfan/4dn/nk/fanc/sth/target_pipe/ Ly49A_pos_rep2
   NK_fanc /bar/cfan/4dn/nk/fanc/sth/target_pipe/ Ly49A_pos_rep3
         ")
  }
  if (is.character(df))
    df <- read.table(df, header = T, as.is = T, sep="\t")
  jsongen <- list()
  jsongen[["allOptions"]] <- lapply(df$collection %>% unique(), function(x) {
    list(value = x, label = x, clearableValue=F)
  })
  jsongen[["allProducts"]] <- split(df, f = factor(df$collection, levels = df$collection %>% unique())) %>%
    lapply(function(x) {
      print(x)
      collection.df <- split(x, f = factor(x$pipedir, levels=x$pipedir %>% unique())) %>% lapply(function(y) {
          pipe.df <- split(y, f = factor(y$sample, levels=y$sample %>% unique())) %>% lapply(function(z) {
            if (nrow(z) != 1)
              stop("something went wrong with the splitting of dataframe. z doesn't have exactly one row")
            url <- target.jsongen(z$pipedir[1], samples = list(z$samples),
                                    prefix = "step3.2_", suffix = "normalized*bigWig", track_type = "bigwig") %>% pull(url)
            # print(url)
              file <- target.jsongen(z$pipedir[1], samples = list(z$samples),
                                     prefix = "QC", suffix = "json", track_type = "bigwig") %>% pull(url)
              #print(file)
              data.frame(sample = z$sample, url = url, assay = "ATAC-seq", file = file) %>% return()
            }) %>% Reduce(rbind,.)

        }) %>% Reduce(rbind,.)
      collection.df$id <- 1:nrow(collection.df)
      return(collection.df)
    })
  if (!is.null(add.to)) {
    if (!file.exists(add.to)) {
      stop(paste0(add.to, " does not exist"))
    }
    timestamp <- Sys.time() %>% as.character() %>% sub(" +", "_", .)
    bk <- paste0("/bar/cfan/software/browser/qATACviewer/frontend/src/bk/data.json.bk.", timestamp)
    
    system(paste0("cp ", add.to, " ", bk))
    json.out <- jsonlite::read_json(path = add.to)
    jsongen$allProducts <- lapply(jsongen$allProducts, function(df) {
      res <- df %>% split(f = 1:nrow(df)) %>% lapply(as.list) %>% 
        `names<-`(NULL)
      return(res)
    })
    json.out$allOptions <- c(jsongen$allOptions, json.out$allOptions)
    json.out$allProducts <- c(jsongen$allProducts, json.out$allProducts)
    json <- json.out %>% jsonlite::toJSON(auto_unbox = T) %>% jsonlite::prettify()
    if (replace == T)
      write(json, add.to)
    if (!is.null(out.json))
      write(json, out.json)
  } else {
    json <- jsongen %>% jsonlite::toJSON() %>% jsonlite::prettify()
    if (!is.null(out.json))
      write(json, out.json)
    if (replace == T) {
      # system(paste0("mv ", "~/software/browser/qATACviewer/frontend/src/data.json ",
      #               "~/software/browser/qATACviewer/frontend/src/data.json.bk",
      #               " && cp ", out.json, " ~/software/browser/qATACviewer/frontend/src/data.json"))
      
      system("rm /bar/cfan/software/browser/qATACviewer/frontend/src/data.json")
      system("ln -s /scratch/qATACViewer_shared/qatacViewer.json /bar/cfan/software/browser/qATACviewer/frontend/src/data.json")
      system(paste0("cp ", out.json, " /scratch/qATACViewer_shared/qatacViewer.json"))
    }
  }

  invisible(jsongen)
}


# target.df <- data.frame(collection = "NK_fanc", pipedir = "~/4dn/nk/fanc/sth/target_pipe/",
#                         samples = c("Ly49A_neg_rep1", "Ly49A_neg_rep2", "Ly49A_neg_rep3",
#                                     "Ly49A_pos_rep1", "Ly49A_pos_rep2", "Ly49A_pos_rep3"))
#
# t <- qatacview.jsongen(target.df)
