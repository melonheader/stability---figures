### INFO: Miscellaneous functions
### DATE: 25.09.2019
### AUTHOR: Artem Baranovskii

# environment
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# force function to run quietly (Hadley Wickam)
# ---------------------------------------------------------------------------- # 
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# lists
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# omit NA lists
# ---------------------------------------------------------------------------- #
na_omit_list <- function(x) { 
  return(x[!sapply(x, function(y) all(is.na(y)))]) 
}
# omit empty lists
# ---------------------------------------------------------------------------- #
empt_omit_list <- function(x) {
  return(x[lapply(x, length) > 0])
}

# matrices
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# get lower triangle of the correlation matrix
# ---------------------------------------------------------------------------- #
get_lower_tri<-function(cormtx){
  cormtx[upper.tri(cormtx)] <- NA
  return(cormtx)
}
# get upper triangle of the correlation matrix
# ---------------------------------------------------------------------------- #
get_upper_tri <- function(cormtx){
  cormtx[lower.tri(cormtx)]<- NA
  return(cormtx)
}
# reorder cor matrix according to hc
# ---------------------------------------------------------------------------- #
reorder_cormtx <- function(cormtx){
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormtx) / 2)
  hc <- hclust(dd)
  cormtx <- cormtx[hc$order, hc$order]
  return(cormtx)
}


# vectors
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
str_reverse <- function(x) {
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
}
  


# tibbles
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# find columns containing a string 
str_detect_cols <- function(table, pattern) {
  names_col <- apply(table, 2, function(x) any(grepl(paste(pattern, collapse = "|"), x, ignore.case = TRUE)) == TRUE) %>%
    which(. == TRUE) %>%
    names(.)
  return(names_col)
}


# percentile vars
# ----------------------------------------------------------------------------------------------- #
fact_percentile <- function(x, y, seq_by = 0.1, fact = FALSE, report_splits = F) {
  #
  if (missing(y)) {
    if (length(x) > 1) {
      y <- quantile(x, probs = seq(0, 1, seq_by))
    } else {
      stop("Provide vector of values with length(x) > 1 or Provide precomputed percentiles (y)")
    }
  }
  #
  looping <- function(x, y) {
    for (i in 2:(length(y) - 1)) {
      if (i == 2) {
        out <- ifelse(x <= y[[i]], 
                      ifelse(report_splits, 
                             paste0(names(y[i]), "; x <= ", round(y[[i]], 2)), 
                             paste0(names(y[i]))), 
                      "Nein")
        names(out) <- names(y[i])
        if (out != "Nein") {
          return(out)
        }
      }
      if (i < (length(y) - 1)) {
        out <- ifelse(x > y[[i]] & x <= y[[i + 1]], 
                      ifelse(report_splits,
                             paste0(names(y[i + 1]), "; ", round(y[[i]], 2), " < x <= ", round(y[[i + 1]], 2)),
                             paste0(names(y[i + 1]))), 
                      ifelse(x > y[[length(y)]], 
                             "Above 100th percentile", 
                             "Nein")
                      )
        names(out) <- names(y[i + 1])
        if (out != "Nein") {
          return(out)
        }
      }
      if (i == (length(y) - 1)) {
        out <- ifelse(x > y[[i]] & x <= y[[i + 1]], 
                      ifelse(report_splits, 
                             paste0(names(y[i + 1]), "; ", "x > ", round(y[[i]], 2)),
                             paste0(names(y[i + 1]))), 
                      ifelse(x > y[[length(y)]], 
                             "Above 100th percentile",
                             "Nein")
                      )
        names(out) <- names(y[i + 1])
        if (out != "Nein") {
          return(out)
        }
      }
    }
  }
  #
  out <- purrr::map(x, 
                    ~ looping(.x, y)) %>% 
    unlist()
  #
  if (fact == TRUE) {
    out <- factor(out, levels = out[names(y[2:length(y)])])
  }
  #
  return(out)
}

# join tables by partial string matching
# ---------------------------------------------------------------------------- #
partial_join <- function(x, y, by_x, pattern_y) {
  require(dplyr, quietly = T, warn.conflicts = F)
  by_x <- enquo(by_x)
  pattern_y <- enquo(pattern_y)
  d <- full_join(mutate(x, dummy = 1), 
                 mutate(y, dummy = 1), by = "dummy")
  d <- d %>% mutate(., take = stringr::str_detect(!!by_x, !!pattern_y))
  d %>% 
    filter(take == T) %>% 
    dplyr::select(-dummy, -take)
}
# spread table keeping values of key in  the names of new columns
# ---------------------------------------------------------------------------- #
myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}
# expand grid keeping only unique pairs
# ---------------------------------------------------------------------------- #
expand_grid_unique <- function(x, y, include.equals=FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i - include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level = 0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}
# normalize numerical vector to 0 to 1 scale
# ---------------------------------------------------------------------------- #
norm_to_one <- function(x, na.rm = FALSE) {
  x_norm <- (x - min(x, na.rm = na.rm)) / (max(x, na.rm = na.rm) - min(x, na.rm = na.rm))
  return(x_norm)
}
# 1-based make.unique()
# ---------------------------------------------------------------------------- #
make.unique.1b = function(x, sep = '.'){
  ave(x, x, FUN = function(a) { if (length(a) > 1) { paste(a, 1:length(a), sep = sep) } else {a} })
}


# memory functions
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- # 
# calculates the memory usage of objects (Vedran Franke)
# ---------------------------------------------------------------------------- #
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing = FALSE, head = FALSE, n = 5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthand (Vedran Franke)
# ---------------------------------------------------------------------------- #
lsos <- function(..., n = 10) {
  rep_tab <- .ls.objects(..., order.by = "Size", decreasing = TRUE, head = TRUE, n = n)
  rep_tab$Size <- paste0(round(rep_tab$Size / 1024^2, 2), " Mb")
  return(rep_tab)
}
# loooped gc (Vedran Franke)
# ---------------------------------------------------------------------------- #
cleanm <- function(n = 10) { 
  for (i in 1:n) {
    gc(verbose = TRUE) 
  }
}