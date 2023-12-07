#' Summarize groups of columns in a matrix
#'
#' @param x A matrix-like object
#' @param groups A factor of group labels equal to ncol(x)
#' @param FUN A function to summarize columns
#' (default value requires sparseMatrixStats)
#' @export
#' 
summarize_groups <- function (
    x = NULL, groups = NULL, FUN = sparseMatrixStats::rowMeans2
) {
  
  stopifnot(
    any(class(x) %in% c("matrix", "dgCMatrix", "dgTMatrix")),
    length(groups) == ncol(x)
  )
  
  if (class(groups) != "factor") {
    warning(paste("Groups are of type", class(groups), 
                  "and will be converted to factor."))
    groups <- factor(groups)
  }
  
  rn <- rownames(x)
  index <- split(1:length(groups), groups)
  
  fun <- function(j) {
    if (length(j) > 1) {
      FUN(x[, j])
    } else {
      x[, j]
    }
  }
  matrix <- sapply(index, fun)
  
  rownames(matrix) <- rn
  
  return(matrix)
}

#' Paste columns of a data.frame into a vector
#'
#' @param x Data.frame
#' @param sep Separator put in between content of different columns
#' @export
#'
combinations <- function(x, sep = "-") {
  result <- apply(x, 1, function(x) paste0(x, collapse = sep))
  result <- factor(result)
  return(result)
}

#' Create unique combinations from a data.frame of factors
#'
#' @param x Data.frame of factors for which to create unique combinations.
#' Alternatively, a list of factor levels.
#' @export
#'
unique_combinations <- function(x, sep = "-") {
  
  stopifnot(class(x) %in% c("data.frame", "list"))
  
  x <- lapply(x, factor)
  x <- lapply(x, levels)
  
  m <- matrix(
    nrow = length(x)+1, ncol = 3,
    dimnames = list(0:length(x), c("x", "y", "z"))
  )
  m["0", ] <- c(1, prod(sapply(x, length)), 1)
  cmb <- data.frame(row.names = 1:prod(sapply(x, length)))
  
  for (i in 1:length(x)) {
    j <- as.character(i-1)
    n <- as.character(i)
    
    m[n, "x"] <- length(x[[i]])
    m[n, "y"] <- m[j, "y"] / m[n, "x"]
    m[n, "z"] <- m[j, "x"] * m[j, "z"]
    
    index <- names(x)[[i]]
    cmb[[index]] <- factor(
      rep(x[[i]], time = m[n, "z"], each = m[n, "y"]), x[[index]]
    )
    
  }
  
  rownames(cmb) <- combinations(cmb)
  
  return(cmb)
}

#' Extract the column names of row maxima
#'
#' @param x Matrix containing numeric values
#' @export
#'
which_rowMax <- function(x) {
  v <- rep(NaN, nrow(x))
  is.max <- x == matrixStats::rowMaxs(x)
  index <- which(!apply(is.max, 1, all))
  v[index] <- apply(is.max[index, ], 1, which)
  
  return(colnames(x)[v])
}

#' Select maximum occuring category in vector
#' 
#' @param v Vector
#' 
#' @export
select_most_frequent_category <- function(v) {
  index <- table(v)
  sample(names(index[index == max(index)]), 1)
}

#' Summarize rows that have overlapping coordinates
#' 
#' @param df Data.frame
#' @param x Column key
#' @param y Column key
#' @param FUN Function to summarize by
#' @param decimals Integer of decimals to round to
#' 
#' @returns Data.frame
#' 
summarize_overlapping_rows <- function(df=NULL, x="x", y="y", FUN=mean,
                                       breaks = 1000) {
  
  stopifnot(
    class(df) == "data.frame",
    class(df[[x]]) == "numeric",
    class(df[[y]]) == "numeric"
  )
  
  # Remove NAs
  df <- df[!is.na(df$col), ]
  
  # Round coordinates
  for (i in c(x, y)) {
    df[[i]] <- as.numeric(cut(df[[i]], breaks))
  }
  df$match <- paste(df$x, df$y, sep = ":")
  
  # Account for split_by column
  if (is.null(df[["wrap"]])) {
    df$wrap <- "1"
  }
  
  # Summarize
  df <- dplyr::summarise(
    dplyr::group_by(df, match, wrap),
    x = mean(x), y = mean(y), col = FUN(col)
  )
  
  # Make data.frame
  df <- as.data.frame(df)
  
  # Return
  return(df)
}