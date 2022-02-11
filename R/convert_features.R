#' Convert feature names
#' 
#' This is a conversion function for gene names. It takes a vector of gene names
#' and returns their corresponding synonyms (e.g. gene IDs) from a dictionary. 
#' 
#' @param names Character to convert
#' @param dict Data frame of gene synonyms (e.g. features.tsv)
#' @param from Column containing the names. Can be index (integer) 
#' or column name (character).
#' @param to Column containing the desired synonyms. Can be index (integer) 
#' or column name (character).
#' @return Character with synonyms
#' @export
#' @examples
#' # Convert human CD3D from 
#' convert_features("CD3D", features, 1, 2)
#' 
convert_names <- function(
  names = NULL,
  dict  = NULL,
  from  = 1,
  to    = 2
) {
  
  stopifnot(
    is.character(names),
    is.data.frame(dict),
    any(is.numeric(from), is.character(from)),
    any(is.numeric(to), is.character(to))
  )
  
  if (!any(from %in% 1:dim(dict)[2], from %in% names(dict))) {
    stop("Column 'from' not found. Please enter a valid column in dict.")
  }
  
  if (!any(to %in% 1:dim(dict)[2], to %in% names(dict))) {
    stop("Column 'to' not found. Please enter a valid column in dict.")
  }
  
  index <- na.omit(
    match(names, dict[[from]])
  )
  
  # Create warning if the output contains missing values
  if (length(attr(index, "na.action")) > 0) {
    warning("NAs detected. Not all input names have been found.")
  }
  
  synonyms <- dict[[to]][index]
  
  return(synonyms)
}