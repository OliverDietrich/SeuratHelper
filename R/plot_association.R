#' Plot association between cell data
#' 
#' This function creates a plot showing the association of elements of metadata
#' or expression of cells stored in a Seurat object.
#' 
#' @param object Seurat object
#' @export
#' @examples
#' # Plot library size against number of detected features
#' plot_association(object, "nCount_RNA", "nFeature_RNA")
#' 
plot_association <- function(
  object = NULL
) {
  
}