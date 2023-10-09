#' Annotate clusters based on maximum AUC score
#' 
#' @param object Seurat object 
#' @param group.by Vector in object's meta.data
#' @param features Features to consider
#' @param assay Assay containing AUCs
#' @param slot Slot to pull data from
#' 
#' @export
annotate_maxAUC <- function(object=NULL, assay = NULL, group.by = NULL,
                            slot = "data") {
  
  stopifnot(
    class(object) == "Seurat",
    group.by %in% names(object@meta.data),
    !is.null(assay)
  )
  
  if (is.null(group.by)) {
    group.by <- Seurat::Idents(object)
  }
  
  # Create labels
  mat <- summarize_groups(slot(object@assays[[assay]], slot), group.by)
  label <- which_rowMax(t(mat))
  names(label) <- levels(object$seurat_clusters)
  object$cell_type <- factor(label[as.character(object$seurat_clusters)])
  
  # Calculate certainty
  mat <- apply(mat, 2, sort, decreasing = TRUE)
  cert <- mat[1,] - mat[2,]
  names(cert) <- levels(object$seurat_clusters)
  object$cell_type_certainty <- cert[as.character(object$seurat_clusters)]
  
  return(object)
}