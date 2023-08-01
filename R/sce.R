#' Convert from Seurat to SingleCellExperiment
#' 
#' This is a conversion function between R objects from class 'Seurat' to
#' 'SingleCellExperiment' to increase interoperability.
#' 
#' @param object Object of class 'Seurat'
#' @param assay Name of assay object to use for the main SingleCellExperiment
#' (all other assays will be stored as altExps, defaults to first assay in list)
#' @return Object of class 'SingleCellExperiment'
#' @export
#' 
to_sce <- function(
  object = NULL,
  assay = NULL
  ) {
  
  stopifnot(
    class(object) == "Seurat"
  )
  
  if (is.null(assay)) {
    assay <- names(object@assays)[1]
  }
  other.assays <- names(object@assays)[!names(object@assays) %in% assay]
  other.assays <- lapply(
    object@assays[other.assays], 
    function (x) {SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = x@counts, data = x@data),
    rowData = x@meta.features
  )})
  
  reductions <- lapply(object@reductions, function(x) {
    x@cell.embeddings
  })
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = slot(object@assays[[assay]], "counts"),
      data   = slot(object@assays[[assay]], "data")
    ), 
    colData = object@meta.data,
    rowData = object@assays[[assay]]@meta.features,
    metadata = object@misc, 
    altExps = other.assays
  )
  
  return(sce)
}


#' Convert from SingleCellExperiment to Seurat
#' 
from_sce <- function(object = NULL) {
  
  stopifnot(
    class(object) == "SingleCellExperiment"
  )
  
}