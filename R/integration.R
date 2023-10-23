#' Perform Mutual Nearest Neighbor Integration
#' 
#' Runs the fastMNN implementation of the MNN algorithm from the batchelor 
#' package.
#' 
#' @param object Seurat object
#' @param features Set of features to use in MNN
#' @param name Name of the dimensional reduction object
#' @param transfer.reconstructed.assay Whether to transfer reconstructed counts
#' into a new assay
#' @param assay Assay of feature counts
#' @param slot Slot for feature counts
#' 
#' @export
#' 
RunMNN <- function(object=NULL,
                   batch=NULL,
                   features=NULL,
                   name="mnn",
                   transfer.reconstructed.assay = FALSE,
                   force.recalc = FALSE,
                   assay=DefaultAssay(object),
                   slot="data",
                   ...
                   ) {
  
  stopifnot(
    !is.null(object),
    batch %in% names(object@meta.data)
  )
  
  # Recalculate integration ----------------------------------------------------
  if (!is.null(object[[name]])) {
    if (force.recalc) {
      warning(paste("Embedding", name, "exists. But will be re-computed."))
    } else {
      stop(paste("Embedding", name, "exists. Aborting."))
    }
  }
  
  # Features -------------------------------------------------------------------
  if (is.null(features)) {
    features <- slot(object@assays[[assay]], "var.features")
  } else if (!all(features %in% rownames(object))) {
    stop("Not all features are present in the Seurat object. Aborting.")
  }
  if (is.logical(features)) {
    stop("No features have been selected. Run FindVariableGenes() or manually
         select features.")
  }
  
  # Run fastMNN ----------------------------------------------------------------
  result <- batchelor::fastMNN(slot(object@assays[[assay]], slot)[features, ], 
                     batch = object@meta.data[[batch]], ...)
  
  # Transfer embedding ---------------------------------------------------------
  object[[name]] <- Seurat::CreateDimReducObject(
    embeddings = result@int_colData$reducedDims$corrected, 
    loadings = result@rowRanges@elementMetadata$rotation,
    assay = assay, key = paste0(name, "_")
    )
  
  # Transfer assay -------------------------------------------------------------
  if (transfer.reconstructed.assay) {
    assay_key <- paste0(name, "_reconstructed")
    object[[assay_key]] <- Seurat::CreateAssayObject(
      data = result@assays@data$reconstructed
    )
  }
  
  return(object)
}