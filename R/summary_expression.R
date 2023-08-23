#' Calculate AUC for marker list
#' 
#' Calculate the area under the curve (AUC) for each gene set for each 
#' observation using the AUCell package
#' 
#' @param object Seurat object
#' @param features A list of vectors of features for expression programs; 
#' each entry should be a vector of feature names
#' @param name Name of the feature set (e.g. celltype) used as assay name
#' @param assay Name of assay to use
#' @param slot
#' @returns SeuratObject with AUC added as 'assay'
#'
AddAUC <- function(
    object=NULL, 
    features=NULL,
    name=NULL,
    assay = NULL,
    slot = "data",
    force = FALSE
    ) {
  
  stopifnot(
    !is.null(object),
    class(name) == "character" & length(name) == 1,
    class(features) == "list"
  )
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  if (name %in% names(ds@assays) & force == FALSE) {
    stop(paste("Assay", name, " already exists. To replace set force = TRUE"))
  }
 
  ranks <- AUCell::AUCell_buildRankings(
    exprMat = slot(object@assays[[assay]], slot)
      )
  aucs <- AUCell::AUCell_calcAUC(gsets, ranks) # 400 - 7000
  aucs <- as.data.frame(t(aucs@assays@data$AUC))
  
  object[[name]] <- Seurat::CreateAssayObject(data = t(aucs))
  
  return(object)
}

#' Heatmap of gene expression
#' 
#' @param object Seurat object
#' @param features Character of genes to plot
#' @param coldata 
#' @param rowdata
#' @param collapse_replicates Whether to average over groups of cells
#' @param title Plot title
#' @param limits Limits of the color scale
#' 
#' @export
#' 
heatmap_expression <- function(
    object = NULL,
    features = NULL,
    coldata = NULL,
    rowdata = NULL,
    rowdata_label = "cluster",
    collapse_replicates = TRUE,
    colors = NULL,
    scale = TRUE,
    assay = NULL,
    slot = "data",
    limits = c(-2, 2),
    title = "Differentially expressed genes",
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 10,
    fontsize_row = fontsize,
    fontsize_col = fontsize
  ) {
  
  stopifnot(
    !is.null(object),
    class(features) == "character"
  )
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  if (is.null(coldata)) {
    message("Coldata not specified. Defaulting to Idents")
    coldata <- Seurat::Idents(ds)
    cann <- data.frame(
      row.names = colnames(object),
      Idents = coldata
    )
  } else if (all(coldata %in% names(object@meta.data))) {
    cann <- object@meta.data[, coldata]
    coldata <- combinations(cann)
  } else {
    stop("Coldata has not been specified correctly. Exiting.")
  }
  
  index <- which(features %in% rownames(slot(object[[assay]], slot)))
  if (length(index) < length(features)) {
    leftout <- features[!features %in% features[index]]
    warning(paste("The following features were not found in the assay:",
                  leftout, "and will be ignored."))
  }
  if (length(index) == 0) {
   stop("No matching features specified.") 
  }
  
  if (is.null(rowdata)) {
    rann <- NA
    gaps_row <- NA
  } else if (class(rowdata) %in% c("character", "factor") &
             length(rowdata)==length(features)) {
    rann <- data.frame(
      row.names = features[index],
      label = rowdata
    )
    names(rann) <- rowdata_label
    gaps_row <- head(as.numeric(cumsum(table(rann[,1]))), -1)
  } else {
    warning(
    "Row annotations do not fit the supplied features and will be ignored..."
    )
    rann <- NA
    gaps_row <- NULL
  }
  
  if (collapse_replicates) {
    cann <- unique_combinations(cann)
    mat <- summarize_groups(slot(object[[assay]], slot)[features[index], ], 
                            coldata)
  } else {
    mat <- slot(object[[assay]], slot)[features[index], order(coldata)]
  }
  
  if (scale) {
    mat <- t(scale(t(mat)))
  }
  
  if (is.null(colors)) {
    ann_colors <- NA
  } else if (class(colors) == "list") {
    # TODO: create input for colors
  }
  
  breaks <- seq(from = min(limits), to = max(limits), length.out = 100)
  
  plot <- pheatmap::pheatmap(
    mat,
    cluster_cols = FALSE,                              
    cluster_rows = FALSE, 
    breaks = breaks,
    annotation_col = cann,
    annotation_row = rann,
    annotation_colors = ann_colors,
    show_colnames = show_colnames,
    show_rownames = show_rownames,
    color  = colorRampPalette(
      rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(length(breaks)), 
    gaps_col = head(as.numeric(cumsum(table(cann[,1]))), -1),
    gaps_row = gaps_row,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    main = title,
    silent = TRUE
    )
  
  return(plot)
}