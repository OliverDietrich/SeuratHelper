#' Plot low-dimensional cell embedding
#' 
#' This function creates a plot of a cell embedding stored in a Seurat object.
#' Cells can be colored by gene expression or metadata.
#' 
#' @param object Seurat object
#' @param embedding Name of the embedding (e.g. 'umap')
#' @param color Gene or metadata to color cells
#' @param pointsize Point size
#' @param dict Dictionary for the convert_features function
#' @param from From argument for convert_features
#' @param to To argument for convert_features
#' @param brush Brushed points (xmin, xmax, ymin, ymax)
#' @param cells Character of cell barcodes to plot
#' @param n.cells Number of cells to plot (randomly sampled)
#' @param assay Assay with expression data (e.g. RNA)
#' @param slot Slot of expression data (counts, data, scale.data)
#' @export
#' @examples 
#' #Plot CD3D on umap embedding in Seurat object
#' plot_embedding(object, "CD3D", "umap")
#' 
plot_embedding <- function(
  object    = NULL,
  color     = "",
  embedding = tail(names(object@reductions), 1),
  pointsize = NULL,
  dim.1     = 1,
  dim.2     = 2,
  dict      = object@misc$features,
  from      = 1,
  to        = 2,
  brush     = NULL,
  cells     = NULL,
  n.cells   = NULL,
  assay     = "RNA",
  slot      = "data"
) {
  
  # Specify conditions that are required for the function to work
  stopifnot(
    class(object) == "Seurat",
    embedding %in% names(object@reductions),
    dim.1 %in% 1:dim(object@reductions[[embedding]]@cell.embeddings)[2],
    dim.2 %in% 1:dim(object@reductions[[embedding]]@cell.embeddings)[2]
  )
  if (is.null(embedding)) {
    stop("No low dimensional embeddings are stored in this Seurat object.")
  }
  
  # Subset data based on cells
  if (!is.null(cells)) {
    object <- subset(object, cells = cells)
  } else if (!is.null(n.cells)) {
    object <- subset(object, cells = sample(colnames(object), n.cells))
  }
  if (is.null(brush)) {
    xlim <- NULL
    ylim <- NULL
  } else {
    xlim <- ggplot2::xlim(brush$xmin, brush$xmax)
    ylim <- ggplot2::ylim(brush$ymin, brush$ymax)
  }
  
  # Compute values for input arguments that are NULL
  if (is.null(pointsize)) {
    pointsize <- 3 / log10( dim(object)[2] - 0.5)
  }
  # Convert gene names to/from synonyms
  if (!color %in% rownames(object) & color %in% dict[[to]]) {
    color <- convert_names(color, dict, to, from)
  }
  if (color %in% rownames(object) & is.null(object@assays[[assay]])) {
    stop(paste("No assay of name", assay, "present in this Seurat object."))
  }
  
  # Retrieve two dimensions from the embedding
  rng <- c(dim.1, dim.2)
  df <- as.data.frame(object@reductions[[embedding]]@cell.embeddings[, rng])
  names(df) <- c("x", "y")
  
  # Retrieve expression/metadata for color
  if (color %in% names(object@meta.data)) {
    df$col <- object@meta.data[[color]]
  } else if (color %in% rownames(object)) {
    df$col <- slot(object@assays[[assay]], slot)[color, ]
  } else {
    warning(paste0("Expression/metadata for '", color, "' not found."))
    df$col <- NaN
  }
  
  # Order cells by color
  df <- df[order(df$col), ]
  
  ggplot2::ggplot(df, ggplot2::aes(x, y, color = col)) +
    ggplot2::geom_point(size = pointsize) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(20) +
    xlim + ylim
}