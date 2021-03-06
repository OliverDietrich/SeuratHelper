#' Plot low-dimensional cell embedding
#' 
#' This function creates a plot of a cell embedding stored in a Seurat object.
#' Cells can be colored by gene expression or metadata.
#' 
#' @param object Seurat object
#' @param embedding Name of the embedding (e.g. 'umap')
#' @param color Gene or metadata to color cells
#' @param label Type of group labels (text or label)
#' @param pointsize Point size
#' @param brush Brushed points (xmin, xmax, ymin, ymax)
#' @param cells Character of cell barcodes to plot
#' @param n.cells Number of cells to plot (randomly sampled)
#' @param assay Assay with expression data (e.g. RNA)
#' @param slot Slot of expression data (counts, data, scale.data)
#' @param dict Dictionary for the convert_features function
#' @param from From argument for convert_features
#' @param to To argument for convert_features
#' @param alpha Transparency of points (0-1)
#' @param color.transform Transformation of values used for color. 
#' Options include log, log2, log1p, log10
#' @param legend.position Position of the color legend.
#' @export
#' @examples 
#' #Plot CD3D on umap embedding in Seurat object
#' plot_embedding(object, "CD3D", "umap")
#' 
plot_embedding <- function(
  object     = NULL,
  color      = "",
  label      = FALSE,
  label.size = 4,
  embedding  = tail(names(object@reductions), 1),
  pointsize  = NULL,
  dim.1      = 1,
  dim.2      = 2,
  brush      = NULL,
  cells      = NULL,
  n.cells    = NULL,
  assay      = "RNA",
  slot       = "data",
  dict       = object@misc$features,
  from       = 1,
  to         = 2,
  alpha      = 1,
  shape      = 20,
  color.transform = "",
  legend.position = "right"
) {
  
  # Specify conditions that are required for the function to work
  stopifnot(
    class(object) == "Seurat",
    length(color) == 1,
    embedding %in% names(object@reductions),
    dim.1 %in% 1:dim(object@reductions[[embedding]]@cell.embeddings)[2],
    dim.2 %in% 1:dim(object@reductions[[embedding]]@cell.embeddings)[2]
  )
  if (is.null(embedding)) {
    stop("No low dimensional embeddings are stored in this Seurat object.")
  }
  
  # Create annotations based on input (some inputs are changed)
  title <- color
  
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
    pointsize <- dplyr::case_when(
      dplyr::between(dim(object)[2],     0,    250) ~ 3,
      dplyr::between(dim(object)[2],   250,   1000) ~ 2,
      dplyr::between(dim(object)[2],  1000,   5000) ~ 1,
      dplyr::between(dim(object)[2],  5000,  50000) ~ 0.1,
      dim(object)[2] > 50000 ~ 0.001
    ) # write function to smooth that
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
  
  # Transform values based on input
  if (class(df$col) %in% c("numeric", "integer")) {
    df$col <- switch (color.transform,
      "log"   = log(df$col),
      "log2"  = log2(df$col),
      "log1p" = log1p(df$col),                      
      "log10" = log10(df$col),
      df$col
    )
  }
  
  # Create color scale & guide elements
  if (class(df$col) %in% c("numeric", "integer")) {
    color_guide <- ggplot2::guide_colorbar(
      barwidth = 1, barheight = 15, ticks = FALSE, frame.colour = "black"
    )
    ann_cols <- viridis::scale_color_viridis(option = "B", direction = -1)
    groups <- FALSE
  } else {
    color_guide <- ggplot2::guide_legend(
      override.aes = list(size = 8), ncol = ceiling(length(unique(df$col))/10)
    )
    ann_cols <- NULL
    groups <- TRUE
  }
  
  # Create group labels
  if (groups & label %in% c("text", "label")) {
    ann <- dplyr::summarise(
      dplyr::group_by(df, col), x = median(x), y = median(y)
    )
    legend.position <- ""
    if (label == "text") {
      groups <- ggplot2::geom_text(
        ggplot2::aes(label = col, col = NULL), ann, size = label.size
        )
    } else {
      groups <- ggplot2::geom_label(
        ggplot2::aes(label = col), ann, size = label.size
        )
    }
  } else {
    groups <- NULL
  }
  
  # Order cells by color
  df <- df[order(df$col), ]
  
  ggplot2::ggplot(df, ggplot2::aes(x, y, color = col)) +
    ggplot2::geom_point(size = pointsize, alpha = alpha, shape = shape) +
    groups +
    ggplot2::coord_fixed() +
    ggplot2::labs(col = NULL, title = title) +
    ggplot2::guides(
      color = color_guide
    ) +
    ann_cols +
    ggplot2::theme_void(20) +
    ggplot2::theme(
      legend.position = legend.position
    ) +
    xlim + ylim
}