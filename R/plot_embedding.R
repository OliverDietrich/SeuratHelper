#' Plot low-dimensional cell embedding
#' 
#' This function creates a plot of a cell embedding stored in a Seurat object.
#' Cells can be colored by gene expression or metadata.
#' 
#' @param object Seurat object
#' @param embedding Name of the embedding (e.g. 'umap')
#' @param color Gene or metadata to color cells
#' @param label Type of group labels (text or label)
#' @param pt.size Point size
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
#' 
plot_embedding <- function(
  object     = NULL,
  color      = "",
  split.by   = NULL, # TODO
  label      = FALSE,
  embedding  = tail(names(object@reductions), 1),
  pt.aggr    = TRUE,
  pt.aggr.breaks = 300,
  pt.size    = 1,
  pt.stroke  = .1,
  dim.1      = 1,
  dim.2      = 2,
  brush      = NULL,
  cells      = NULL,
  n.cells    = NULL,
  assay      = Seurat::DefaultAssay(object),
  slot       = "data",
  theme.size = 15, # TODO
  label.size = 6,
  split.by.rows = NULL,
  split.by.max = 20,
  dict       = object@misc$features,
  from       = 1,
  to         = 2,
  alpha      = 1,
  shape      = 20,
  color.transform = "",
  pl.title   = NULL,
  legend.position = "right",
  legend.rows = NULL
) {
  
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
  if (is.null(pl.title)) {
    title <- color
  } else {
    title <- pl.title
  }
  
  # Subset data ----------------------------------------------------------------
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
  
  # Convert gene names to/from synonyms ----------------------------------------
  if (!color %in% rownames(object) & color %in% dict[[to]]) {
    color <- convert_names(color, dict, to, from)
  }
  if (color %in% rownames(object) & is.null(object@assays[[assay]])) {
    stop(paste("No assay of name", assay, "present in this Seurat object."))
  }
  
  # Retrieve two dimensions from the embedding ---------------------------------
  rng <- c(dim.1, dim.2)
  df <- as.data.frame(object@reductions[[embedding]]@cell.embeddings[, rng])
  names(df) <- c("x", "y")
  
  # Color ----------------------------------------------------------------------
  if (color %in% names(object@meta.data)) {
    df$col <- object@meta.data[[color]]
  } else if (color %in% rownames(slot(object@assays[[assay]], slot))) {
    df$col <- slot(object@assays[[assay]], slot)[color, ]
  } else if (color %in% unlist(lapply(object@assays, rownames))) {
    assay <- stringr::str_remove_all(
      names(which(unlist(lapply(object@assays, rownames)) == color)), 
      "\\d"
      )
    message(paste("Expression/metadata found in different assay.",
                   "Switching to", assay))
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
  
  # Color scale & guides -------------------------------------------------------
  if (is.null(legend.rows)) {
    legend.rows <- ceiling(length(unique(df$col))/10)
  }
  if (class(df$col) %in% c("numeric", "integer", "array")) {
    color_guide <- ggplot2::guide_colorbar(
      barwidth = 1, barheight = 15, ticks = FALSE, frame.colour = "black"
    )
    ann_cols <- viridis::scale_color_viridis(option = "B", direction = -1)
    groups <- FALSE
  } else {
    color_guide <- ggplot2::guide_legend(
      override.aes = list(size = 8), ncol = legend.rows
    )
    ann_cols <- NULL
    groups <- TRUE
  }
  
  # Faceting -------------------------------------------------------------------
  if (!is.null(split.by)) {
    if (length(split.by) > 1) {
      message("More than one category present for split.by, taking the first.")
      split.by <- split.by[1]
    } 
    if (split.by %in% names(object@meta.data)) {
      n_split <- length(unique(object[[split.by]]))
      if (n_split > split.by.max) {
        stop(paste0("More categories for split.by than", split.by.max))
      }
      if (is.null(split.by.rows)) {
        split.by.rows <- round(sqrt(n_split))
      }
      df[["wrap"]] <- object@meta.data[, split.by]
      wrap <- ggplot2::facet_wrap(~wrap, nrow = split.by.rows)
    } else {
      message(paste0(
        "Plot will not be split.by '", split.by, 
        "' which is not present in meta.data"
      ))
      wrap <- NULL
    }
  } else {
    wrap <- NULL
  }
  
  # Create group labels --------------------------------------------------------
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
  
  # Order/summarize points -----------------------------------------------------
  if (pt.aggr) {
    # Distribute equally and determine color
    if (all(is.nan(df$col))) {
      df$col <- 1
      FUN <- sum
      title <- "Density (cells/dot)"
    } else if (class(df$col) %in% c("numeric", "integer", "array")) {
      FUN <- mean
    } else {
      FUN <- select_most_frequent_category
    }
    df <- summarize_overlapping_rows(df, breaks = pt.aggr.breaks, FUN=FUN)
  } else {
    # Order by value
    df <- df[order(df$col), ]
  }
  
  # Plot -----------------------------------------------------------------------
  plot <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = col)) +
    ggplot2::geom_point(size = pt.size, alpha = alpha, shape = shape,
                        stroke = pt.stroke) +
    groups +
    ggplot2::coord_fixed() +
    ggplot2::labs(col = NULL, title = title) +
    ggplot2::guides(
      color = color_guide
    ) +
    ann_cols +
    ggplot2::theme_void(theme.size) +
    ggplot2::theme(
      legend.position = legend.position,
      title = ggplot2::element_text(vjust = .5),
      panel.border = ggplot2::element_rect(size = .5, fill=NA)
    ) +
    xlim + ylim +
    wrap
  
  return(plot)
}

#' Plot cell type markers on embedding
#' 
#' @param object SeuratObject
#' @param markers Character vector of marker genes
#' @param embedding Embedding name (slot in reductions)
#' @param nrow Number of rows
#' @param pt.size Point size
#' 
#' @returns plot
#' @export
#'
plot_markers_embedding <- function(object, markers=NULL,
                                   embedding=tail(names(object@reductions), 1),
                                   pt.aggr    = TRUE,
                                   pt.aggr.breaks = 300,
                                   pt.size    = .5,
                                   pt.stroke  = .1,
                                   pt.shape=16,
                                   nrow=floor(sqrt(length(markers))),
                                   assay = NULL,
                                   slot = "data",
                                   scale = TRUE,
                                   markers.max = 100,
                                   pl.title = NULL,
                                   col.title = NULL
                                   ) {
  
  stopifnot(
    class(object) == "Seurat"
  )
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  if (is.null(markers)) {
    warning("No markers specified. Defaulting to whole assay.")
    markers <-rownames(object@assays[[assay]])
    if (length(markers) > markers.max) {
      stop(paste("To many features in assay:", assay))
    }
  }
  
  # Fetch data -----------------------------------------------------------------
  df <- data.frame(
    x = object@reductions[[embedding]]@cell.embeddings[, 1],
    y = object@reductions[[embedding]]@cell.embeddings[, 2]
  )
  for (i in markers) {
    if (i %in% rownames(slot(object[[assay]], slot))) {
      df[[i]] <- slot(object[[assay]], slot)[i, ]
    }
  }
  
  # Re-shape -------------------------------------------------------------------
  df <- tidyr::gather(df, "gene", "count", -x, -y)
  df$gene <- factor(df$gene, unique(df$gene))
  
  # Order/summarize points -----------------------------------------------------
  if (pt.aggr) {
    df$x <- round(df$x, 2)
    df$y <- round(df$y, 2)
    df <- dplyr::summarise(dplyr::group_by(df, gene, x, y), count = mean(count))
  }
  
  # Scale ----------------------------------------------------------------------
  if (scale) {
    df <- dplyr::group_by(df, gene)
    df <- dplyr::mutate(df, count = scale(count)[,1])
    color_scale <- ggplot2::scale_color_distiller(palette = "RdBu")
    plot_title <- "Normalized and scaled gene expression"
    col_title <- "z-score"
  } else {
    plot_title <- "Normalized gene expression"
    col_title <- "logcount"
    color_scale <- viridis::scale_color_viridis(option = "A", direction = -1)
  }
  
  if (!is.null(pl.title)) {
    plot_title <- pl.title
  }
  if (!is.null(col.title)) {
    col_title <- col.title
  }
  
  # Define range of scale ------------------------------------------------------
  rng <- c(-2, 2)
  df$count[df$count < min(rng)] <- min(rng)
  df$count[df$count > max(rng)] <- max(rng)
  
  # Plot
  plot <- ggplot2::ggplot(df, ggplot2::aes(x, y, col = count)) +
    ggplot2::geom_point(size = pt.size, stroke = pt.stroke, shape = pt.shape) +
    ggplot2::facet_wrap(~gene, nrow = nrow) +
    color_scale +
    ggplot2::theme_void(20) +
    ggplot2::coord_fixed() +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, ticks = F)
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(size = .5, fill=NA)
    ) +
    ggplot2::labs(title = plot_title,
                  subtitle = paste("Embedding:", embedding), 
                  col = col_title)
  
  return(plot)
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
