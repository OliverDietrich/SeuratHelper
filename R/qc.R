#' Add QC metrics
#' 
#' @param object SeuratObject
#' @returns SeuratObject
#' @export
#' 
add_qc_metrics <- function(object, assay = Seurat::DefaultAssay(object)) {
  
  stopifnot(
    class(object) == "Seurat"
  )
  
  # Library size
  object$libsize <- Matrix::colSums(object[[assay]]@counts)
  
  # Features
  object$ngenes <- Matrix::colSums(object[[assay]]@counts > 0)
  
  # Mitochondrial genes
  mt.genes <- rownames(object)[grep("^MT-", rownames(object))]
  object$percent.mt <- round(
    Matrix::colSums(object[[assay]][mt.genes, ]) / object$libsize, 3
  ) * 100
  
  # Ribosomal genes
  rp.genes <- rownames(object)[grep("^RP", rownames(object))]
  object$percent.rp <- round(
    Matrix::colSums(object[[assay]][rp.genes, ]) / object$libsize, 3
  ) * 100
  
  return(object)
}

#' Plot QC metrics
#' 
#' @param object SeuratObject
#' @param x Metadata slot for x-axis (default: libsize)
#' @param metrics Metadata slots for y-axis
#' @param split_by Metadat slot to split plot vertically
#' @param split_by_max Maximum number of allowed categories for 'split_by'
#' @param color_by Metadata slot to color points by
#' @param pt.bin Boolean, bin points to handle overplotting (uses bin2d)
#' @param nbins Number of bins used in bin2d
#' @param pt.size Point size
#' @param theme.size Theme size
#' 
#' @returns plot
#' @export
#' 
plot_qc_metrics <- function(object, x = "libsize", 
                            metrics = c("ngenes", "percent.mt", "percent.rp"),
                            split_by = NULL, split_by_max = 20,
                            color_by = NULL,
                            pt.bin = TRUE, nbins = 100,
                            pt.size = .5, theme.size = 15
                            ) {
  
  stopifnot(
    class(object) == "Seurat",
    x %in% names(object@meta.data),
    metrics %in% names(object@meta.data)
  )
  
  # Create data frame
  df <- data.frame(
    row.names = colnames(object),
    x = object@meta.data[[x]]
  )
  for (i in metrics) {
    df[[i]] <- object@meta.data[[i]]
  }
  if (is.null(split_by)) {
    df$split <- ""
  } else if (length(unique(object@meta.data[[split_by]])) <= split_by_max) {
    df$split <- factor(object@meta.data[[split_by]])
  } else if (length(unique(object@meta.data[[split_by]])) > split_by_max) {
    stop(paste(split_by, "has too many categories."))
  } else {
    stop(paste(split_by, "is not part of the objects metadata."))
  }
  
  # Bin points
  if (!is.null(color_by)) {
    df$col <- object@meta.data[[color_by]]
    main <- ggplot2::geom_point(ggplot2::aes(color = col), size = pt.size)
  } else if (pt.bin) {
    df$col <- 0
    main <- ggplot2::geom_bin2d(bins = nbins)
  } else {
    df$col<- ""
    main <- ggplot2::geom_point(size = pt.size)
  }
  
  # Colorscale
  if (class(df$col) %in% c("character", "factor", "boolean")) {
    colorscale <- ggplot2::scale_color_discrete()
    guides <- ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(
        size = 5
      ))
    )
  } else if (all(df$col == 0)) {
    colorscale <- viridis::scale_fill_viridis()
    guides <- ggplot2::guides(
      fill = ggplot2::guide_colorbar(barwidth=1, barheight=10, ticks=FALSE)
    )
  } else {
    colorscale <- viridis::scale_color_viridis()
    guides <- ggplot2::guides(
      color = ggplot2::guide_colorbar(barwidth=1, barheight=10, ticks=FALSE)
    )
  }
  
  # Collapse metrics
  df <- tidyr::gather(df, "metric", "value", -x, -split, -col)
  
  plot <- ggplot2::ggplot(df, ggplot2::aes(x, value)) +
    main +
    ggplot2::facet_grid(metric~split, scales = "free_y") +
    ggplot2::theme_classic(theme.size) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.major = ggplot2::element_line(),
      panel.grid.minor = ggplot2::element_line(),
      axis.text.x = ggplot2::element_text(angle=45,hjust=1,vjust=1)
    ) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_y_continuous(trans = "log10") +
    colorscale +
    guides +
    ggplot2::labs(col = color_by, x = x, y = NULL)
  
  return(plot)
}
