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
  
  # Percent MT
  mt.genes <- rownames(object)[grep("^MT-", rownames(object))]
  object$percent.mt <- round(
    Matrix::colSums(object[[assay]][mt.genes, ]) /
      object$libsize, 3
  ) * 100
  
  return(object)
}

#' Plot QC metrics
#' 
#' @param object SeuratObject
#' @returns plot
#' @export
#' 
plot_qc_metrics <- function(object, libsize = "libsize", 
                            metrics = c("ngenes", "percent.mt"),
                            quality = NULL, color_by = NULL,
                            pt.bin = TRUE, nbins = 100,
                            pt.size = .5, theme.size = 15
                            ) {
  
  stopifnot(
    class(object) == "Seurat",
    libsize %in% names(object@meta.data),
    metrics %in% names(object@meta.data)
  )
  
  # Create data frame
  df <- data.frame(
    row.names = colnames(object),
    x = object@meta.data[[libsize]]
  )
  for (i in metrics) {
    df[[i]] <- object@meta.data[[i]]
  }
  if (is.null(quality)) {
    df$qc <- ""
  } else {
    df$qc <- object@meta.data[[quality]]
  }
  
  # Bin points
  if (!is.null(color_by)) {
    df$col <- object@meta.data[[color_by]]
    main <- ggplot2::geom_point(ggplot2::aes(color = col), size = pt.size)
  } else if (pt.bin) {
    df$col <- ""
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
  } else {
    colorscale <- viridis::scale_color_viridis()
    guides <- ggplot2::guides(
      color = ggplot2::guide_colorbar(barwidth=1, barheight=10, ticks=FALSE)
    )
  }
  
  # Collapse metrics
  df <- tidyr::gather(df, "metric", "value", -x, -qc, -col)
  
  plot <- ggplot2::ggplot(df, ggplot2::aes(x, value)) +
    main +
    ggplot2::facet_grid(metric~qc, scales = "free") +
    ggplot2::theme_classic(theme.size) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.major = ggplot2::element_line(),
      panel.grid.minor = ggplot2::element_line(),
    ) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_y_continuous(trans = "log10") +
    colorscale +
    guides +
    ggplot2::labs(col = color_by)
  
  return(plot)
}
