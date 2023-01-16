#' Add QC metrics
#' 
#' @param object SeuratObject
#' 
#' @returns SeuratObject
#' 
add_qc_metrics <- function(object) {
  
  stopifnot(
    class(object) == "Seurat"
  )
  
  # Library size
  object$libsize <- ds$nCount_RNA
  
  # Features
  object$ngenes <- object$nFeature_RNA
  
  # Percent MT
  mt.genes <- rownames(object)[grep("^MT-", rownames(object))]
  object$percent.mt <- round(
    Matrix::colSums(object@assays$RNA@counts[mt.genes, ]) /
      object$libsize, 3
  ) * 100
  
  return(object)
}

#' Plot QC metrics
#' 
#' @param object SeuratObject
#' 
#' @returns plot
#' 
plot_qc_metrics <- function(object, libsize = "libsize", 
                            metrics = c("ngenes", "percent.mt"),
                            quality = NULL, color_by = NULL
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
    df$qc <- ds@meta.data[[quality]]
  }
  
  if (!is.null(color_by)) {
    df$color_by <- ds@meta.data[[color_by]]
    main <- ggplot2::geom_point(ggplot2::aes(color = color_by), size = .25)
  } else {
    main <- ggplot2::geom_bin2d(bins = 200)
  }
  
  # Collapse metrics
  df <- tidyr::gather(df, "metric", "value", -x, -qc, -color_by)
  
  plot <- ggplot2::ggplot(df, ggplot2::aes(x, value)) +
    main +
    ggplot2::facet_grid(metric~qc, scales = "free") +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_y_continuous(trans = "log10") +
    viridis::scale_fill_viridis()
  
  return(plot)
}
