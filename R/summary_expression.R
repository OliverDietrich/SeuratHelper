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
#' @param slot Name of slot to use
#' @param force.recalc Wheter to force recalculation
#' @returns SeuratObject with AUC added as 'assay'
#'
#' @export
#'
AddAUC <- function(
    object=NULL, 
    features=NULL,
    name=NULL,
    assay = NULL,
    slot = "data",
    force.recalc = FALSE
    ) {
  
  stopifnot(
    !is.null(object),
    class(name) == "character" & length(name) == 1,
    class(features) == "list"
  )
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  if (name %in% names(object@assays) & force.recalc == FALSE) {
    stop(paste("Assay", name, " already exists. 
               To replace set force.recalc = TRUE"))
  }
 
  ranks <- AUCell::AUCell_buildRankings(
    exprMat = slot(object@assays[[assay]], slot)
      )
  aucs <- AUCell::AUCell_calcAUC(features, ranks) # 400 - 7000
  aucs <- as.data.frame(t(aucs@assays@data$AUC))
  
  object[[name]] <- Seurat::CreateAssayObject(data = t(aucs))
  
  return(object)
}

#' Heatmap of gene expression
#' 
#' @param object Seurat object
#' @param features Character of genes to plot
#' @param coldata Annotation for columns
#' @param rowdata Annotation for rows
#' @param rowdata_label Name of the row annotation
#' @param collapse_replicates Whether to average over groups of cells
#' @param coldata_group_max Maxiumum number of categories for each 
#' column annotation
#' @param cells Cells to plot
#' @param heatmap_colors Name of color scale to use (default: "RdBu")
#' @param heatmap_color_dir Direction of heatmap color scale
#' @param annotation_colors List of color assignments for annotations
#' @param scale Whether to scale expression data
#' @param assay Which assay to take expression data from
#' @param slot Which slot of the assay to use ("data", "counts")
#' @param limits Limits of the color scale
#' @param title Plot title
#' @param show_rownames Whether to show rownames
#' @param show_colnames Whether to show colnames
#' @param ... Other parameters passed to pheatmap
#' @export
heatmap_expression <- function(
    object = NULL,
    features = NULL,
    coldata = NULL,
    rowdata = NULL,
    rowdata_label = "group",
    collapse_replicates = TRUE,
    coldata_group_max = 250,
    cells = colnames(object),
    heatmap_colors = "RdBu",
    heatmap_color_dir = -1,
    annotation_colors = NULL,
    scale = TRUE,
    assay = NULL,
    slot = "data",
    limits = c(-2, 2),
    title = "Differentially expressed genes",
    show_rownames = FALSE,
    show_colnames = FALSE,
    ...
  ) {
  
  stopifnot(
    !is.null(object),
    class(features) == "character"
  )
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  # Subset by cells ------------------------------------------------------------
  object <- subset(object, cells = cells)
  
  # Column annotations ---------------------------------------------------------
  if (is.null(coldata)) {
    message("Coldata not specified. Defaulting to Idents")
    coldata <- Seurat::Idents(object)
    cann <- data.frame(
      row.names = colnames(object),
      Idents = coldata
    )
  } else if (all(coldata %in% names(object@meta.data))) {
    if (length(coldata) == 1) { 
      # Select single column
      cann <- data.frame(
        row.names = colnames(object),
        label = object[[coldata]]
      )
      names(cann) <- coldata
    } else { 
      # Select multiple columns
      cann <- object@meta.data[, coldata]
    }
  } else {
    stop("Coldata has not been specified correctly. Exiting.")
  }
  
  # Row data & annotations -----------------------------------------------------
  
  index <- which(features %in% rownames(slot(object[[assay]], slot)))
  if (length(index) < length(features)) {
    leftout <- features[!features %in% features[index]]
    warning(paste("The following features were not found in the assay:",
                  paste(leftout, collapse = ", "), ". Will be removed..."))
  }
  if (length(index) == 0) {
    stop("No matching features specified.") 
  }
  
  if (is.null(rowdata)) {
    # No rowdata specified
    rann <- NA
    gaps_row <- NULL
  } else if (class(rowdata) %in% c("character", "factor") &
             length(rowdata)==length(features)) {
    # One vector of labels
    rann <- data.frame(
      row.names = make.unique(features[index], sep = "-"),
      label = factor(rowdata[index], unique(rowdata[index]))
    )
    names(rann) <- rowdata_label
    rann$duplicated <- as.character(duplicated(features[index]))
    # gaps_row <- head(as.numeric(cumsum(table(rann[[rowdata_label]]))), -1)
    gaps_row <- NULL
  } else {# TODO: add other possibilities (e.g. multiple annotations) #####
    warning(
    "Row annotations do not fit the supplied features and will be ignored..."
    )
    rann <- NA
    gaps_row <- NULL
  }
  
  # Summarize count data -------------------------------------------------------
  
  # Create groups
  ind <- unlist(lapply(lapply(cann, unique), length))
  cat_in <- coldata[ind <= coldata_group_max]
  cat_out <- coldata[ind > coldata_group_max]
  lvls <- row.names(unique_combinations(cann[, cat_in]))
  groups <- combinations(cann[, cat_in])
  lvls <- lvls[lvls %in% groups]
  groups <- factor(groups, lvls)
  
  
  # Summarize across replicates (e.g. cells)
  if (collapse_replicates) {
    mat <- summarize_groups(slot(object[[assay]], slot)[features[index], ], 
                            groups
    )
    smry <- unique_combinations(cann[, cat_in])
    cann$cells <- 1
    temp <- data.frame(group = groups, cells = 1)
    for (i in cat_out) {
      temp$num <- cann[[i]]
      x <- dplyr::summarise(dplyr::group_by(temp, group), num=mean(num))
      smry[[i]] <- x$num[match(rownames(smry), x$group)]
    }
    x <- dplyr::summarise(dplyr::group_by(temp, group), cells=sum(cells))
    smry[["log10_cells"]] <- log10(x$cells[match(rownames(smry), x$group)])
    ind <- rownames(smry) %in% colnames(mat)
    smry <- subset(smry, ind)
    cann <- smry
  } else {
    mat <- as.matrix(
      slot(object[[assay]], slot)[features[index], order(groups)]
    )
  }
  
  # Color encoding -------------------------------------------------------------
  
  # Scaling
  if (scale) {
    mat <- t(scale(t(mat)))
  }
  
  # Annotation colors
  if (is.null(annotation_colors)) {
    ann_colors <- NA
  } else if (class(colors) == "list") {
    # TODO: create input for colors
    ann_colors <- NA
  }
  
  # Breaks
  breaks <- seq(from = min(limits), to = max(limits), length.out = 100)
  
  # Color scale
  rcbrew <- RColorBrewer::brewer.pal.info
  if (class(heatmap_colors) != "character") {
    stop("Heatmap colors must be specified as a character vector.")
  } 
  if (length(heatmap_colors) > 1) {
    warning("Multiple heatmap colors specified. 
            Trying to generate colorRampPalette.")
  } else if (length(heatmap_colors) < 1) {
    warning("No heatmap colors have been specified. Reverting to default.")
    heatmap_colors <- "RdBu"
  } else if (heatmap_colors %in% rownames(rcbrew)) {
    print("Selecting color values from RColorBrewer...")
    heatmap_colors <- RColorBrewer::brewer.pal(n = 9, name = heatmap_colors)
  }
  if (heatmap_color_dir == -1) {
    heatmap_colors <- rev(heatmap_colors)
  }
  color_scale <- colorRampPalette(heatmap_colors)(length(breaks))
  
  # Plot -----------------------------------------------------------------------
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
    color  = color_scale, 
    gaps_col = head(as.numeric(cumsum(table(cann[,1]))), -1),
    gaps_row = gaps_row,
    main = title,
    silent = TRUE, 
    ...
    )
  
  return(plot)
}