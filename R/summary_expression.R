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
#' @param order_rows Whether to order rows by group
#' @param remove_duplicates Whether to remove duplicated rows
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
#' @param gaps_row Numeric vector specifying row gaps
#' @param ... Other parameters passed to pheatmap
#' @export
heatmap_expression <- function(
    object = NULL,
    features = NULL,
    coldata = NULL,
    rowdata = NULL,
    rowdata_label = "group",
    order_rows = FALSE,
    remove_duplicates = FALSE,
    collapse_replicates = TRUE,
    coldata_group_max = 100,
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
    gaps_row = NULL,
    ...
  ) {
  
  stopifnot(
    !is.null(object),
    class(features) == "character",
    cells %in% colnames(object),
    is.logical(order_rows),
    is.logical(remove_duplicates),
    is.logical(collapse_replicates),
    is.logical(scale)
  )
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  # Subset by cells ------------------------------------------------------------
  object <- subset(object, cells = cells)
  
  # Column annotations ---------------------------------------------------------
  if (is.null(coldata)) {
    message("Coldata not specified. Defaulting to Idents")
    coldata <- factor(Seurat::Idents(object))
    cann <- data.frame(
      row.names = colnames(object),
      Idents = coldata
    )
  } else if (all(coldata %in% names(object@meta.data))) {
    cann <- data.frame(row.names = colnames(object))
    for (i in coldata) {
      x <- object@meta.data[[i]]
      if (length(unique(x)) <= coldata_group_max) {
        cann[[i]] <- factor(x)
      } else {
        cann[[i]] <- x
      }
    }
  } else {
    stop("Coldata has not been specified correctly. Exiting.")
  }
  
  # Row data & annotations -----------------------------------------------------
  
  # Detect features in assay
  index <- which(features %in% rownames(slot(object[[assay]], slot)))
  if (length(index) == 0) {
    stop("No matching features specified.") 
  } else if (length(index) < length(features)) {
    leftout <- features[!features %in% features[index]]
    warning(paste("The following features were not found in the assay:",
                  paste(leftout, collapse = ", "), ". Will be removed..."))
  }
  
  if (is.null(rowdata)) {
    # No rowdata specified
    rann <- NA
    features <- features[index]
  } else if (class(rowdata) %in% c("character", "factor") &
             length(rowdata)==length(features)) {
    # One vector of labels
    features <- features[index]
    rowdata <- rowdata[index]
    # Vector levels
    if (class(rowdata) == "character") {
      rowdata <- factor(rowdata, unique(rowdata))
    } else {
      rowdata <- factor(rowdata)
    }
    # Create annotation data.frame
    if (remove_duplicates) {
      ind <- which(!duplicated(features))
      features <- features[ind]
      rann <- data.frame(row.names = features, label = rowdata[ind])
      names(rann) <- rowdata_label
    } else {
      rann <- data.frame(row.names = make.unique(features, sep = "-"), 
                         label = rowdata)
      names(rann) <- rowdata_label
      rann$duplicated <- factor(duplicated(features))
    }
  } else if (class(rowdata) == "data.frame" & 
             nrow(rowdata) == length(features)) {
    # Data frame with multiple annotations
    stop("Data.frame for row annotations is not supported yet!")
  } else {
    warning(
    "Row annotations do not fit the supplied features and will be ignored..."
    )
    rann <- NA
  }
  
  # Summarize count data -------------------------------------------------------
  
  # Create groups
  ind <- unlist(lapply(lapply(cann, unique), length))
  cat_in <- names(cann)[ind <= coldata_group_max]
  cat_out <- names(cann)[ind > coldata_group_max]
  lvls <- row.names(unique_combinations(subset(cann, select=cat_in)))
  groups <- combinations(subset(cann, select = cat_in))
  lvls <- lvls[lvls %in% groups]
  groups <- factor(groups, lvls)
  
  
  # Summarize across replicates (e.g. cells)
  if (collapse_replicates) {
    mat <- summarize_groups(slot(object[[assay]], slot)[features, ], 
                            groups
    )
    smry <- unique_combinations(subset(cann, select=cat_in))
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
      slot(object[[assay]], slot)[features, order(groups)]
    )
  }
  
  # Color encoding -------------------------------------------------------------
  
  # Scaling
  if (scale) {
    mat <- t(scale(t(mat)))
  }
  
  # Annotation colors
  ann_colors <- list()
  # Column annotations
  for (i in names(cann)) {
    x <- unique(cann[[i]])
    if (length(x) < 3) {
      ann_colors[[i]] <- c("white", "black", "grey")[1:length(x)]
      names(ann_colors[[i]]) <- levels(cann[[i]])
    } else if (length(x) <= 9) {
      ann_colors[[i]] <- RColorBrewer::brewer.pal(length(x), "Set1")
      names(ann_colors[[i]]) <- levels(cann[[i]])
    } else if (length(x) <= coldata_group_max) {
      ann_colors[[i]] <- NULL
    } else {
      ann_colors[[i]] <- RColorBrewer::brewer.pal(4, "Greens")
    }
  }
  # Row annotations
  for (i in names(rann)) {
    x <- unique(rann[[i]])
    if (length(x) < 3) {
      ann_colors[[i]] <- c("white", "black", "grey")[1:length(x)]
      names(ann_colors[[i]]) <- levels(rann[[i]])
    } else if (length(x) <= 9) {
      ann_colors[[i]] <- RColorBrewer::brewer.pal(length(x), "Set2")
      names(ann_colors[[i]]) <- levels(rann[[i]])
    } else if (length(x) <= coldata_group_max) {
      ann_colors[[i]] <- NULL
    } else {
      ann_colors[[i]] <- RColorBrewer::brewer.pal(4, "Greens")
    }
  }
  # Add manually specified colors
  if (is.list(annotation_colors)) {
    for (i in names(annotation_colors)) {
      if (i %in% names(ann_colors) & 
          all(names(annotation_colors[[i]]) %in% names(ann_colors[[i]]))) {
        ann_colors[[i]] <- annotation_colors[[i]]
      }
    }
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
  
  # Order ----------------------------------------------------------------------
  
  # Rows
  if (order_rows) {
    index <- order(rann[rowdata_label])
    mat <- mat[features[index],]
  }
  
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