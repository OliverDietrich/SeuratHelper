#' Test differential expression
#' 
#' This function performs differential expression analysis on a Seurat object
#' 
#' @param object Seurat object
#' @param groups Vector of group assignments
#' @param block Factor of group assignments to block
#' @param min.lfc Numeric scalar indicating the minimum difference in means to
#' calculate a p.value
#' @param assay Assay with expression data (e.g. RNA)
#' @param slot Slot of expression data (counts, data, scale.data)
#' 
test_diffexp <- function(
  object  = NULL,
  groups  = NULL,
  block   = NULL,
  min.lfc = 0.25,
  assay   = "RNA",
  slot    = "data",
  dict       = object@misc$features,
  from       = 1,
  to         = 2
) {
  
  # Check inputs
  stopifnot(
    class(object) == "Seurat",
    !is.null(groups),
    groups %in% names(object@meta.data)
  )
  
  # Convert group key to vector
  groups <- object@meta.data[[groups]]
  
  # Convert gene names to/from synonyms
  index <- rownames(slot(ds@assays[[assay]], slot))
  if (all(index %in% features[[from]])) {
    ids <- convert_names(index, dict, from, to)
  } else {
    ids <- index
  }
  
  # Compute group means
  group_means <- data.frame(row.names = index)
  for (i in levels(groups)) {
    j <- which(groups == i)
    group_means[[i]] <- Matrix::rowMeans(slot(ds@assays[[assay]], slot)[, j])
  }
  
  # Compute gene-level summary statistics
  markers <- data.frame(
    row.names   = index,
    name        = ids,
    max_logFC   = as.numeric(diff(apply(group_means, 1, range)))
  )
  
  index <- which(markers$fold_change > min.lfc)
  if (length(index) < 1) {
    warning("No genes exceed minimum log-fold change. Output empty.")
  }
  
  markers <- head(markers, 10) # REMOVE once ready
  
  return(markers)
}

if (FALSE) {
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Set values
  gene <- "MT-ND1"
  id <- convert_names(gene, features, 2, 1)
  assay <- "RNA"
  slot <- "data"
  n <- colnames(ds)
  p <- rownames(ds@assays$RNA@data)
  groups <- ds$seurat_clusters
  block <- ds$patient
  # block <- factor(sample(1:8, length(groups), replace = T))
  design <- data.frame(groups, block)
  sep <- "_of_"
  index <- split(n, design, sep = sep)
  group_centers <- data.frame(row.names = p)
  group_vars <- data.frame(row.names = p)
  FUN <- sparseMatrixStats::rowMedians
  FUN <- sparseMatrixStats::rowMeans2
  
  # Plot expression differences between groups by block
  df <- data.frame(
    x = design$groups,
    y = ds@assays$RNA@data[id, ],
    block = design$block
  )
  ggplot2::ggplot(df, ggplot2::aes(x,y, fill = block)) +
    ggplot2::geom_point(position = "jitter", size = 0.1) +
    ggplot2::geom_violin(
      fill = "antiquewhite", scale = "width", draw_quantiles = c(0.5)
    ) +
    ggplot2::geom_violin(scale = "width", draw_quantiles = c(0.5)) +
    ggplot2::labs(title = gene, x = NULL, y = NULL) +
    ggplot2::theme_classic(20)
  
  # Summarize across block (x ~ groups + block)
  for (i in names(index)) {
    if (length(index[[i]]) > 1)  {
      group_centers[[i]] <- FUN(ds@assays$RNA@data[, index[[i]] ])
    } else if (length(index[[i]]) == 1) {
      group_centers[[i]] <- ds@assays$RNA@data[, index[[i]] ]
    }
  }
  
  g <- names(group_centers)
  df <- as.data.frame(stringr::str_split(g, "_of_", simplify = TRUE))
  names(df) <- names(design)
  for (i in names(df)) {
    df[[i]] <- factor(df[[i]], levels(design[[i]]))
  }
  
  df$centers <- as.numeric(group_centers[id, ])
  
  ggplot2::ggplot(
    df, ggplot2::aes(x = groups, y = centers, col = block)
  ) + 
    ggplot2::geom_boxplot(col = "black", outlier.size = -1) +
    ggplot2::geom_point(position = "jitter") +
    ggplot2::labs(title = gene) +
    ggplot2::theme_classic(20)
  
  df$centers <- NULL
  
  # Summarize across groups
  group_centers <- t(
    apply(group_centers, 1, function (x) sapply(split(x, df$groups), median))
  )
  
  # Compute gene-level summary statistics
  markers <- data.frame(
    row.names   = p,
    name        = convert_names(p, features, 1, 2),
    max_logFC   = as.numeric(diff(t(sparseMatrixStats::rowRanges(group_centers)))),
    avg_logFC   = rowMeans(
      group_centers - sparseMatrixStats::rowMedians(group_centers)
    )
  )
  
  plot(markers$max_logFC, markers$avg_logFC)
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Groupwise testing
  
  groups <- ds$seurat_clusters # must be factor
  block <- ds$patient
  # block <- factor(sample(1:8, length(groups), replace = T))
  min_lfc <- 0.1
  
  n <- colnames(ds)
  p <- rownames(ds@assays$RNA@data)
  design <- data.frame(row.names = n, groups = groups)
  design$block <- block
  sep <- "_of_"
  
  assay <- "RNA"
  slot <- "data"
  
  # Summarize across groups and block (x ~ groups + block)
  index <- split(n, design, sep = sep)
  index <- index[which(sapply(index, length) > 0)]
  
  g <- design
  row.names(g) <- NULL
  g <- unique(g)
  g$index <- stringr::str_c(g$groups, g$block, sep = sep)
  
  group_centers <- data.frame(row.names = p)
  FUN <- sparseMatrixStats::rowMedians
  for (i in levels(g$groups)) {
    j <-g$index[g$groups == i]
    
    tmp <- data.frame(row.names = p)
    for (ii in j) {
      if (length(index[[ii]]) > 1) {
        tmp[[ii]] <- FUN(
          slot(ds@assays[[assay]], slot)[, index[[ii]] ]
        )
        group_centers[[i]] <- rowMedians(as.matrix(tmp))
      } else if (length(index[[ii]]) == 0) {
        group_centers[[i]] <- slot(ds@assays[[assay]], slot)[, index[[ii]] ]
      }
    }
  }
  
  markers <- data.frame(
    row.names = p,
    name      = convert_names(p, features, 1, 2),
    mean      = sparseMatrixStats::rowMeans2(slot(ds@assays[[assay]], slot)),
    max_logFC = as.numeric(diff(t(rowRanges(as.matrix(group_centers))))),
    cluster   = factor(max.col(group_centers), levels(group_design$groups))
  )
  ids <- rownames(markers[markers$max_logFC > min_lfc, ])
  
  group_centers <- list()
  FUN <- sparseMatrixStats::rowMedians
  for (i in levels(g$groups)) {
    j <-g$index[g$groups == i]
    
    group_centers[[i]] <- data.frame(row.names = ids)
    for (ii in j) {
      if (length(index[[ii]]) > 1) {
        group_centers[[i]][[ii]] <- FUN(
          slot(ds@assays[[assay]], slot)[ids, index[[ii]] ]
        )
      } else if (length(index[[ii]]) == 0) {
        group_centers[[i]] <- slot(ds@assays[[assay]], slot)[, index[[ii]] ]
      }
    }
  }
  
  # Plot
  markers$model <- predict(loess(max_logFC ~ mean, markers, span = 0.02))
  markers$corrFC <- markers$max_logFC - markers$model
  plot <- ggplot2::ggplot(
    data = markers[markers$max_logFC > 0.1, ], 
    ggplot2::aes(mean, max_logFC, col = cluster, label = name)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(y = model), col = "navy") +
    ggplot2::theme_classic(20) +
    ggrepel::geom_label_repel(data = markers[markers$corrFC > 1.9, ])
  
  plot
  
  plotly::ggplotly(plot)
  
  # Benchmark runtime
  {
    time_1 <- Sys.time()
    markers <- scran::findMarkers(ds@assays$RNA@data, groups = ds$seurat_clusters, block = ds$patient)
    time_2 <- Sys.time()
    time_2 - time_1
  }
  
}

