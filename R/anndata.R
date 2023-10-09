#' Get layer from AnnData
#' 
#' @param file Path to h5ad file of type AnnData
#' @param assay Assay name (default: X)
#' 
#' @returns Sparse matrix
#' 
#' @export
read_layer_h5ad <- function(file = NULL, name = "X") {
  
  stopifnot(
    endsWith(file, "h5ad"),
    name %in% rhdf5::h5ls(file)$name
  )
  
  if (name != "X") {
    name <- paste0("layers/", name)
  }
  
  X <- Matrix::sparseMatrix(
    i = rhdf5::h5read(file, paste0(name, "/indices")) + 1,
    p = rhdf5::h5read(file, paste0(name, "/indptr")),
    x = as.integer(rhdf5::h5read(file, paste0(name, "/data"))),
    dims = c(
      length(rhdf5::h5read(file, "var/_index")), 
      length(rhdf5::h5read(file, "obs/_index"))
    )
  )
  
  return(X)
}

#' Get slot from AnnData
#' 
#' @param file Path to h5ad file of type AnnData
#' @param slot Slot name (default: obs)
#' 
#' @return Data.frame
#' 
#' @export
read_slot_h5ad <- function(
    file = NULL,
    slot = "obs"
) {
  
  stopifnot(
    endsWith(file, "h5ad"),
    slot %in% rhdf5::h5ls(file)$name
  )
  
  index <- match(slot, rhdf5::h5ls(file)$name)
  path <- paste0(rhdf5::h5ls(file)$group[index], "/", slot)
  data <- rhdf5::h5read(file, path, read.attributes = TRUE)
  
  # Convert lists into factors
  for (i in names(data)) {
    if (class(data[[i]]) == "list") {
      print(paste("Converting", i, "to factor"))
      index <- data[[i]]$codes + 2
      lvls <- c(NA, data[[i]]$categories)
      names(lvls) <- sort(unique(index))
      data[[i]] <- factor(lvls[index], lvls)
    }
  }
  
  # Convert to encoding type
  df <- data.frame(
    row.names = data[[attr(data, "_index")]]
  )
  for (i in attr(data, "column-order")) {
    df[[i]] <- data[[i]]
  }
  
  return(df)
}

#' Read obsm from AnnData
#' 
#' @param file Path to h5ad file of type AnnData
#' 
#' @returns List of matrices
#' 
#' @export
read_obsm_h5ad <- function(file) {
  
  stopifnot(
    endsWith(file, "h5ad")
  )
  
  data <- rhdf5::h5read(file, "obsm", read.attributes = TRUE)
  
  obsm <- list()
  for (i in names(data)) {
    
    if (all(class(data[[i]]) == "list")) {
      obsm[[i]] <- data.frame(
        row.names = data[[i]][["_index"]]
      )
      for (j in names(data[[i]])) {
        obsm[[i]][[j]] <- data[[i]][[j]]
      }
      obsm[[i]][["_index"]] <- NULL
    }
    
    if (any(class(data[[i]]) == "matrix")) {
      obsm[[i]] <- t(data[[i]])
    }
  }
  attributes(obsm) <- attributes(data)
  
  return(obsm)
}

#' Read anndata object (.h5ad) and convert to Seurat object on disk
#' 
#' @param file Path to h5ad object
#' @returns Seurat object
#' @export
#' 
read_h5ad <- function(
  file = NULL
  ) {
  
  if (!stringr::str_detect(file, "h5ad")) stop("Please supply object.h5ad") 
  
  print("Begin conversion")
  
  # Retrieve row data (var)
  var <- as.data.frame(
    dplyr::bind_cols(rhdf5::h5read(file, "var"))
  )
  rownames(var) <- var[["_index"]]
  var <- var[, -1]
  
  # Retrieve col data (obs)
  obs <- as.data.frame(
    dplyr::bind_cols(rhdf5::h5read(file, "obs")[-1]) # fix level import
  )
  rownames(obs) <- obs[["_index"]]
  obs <- obs[, -1]
  
  # Retrieve (normalized) count matrix
  data <- as(rhdf5::h5read(file, "X"), "dgTMatrix")
  rownames(data) <- rownames(var)
  colnames(data) <- rownames(obs)
  
  print("Create Seurat object")
  # Create Seurat object
  ds <- Seurat::CreateSeuratObject(
    counts = data, data = data, meta.data = obs
  )
  
  # Add meta features (vars)
  ds@assays$RNA@meta.features <- var
  
  # Add dimensional reduction objects (obsm)
  for (i in names(rhdf5::h5read(file, "obsm"))) {
    print(i)
    emb <- t(rhdf5::h5read(file, "obsm")[[i]])
    rownames(emb) <- rownames(obs)
    ds@reductions[[i]] <- Seurat::CreateDimReducObject(
      embeddings = emb, key = i
    )
  }
  
  print("Done.")
  return(ds)
}

#' Convert Seurat to h5ad
#' 
#' Save on disk Seurat object as anndata object (.h5ad)
#' 
#' @param object Seurat object
#' @param file Path to h5ad object (e.g. "~/Downloads/pbmc3k.h5ad")
#' @export
#' 
write_h5ad <- function(
  object = NULL, 
  file   = NULL
  ) {
  
  if (class(object) != "Seurat") stop("Please supply object.h5ad") 
  if (!stringr::str_detect(file, "h5ad")) stop("Please supply object.h5ad") 
  
  # Create file
  if (file.exists(file)) stop("File already exists.")
  rhdf5::h5createFile(file)
  
  print("Begin conversion")
  
  # Add count matrix
  rhdf5::h5write(
    obj = Matrix::as.matrix(ds@assays$RNA@data), file = file, name = "X"
  )
  
  # Add coldata (obs)
  df <- ds@meta.data
  df <- cbind(data.frame("gene_name" = rownames(df)), df)
  rhdf5::h5write(obj = df, file = file, name = "obs")
  
  # Add meta features (var)
  df <- ds@assays$RNA@meta.features
  df <- cbind(data.frame("gene_name" = rownames(df)), df)
  rhdf5::h5write(obj = df, file = file, name = "var")
  
  # Add dimensional reductions (obsm)
  rhdf5::h5createGroup(file, "obsm")
  for (i in names(ds@reductions)) {
    j <- paste0("obsm/", i)
    rhdf5::h5write(
      obj = t(ds@reductions[[i]]@cell.embeddings), file = file, name = j
    )
  }
  print("Done.")
}
