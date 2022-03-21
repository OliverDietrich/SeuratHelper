################################################################################
# Helper functions for interoperability between Seurat and anndata

#' Read anndata object (.h5ad) and convert to Seurat object on disk
#' 
#' @param file Path to h5ad object
#' @returns Seurat object
#' @export
#' @example
#' ds <- read_h5ad("~/Downloads/pbmc3k.h5ad)
#' 
read_h5ad <- function(file) {
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
#' @example 
#' write_h5ad(ds, "~/Documents/pbmc3k.h5ad")
#' 
write_h5ad <- function(object, file) {
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

# end of document
################################################################################
