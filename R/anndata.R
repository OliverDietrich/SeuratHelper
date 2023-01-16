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
