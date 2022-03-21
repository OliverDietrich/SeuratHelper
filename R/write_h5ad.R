#' Convert Seurat to h5ad
#' 
#' Save on disk Seurat object as anndata object (.h5ad)
#' 
#' @param object Seurat object
#' @param file Path to h5ad object (e.g. "~/Downloads/pbmc3k.h5ad")
#' @export
#' @examples
#' write_h5ad(ds, "~/Documents/pbmc3k.h5ad")
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
