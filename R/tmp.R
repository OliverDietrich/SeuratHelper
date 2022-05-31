# TRIAL ZONE

path <- list(
  pbmc3k = "/home/oliver/Projects/Robert_Strasser/data/pbmc3k.h5ad",
  vieira = "/home/oliver/Projects/Robert_Strasser/data/vieira19.h5ad",
  lukassen = "/home/oliver/Projects/Robert_Strasser/data/lukassen20.h5ad",
  to = "/home/oliver/Downloads/object.h5ad"
)

from <- path$vieira
to <- "/home/oliver/Downloads/object.h5ad"

rhdf5::h5ls(from, all = T)

# Fetch coldata from h5ad
obs <- rhdf5::h5read(from, "obs", read.attributes = TRUE)
if (class(obs) == "list") {
  obs_attr <- attributes(obs)
  lvls <- rhdf5::h5read(from, "obs")[["__categories"]]
  obs_names <- rhdf5::h5read(from, "obs")[[obs_attr[["_index"]]]]
  obs <- as.data.frame(obs[obs_attr$`column-order`], row.names = obs_names
  )
  for (i in names(lvls)) {
    lv <- lvls[[i]]
    names(lv) <- sort(unique(obs[[i]])) + 1
    obs[[i]] <- factor(as.character(lv[obs[[i]] + 1]), levels = lvls[[i]])
  }
} else {
  if (class(obs) == "data.frame") {
    rownames(obs) <- obs$index
    obs <- obs[, -which(names(obs) == "index")]
  } else {
    stop("Unknown class of 'obs'")
  }
}


# Create new h5ad file
file.remove(to)
rhdf5::h5createFile(to)

# Add data.frame complete
lvls <- list()
for (i in names(obs)) {
  if (class(obs[[i]]) == "factor") {
    lvls[[i]] <- levels(obs[[i]])
    obs[[i]] <- as.character(obs[[i]])
  }
}
obs <- cbind(data.frame(index = row.names(obs)), obs)
rhdf5::h5write(obs, to, "obs")



################################################################################
# Convert obs data frame to list
nobs <- list(
  lapply(obs, levels)[!sapply(lapply(obs, levels), is.null)], obs_names
)
names(nobs) <- c("__categories", obs_attr[["_index"]])
obs <- as.list(obs)
for (i in names(obs)) {
  obs[[i]] <- as.double(obs[[i]]) -1
  attr(obs[[i]], "dim") <- length(obs[[i]])
}
obs_attr$`column-order` <- names(obs)
obs <- c(nobs, obs)
attributes(obs) <- obs_attr
