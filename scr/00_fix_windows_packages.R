############################################################
# Windows Fix: Seurat VlnPlot compatibility
# Run ONCE in clean R session
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(pkgType = "binary")

# Close RStudio before running if possible

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", type = "binary")
}

lib <- .libPaths()[1]

# Remove leftover lock directories
locks <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
if (length(locks) > 0) {
  unlink(locks, recursive = TRUE, force = TRUE)
  message("Removed lock folders.")
}

# Remove ggplot2 if installed
if ("ggplot2" %in% rownames(installed.packages())) {
  remove.packages("ggplot2")
}

# Install stable version
remotes::install_version(
  "ggplot2",
  version = "3.5.2",
  repos = "https://cloud.r-project.org",
  upgrade = "never"
)

install.packages(c("patchwork", "scales", "gtable"), type = "binary")

message("✅ Fix completed. Restart RStudio.")
