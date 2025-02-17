.libPaths(c("/data/ebaird/R-packages", .libPaths()))
library(Seurat)
library(ggplot2)
#library(rhdf5)
library(Matrix)
library(hdf5r)
library(sceasy)
library(reticulate)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set the output PDF file
pdf(file.path(output_dir, "Rplots.pdf"))

# Load CellBender output .h5 file and convert to seurat object
Read_CellBender_h5_Mat <- function(
  file_name,
  use.names = TRUE,
  unique.features = TRUE
) {
  # Check hdf5r installed
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    cli_abort(message = c("Please install hdf5r to read HDF5 files",
                          "i" = "`install.packages('hdf5r')`")
    )
  }
  # Check file
  if (!file.exists(file_name)) {
    stop("File not found")
  }

  if (use.names) {
    feature_slot <- 'features/name'
  } else {
    feature_slot <- 'features/id'
  }

  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")

  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]


  sparse.mat <- sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )

  if (unique.features) {
    features <- make.unique(names = features)
  }

  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")

  infile$close_all()

  return(sparse.mat)
}

sparse.mat <- Read_CellBender_h5_Mat(input_file, use.names=TRUE, unique.features=TRUE)

#clustering analysis
obj <- CreateSeuratObject(counts = sparse.mat, min.cells = 3, min.features = 200)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
cbres <- NormalizeData(obj)
cbres <- FindVariableFeatures(cbres, selection.method = "vst", nfeatures = 2000)
cbres <- ScaleData(cbres)
cbres <- FindVariableFeatures(cbres)
cbres <- RunPCA(cbres)
cbres <- FindNeighbors(cbres, dims = 1:10)
cbres <- FindClusters(cbres, resolution = 0.5)
cbres <- RunUMAP(cbres, dims = 1:9)
DimPlot(cbres, reduction = "umap")
#cbres <- FindNeighbors(cbres, k.param = 50)

pca_plot <- DimPlot(cbres, reduction = "pca")
umap_plot <- DimPlot(cbres, reduction = "umap")

# Print the plots
print(pca_plot)
print(umap_plot)

# Close the PDF device
dev.off()