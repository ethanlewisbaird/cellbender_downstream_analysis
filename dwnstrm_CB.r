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
sample_id <- args[3]

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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

# Clustering analysis
obj <- CreateSeuratObject(counts = sparse.mat, min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 45000 & percent.mt < 5)
vln_plot <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = file.path(output_dir, paste0(sample_id, "_vln_plot.pdf")), plot = vln_plot)

cbres <- NormalizeData(obj)
cbres <- FindVariableFeatures(cbres, selection.method = "vst", nfeatures = 2000)
cbres <- ScaleData(cbres)
cbres <- RunPCA(cbres)
pca_plot <- DimPlot(cbres, reduction = "pca")
ggsave(filename = file.path(output_dir, paste0(sample_id, "_pca_plot.pdf")), plot = pca_plot)

ElbowPlot <- ElbowPlot(cbres)
ggsave(filename = file.path(output_dir, paste0(sample_id, "_elbow_plot.pdf")), plot = ElbowPlot)
cbres <- FindNeighbors(cbres, dims = 1:30)
cbres <- FindClusters(cbres, resolution = 0.5)
cbres <- RunUMAP(cbres, dims = 1:30)
umap_plot <- DimPlot(cbres, reduction = "umap")
ggsave(filename = file.path(output_dir, paste0(sample_id, "_umap_plot.pdf")), plot = umap_plot)

# Save the Seurat object
#saveRDS(cbres, file.path(output_dir, paste0(sample_id, "_seurat_obj.rds")))