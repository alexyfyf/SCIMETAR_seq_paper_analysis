## refered to https://github.com/satijalab/azimuth-references/blob/master/human_motorcortex/scripts/export.R
## and https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format
 
library(Seurat)
library(SummarizedExperiment)

options(future.globals.maxSize = 1000 * 1024^2)

ref_se <- readRDS('../../R/data/dev_cell_2022/ref_se.rds')

library(Seurat)
library(Azimuth)

# 1. Convert SummarizedExperiment to Seurat (if not already done)
counts_matrix <- assay(ref_se, "counts")
meta_data <- as.data.frame(colData(ref_se))
ref_seurat <- CreateSeuratObject(counts = counts_matrix, meta.data = meta_data)

# 2. Pre-process the Reference
# Note: Azimuth references typically use SCTransform
ref_seurat <- SCTransform(ref_seurat, verbose = FALSE)
ref_seurat <- RunPCA(ref_seurat, verbose = FALSE)

# 3. Run UMAP and SAVE THE MODEL
# 'return.model = TRUE' is CRITICAL. It allows query cells to be projected onto this map.
ref_seurat <- RunUMAP(ref_seurat, dims = 1:30, reduction = "pca", return.model = TRUE)

# 4. Create the Azimuth Compatible Object
# We map your specific metadata columns to the 'metadata' argument
ref_seurat <- AzimuthReference(
  object = ref_seurat,
  refUMAP = "umap",       # Name of your UMAP reduction
  refDR = "pca",          # Name of your DR (we used "pca" above)
  refAssay = "SCT",       # The assay we processed
  dims = 1:50,            # Dimensions to use
  k.param = 31,           # Standard neighbor count
  metadata = c("cell", "immunophenotypes", "RNA_clusters", "developmental_stages"), # Your annotation columns
  reference.version = "1.0.0"
)

# 5. Define output directory
ref_dir <- "ref_azimuth/"
dir.create(ref_dir, showWarnings = FALSE)

# 6. Save the Annoy Index (Neighbor Graph)
# This extracts the calculated neighbors from the object and saves them
SaveAnnoyIndex(
  object = ref_seurat[["refdr.annoy.neighbors"]], 
  file = file.path(ref_dir, "idx.annoy")
)

# 7. Save the Reference Object
saveRDS(
  object = ref_seurat, 
  file = file.path(ref_dir, "ref.Rds")
)

print(paste("Reference saved to:", ref_dir))
