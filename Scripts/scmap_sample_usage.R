#######################################
##### scID sample usage ######
#######################################
library(SingleCellExperiment)
library(scmap)
library(Matrix)
set.seed(3456)
DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"

# Load reference meta
reference_meta <- readRDS(file = file.path(DATA_DIR, 'ref_meta.rds'))
# load reference data
reference_norm <- readRDS(file = file.path(DATA_DIR, 'ref_norm.rds'))
# replace Group with the right column name in meta data
ref_ann <-  as.data.frame(reference_meta$Group)
colnames(ref_ann) <- "celltype"
# Create sce object
reference <- SingleCellExperiment(assays = list(normcounts = as.matrix(reference_norm)), colData = ref_ann)
logcounts(reference) <- normcounts(reference)
rowData(reference)$feature_symbol <- rownames(reference)
reference <- reference[!duplicated(rownames(reference)), ]
reference <- selectFeatures(reference, suppress_plot = FALSE)

#scmap-cell
reference <- indexCell(reference)

#query data
query_norm <- readRDS(file = file.path(DATA_DIR, "query_norm.rds"))
query <- SingleCellExperiment(assays = list(normcounts = as.matrix(query_norm)))
logcounts(query) <- normcounts(query)
rowData(query)$feature_symbol <- rownames(query)

# Prediction
scmapCell_results <- scmapCell(
  projection = query,
  list(
    ref = metadata(reference)$scmap_cell_index
  )
)

scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(reference)$celltype)
  ), threshold = -Inf
)

result <- scmapCell_clusters$scmap_cluster_labs
row.names(result) <-  colnames(scmapCell_results$ref$cells)
saveRDS(result, file = file.path(RESULT_DIR, 'scmap_prediction.rds'))