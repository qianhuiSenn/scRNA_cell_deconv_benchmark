#### sample code for data preprocessing
#### example using Seurat
raw_count <- readRDS("path/to/data")
meta <- readRDS("path/to/data")
library(Seurat)
library(dplyr)
# 1. Cell and gene filtering (for both query and ref)
# 2. Normalization
data_tmp <- CreateSeuratObject(raw_count, min.cells = 3, min.features = 200, project = "example")
data_tmp <- NormalizeData(data_tmp, normalization.method = "LogNormalize", scale.factor=10000)
count <- as.matrix(GetAssayData(object = data_tmp, slot = "counts"))
norm <- as.matrix(GetAssayData(object = data_tmp))
saveRDS(count, "path/to/result/count.rds")
saveRDS(norm, "path/to/result/norm.rds")
# 3. Psudo-bulk reference matrix (for reference norm matrix)
norm_t <- t(norm)
merge_ref <- merge(norm_t, meta, by = 0)
merge_ref <- data.frame(merge_ref, row.names = 1, check.names = F, check.rows = F)
# replace orig_id with the column name of your meta data
bulk_ref <- aggregate(merge_ref[, 1:nrow(norm)], list(merge_ref$orig_id), mean)
bulk_ref <- data.frame(bulk_ref, row.names = 1, check.names = F, check.rows = F)
bulk_ref <- t(bulk_ref)
saveRDS(bulk_ref, "path/to/result/bulk_ref.rds")
# 4. Marker genes selection (from reference matrix)
data_tmp <- AddMetaData(data_tmp, meta)
#check the correct column name for the following command
Idents(data_tmp) <- "celltype"
data_tmp_m <- FindAllMarkers(data_tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data_tmp_marker <- data_tmp_m %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
saveRDS(data_tmp_marker, 'path/to/result/top_10_markers.rds')


# Batch correction
reference_count <- readRDS(file = file.path(DATA_DIR, 'ref_count.rds'))
query_count <- readRDS(file = file.path(DATA_DIR, "query_count.rds"))
# Use CCA for alignment
reference_seurat <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "batch")
query_seurat <- CreateSeuratObject(counts = query_count, min.cells = 0, min.features = 0, project = "batch")
#standard pipeline
reference_seurat <- NormalizeData(object = reference_seurat)
reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = 'vst', nfeatures = 10000)

query_seurat <- NormalizeData(object = query_seurat)
query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = 'vst', nfeatures = 10000)

##integration/batch
int_list <- list(reference_seurat,query_seurat)
int.anchors <- FindIntegrationAnchors(object.list = int_list, dims = 1:30, anchor.features = 10000)
int.integrated <- IntegrateData(anchorset = int.anchors, dims = 1:30)
DefaultAssay(int.integrated) <- "integrated"
data_tmp <- as.matrix(GetAssayData(int.integrated))

##separate data_tmp into batch corrected ref and query
ref_batch <- data_tmp[,which(colnames(data_tmp) %in% colnames(reference_count))]
qeury_batch <- data_tmp[,which(colnames(data_tmp) %in% colnames(query_count))]
