#### sample code for data preprocessing
#### example using Seurat
raw_count <- readRDS("path/to/data")
meta <- readRDS("path/to/data")
library(Seurat)
library(dplyr)
# 1. Cell and gene filtering
# 2. Normalization
data_tmp <- CreateSeuratObject(raw_count, min.cells = 3, min.features = 200, project = "example")
data_tmp <- NormalizeData(data_tmp, normalization.method = "LogNormalize", scale.factor=10000)
count <- as.matrix(GetAssayData(object = data_tmp, slot = "counts"))
norm <- as.matrix(GetAssayData(object = data_tmp))
saveRDS(count, "path/to/result/count.rds")
saveRDS(norm, "path/to/result/norm.rds")
# 3. Psudo-bulk reference matrix
norm_t <- t(norm)
merge_ref <- merge(norm_t, meta, by = 0)
merge_ref <- data.frame(merge_ref, row.names = 1, check.names = F, check.rows = F)
# replace orig_id with the column name of your meta data
bulk_ref <- aggregate(merge_ref[, 1:nrow(norm)], list(merge_ref$orig_id), mean)
bulk_ref <- data.frame(bulk_ref, row.names = 1, check.names = F, check.rows = F)
bulk_ref <- t(bulk_ref)
saveRDS(bulk_ref, "path/to/result/bulk_ref.rds")
# 4. Marker genes selection
data_tmp <- AddMetaData(data_tmp, meta)
#check the correct column name for the following command
Idents(data_tmp) <- "celltype"
data_tmp_m <- FindAllMarkers(data_tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data_tmp_marker <- data_tmp_m %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
saveRDS(data_tmp_marker, 'path/to/result/top_10_markers.rds')




