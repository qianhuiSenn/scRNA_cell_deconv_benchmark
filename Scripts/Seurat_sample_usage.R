################################
##### Seurat sample usage ######
################################

library(Seurat)
library(Matrix)
library(ggplot2)
library(knitr)
library(dplyr)
set.seed(3456)
DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"

#load reference
reference_count <- readRDS(file = file.path(DATA_DIR, 'ref_count.rds'))
reference_meta <- readRDS(file = file.path(DATA_DIR, 'ref_meta.rds'))
reference_seurat <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "example")
reference_seurat <- AddMetaData(reference_seurat, reference_meta)

#load query
query_count <- readRDS(file = file.path(DATA_DIR, "query_count.rds"))
query_seurat <- CreateSeuratObject(counts = query_count, min.cells = 0, min.features = 0, project = "example")

#standard pipeline
reference_seurat <- NormalizeData(object = reference_seurat)
reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = 'vst', nfeatures = 2000)
reference_seurat <- ScaleData(reference_seurat)

query_seurat <- NormalizeData(object = query_seurat)
query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = 'vst', nfeatures = 2000)
query_seurat <- ScaleData(query_seurat)

##prediction###
sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = query_seurat,
                                   dims = 1:30)
##replace Group with the actual column name from meta
predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$Group,
                            dims = 1:30)
query_seurat <- AddMetaData(object = query_seurat, metadata = predictions)

saveRDS(query_seurat$predicted.id, file = file.path(RESULT_DIR, 'seurat_prediction.rds'))
