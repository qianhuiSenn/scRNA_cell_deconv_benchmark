################################
##### CHETAH sample usage ######
################################
library(CHETAH)
set.seed(3456)
DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"

# Load meta
ref_meta <- readRDS(file = file.path(DATA_DIR, 'ref_meta.rds'))

## Make SingleCellExperiments
### reference need to be normalized
reference_ct <- readRDS(file = file.path(DATA_DIR, 'ref_norm.rds'))
### replace Group with correct column name in meta
ref_ct <- as.character(ref_meta$Group)
reference <- SingleCellExperiment(assays = list(counts = reference_ct),
                                  colData = DataFrame(celltypes = ref_ct))

###optional normalized for query
query_ct <- readRDS(file = file.path(DATA_DIR, "query_count.rds"))
input <- SingleCellExperiment(assays = list(counts = query_ct))

## Run CHETAH, minimize unassign by set threshold as -Inf by thresh = -Inf, otherwise use default parameter
input <- CHETAHclassifier(input = input, ref_cells = reference)
## Extract celltypes
CHETAH_pred <- as.data.frame(input$celltype_CHETAH)
colnames(CHETAH_pred) <- 'CHETAH_pred'
saveRDS(CHETAH_pred, file = file.path(RESULT_DIR, 'CHETAH_prediction.rds'))
