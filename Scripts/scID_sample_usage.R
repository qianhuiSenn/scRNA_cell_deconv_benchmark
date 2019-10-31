#######################################
##### scID sample usage ######
#######################################
library(scID)
library(MAST) #dependency of scID
library(biomod2) #dependency of scID if estimate weight from target
set.seed(3456)

DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"

reference_meta <- readRDS(file = file.path(DATA_DIR, 'ref_meta.rds'))

#require both ref and query to be normalized
query_norm <- readRDS(file = file.path(DATA_DIR, 'query_norm.rds'))
reference_norm <- readRDS(file = file.path(DATA_DIR, 'ref_norm.rds'))

## Replace Group with correct column name in meta
cell_type <- reference_meta$Group
names(cell_type) <- row.names(reference_meta)
reference_clusters <- cell_type

# Prediction
scID_output <- scid_multiclass(target_gem = query_norm, reference_gem = reference_norm, reference_clusters = reference_clusters, 
                               only_pos = FALSE, logFC = 0.5, estimate_weights_from_target = F, normalize_reference = F)
pred <- as.data.frame(scID_output$labels)
colnames(pred) <- 'scID_pred'

saveRDS(pred, file = file.path(RESULT_DIR, 'scID_5000_prediction.rds'))