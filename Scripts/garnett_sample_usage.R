#######################################
##### garnett sample usage ######
#######################################
library(monocle)
library(garnett)
set.seed(3456)

DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"

#####training######
reference_count <- readRDS(file = file.path(DATA_DIR, 'ref_count.rds'))
gene_id <- row.names(reference_count)
fd <- data.frame(gene_id, gene_id, row.names = 1)
colnames(fd) <- c('gene_short_name')
pd <- readRDS(file = file.path(DATA_DIR, 'ref_meta.rds'))
reference_cds <- newCellDataSet(reference_count, phenoData = new("AnnotatedDataFrame", data = pd),
                                featureData = new("AnnotatedDataFrame", data = fd),
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
reference_cds <- estimateSizeFactors(reference_cds)

marker_file_path <- file.path(DATA_DIR, 'celltype_gene_marker.txt')

# suggestion: do marker check to remove or add gene
# marker_check <- check_markers(reference_cds, marker_file_path, db = 'none')
reference_classifier <- train_cell_classifier(cds = reference_cds,
                                              marker_file = marker_file_path,
                                              db='none',
                                              num_unknown = 50)
feature_genes <- get_feature_genes(reference_classifier)
print(head(feature_genes))

####Prediction####
query_count <- readRDS(file = file.path(DATA_DIR, 'query_count.rds'))
gene_id <- row.names(query_count)
fd <- data.frame(gene_id, gene_id, row.names = 1)
colnames(fd) <- c('gene_short_name')
# pd <- readRDS(file = file.path(DATA_DIR, 'query_meta.rds'))
query_cds <- newCellDataSet(query_count,
                            featureData = new("AnnotatedDataFrame", data = fd),
                            lowerDetectionLimit = 0.5,
                            expressionFamily = negbinomial.size())
query_cds <- estimateSizeFactors(query_cds)
query_cds <- classify_cells(query_cds, reference_classifier,
                            db = 'none',
                            cluster_extend = TRUE)
saveRDS(pData(query_cds)[,c('cluster_ext_type'),drop = F], file = file.path(RESULT_DIR, 'garnett_prediction.rds'))
