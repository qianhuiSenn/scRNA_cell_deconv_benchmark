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
reference_cds <- newCellDataSet(reference_count,
                                featureData = new("AnnotatedDataFrame", data = fd),
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
reference_cds <- estimateSizeFactors(reference_cds)

marker_file_path <- file.path(DATA_DIR, 'celltype_gene_marker.txt') #if available

##alternatively if not##
##########DE Table format#########
## cluster   gene  ###
## Group1   GeneA ###
## Group1   GeneB ###
## ......  ...... ###
#####################
#####how to make a marker txt from a output of DE table#####
marker <- readRDS(file = "path/to/your/DE/table.rds")
fileConn <- file("celltype_gene_marker.txt")
group_list <- unique(marker$cluster)
cell_num <- length(group_list)
temp <- NULL
for(i in 1:cell_num){
  gene_list <- marker[which(marker$cluster == group_list[i]), ]$gene
  anno1 <- paste('>', group_list[i], sep = '')
  anno2 <- paste("expressed: ", paste(gene_list, collapse=","), sep = '')
  temp <- c(temp, anno1, anno2)
}
writeLines(temp, fileConn)
close(fileConn)
# suggestion: do marker check to remove or add gene
marker_check <- check_markers(reference_cds, marker_file_path, db = 'none')
gene_exclude <- marker_check[which(marker_check$summary != "Ok"), ]$gene_id
marker_after_exclude <- marker[which(!(marker$gene %in% gene_exclude)), ]

fileConn <- file("abc.txt")
group_list <- unique(marker_after_exclude$cluster)
cell_num <- length(group_list)
temp <- NULL
for(i in 1:cell_num){
  gene_list <- marker_after_exclude[which(marker_after_exclude$cluster == group_list[i]), ]$gene
  anno1 <- paste('>', group_list[i], sep = '')
  anno2 <- paste("expressed: ", paste(gene_list, collapse=","), sep = '')
  temp <- c(temp, anno1, anno2)
}
writeLines(temp, fileConn)
close(fileConn)

###make classifier
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
