#######################################
##### SCINA sample usage ######
#######################################
library(SCINA)
library(dplyr)
DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"
set.seed(3456)

true_marker <- readRDS(file = file.path(DATA_DIR, 'marker_top_10.rds'))
##########Table format##########
## cluster   gene  ###
## Group1   GeneA ###
## Group1   GeneB ###
## ......  ...... ###
#####################
##SUGGESTION FOR THE MARKER GENE LIST MAKING
## MODIFIY THIS IF YOU HAVE MORE THAN 5 GROUPS OR DIFFERENT CELL TYPE NAMES
marker_gene_list <- list(
  Group1 = true_marker[which(true_marker$cluster == "Group1"), ]$gene,
  Group2 = true_marker[which(true_marker$cluster == "Group2"), ]$gene,
  Group3 = true_marker[which(true_marker$cluster == "Group3"), ]$gene,
  Group4 = true_marker[which(true_marker$cluster == "Group4"), ]$gene,
  Group5 = true_marker[which(true_marker$cluster == "Group5"), ]$gene
)

###Alternative use
group_list <- unique(true_marker$cluster)
cell_num <- length(group_list)
marker_gene_list <- list()
for (i in 1:cell_num){
  cell <- group_list[i]
  gene <- marker_list[which(marker_list$cluster == cell), ]$gene
  marker_gene_list[[i]] <- gene
}
names(marker_gene_list) <- group_list

###load query###
query <- readRDS(file = file.path(DATA_DIR,'query_norm.rds'))

#predicition
#set allow_unknown = F if you are certain the cell type in query
results = SCINA(query, marker_gene_list, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=T, log_file='SCINA.log')

SCINA_pred <- data.frame(cell_id = colnames(query), SCINA_pred = results$cell_labels, row.names = 1)
SCINA_pred$SCINA_pred <- as.character(SCINA_pred$SCINA_pred)
saveRDS(SCINA_pred,  file = file.path(RESULT_DIR, 'SCINA_prediction.rds'))
