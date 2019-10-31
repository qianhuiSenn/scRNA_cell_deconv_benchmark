#############################################################
#### example generating five-fold cross validation index ####
#############################################################
library(caret)
set.seed(3456)
##read in meta data###
meta_data <- readRDS(file = "/path/to/your/meta/file/")
#replace id_extract with correct column name##
class(meta_data$id_extract)
folds <- createFolds(meta_data$id_extract, 5)
str(folds)
split_up <- lapply(folds, function(ind, dat) dat[ind, ], dat = meta_data)
unlist(lapply(split_up, nrow))
table(split_up[[1]]$id_extract)
saveRDS(split_up, file = "sample_five_fold_list.rds")