#####Here is the illustration of pulling results from cross-dataset experiments######

RESULT_DIR <- '/PATH/TO/RESULT'
dataset_list <- c('pbmc', 'pancreas', 'tabula_lung', 'tabula_full', 'simulation')

for (dataset in dataset_list){
  if (dataset == 'pbmc'){
    orig_id <- readRDS('/path/to/pbmc/query/meta/pbmc3k_meta.rds')
    orig_id <- orig_id[,'orig_id', drop = F]
    orig_id <- orig_id[order(row.names(orig_id)), , drop = F]
    DATA_DIR <- '/path/to/prediction/from/pbmc'
  }else if(dataset == 'pancreas'){
    orig_id <- readRDS('/path/to/pancreas/query/meta/pancreas_fluidigmc1_metadata.rds')
    orig_id <- orig_id[,'celltype', drop = F]
    orig_id <- orig_id[order(row.names(orig_id)), , drop = F]
    DATA_DIR <- '/path/to/prediction/from/pancreas'
  }else if(dataset == 'tabula_lung'){
    orig_id <- readRDS('/path/to/tabula_lung/query/meta/droplet_lung_meta.rds')
    orig_id <- orig_id[,'cell_ontology_class', drop = F]
    orig_id <- orig_id[order(row.names(orig_id)), , drop = F]
    DATA_DIR <- '/path/to/prediction/from/tabula_muris_lung'
  }else if(dataset == 'tabula_full'){
    orig_id <- read.csv('/path/to/tabula_full/query/meta/drop_meta_subset.csv', row.names = 1)
    orig_id <- orig_id[,'cell_ontology_class', drop = F]
    orig_id <- orig_id[order(row.names(orig_id)), , drop = F]
    DATA_DIR <- '/path/to/prediction/from/tabula_muris_full'
  }else{
    orig_id <- read.csv('/path/to/simulation/query/meta/cell_info_2.csv', row.names = 1)
    orig_id <- orig_id[,'Group', drop = F]
    orig_id <- orig_id[order(row.names(orig_id)), , drop = F]
    DATA_DIR <- '/path/to/prediction/from/Simulation'
  }
  if(dataset != 'tabula_full'){
    method_list <- c('RPC', 'CP', 'singleCellNet', 'scID', 'SCINA', 'CHETAH', 'seurat', 'singler', 'scmap', 'garnett')
    for (method in method_list){
      prediction <- readRDS(file = file.path(DATA_DIR, paste(dataset, '_', method, '_', 'prediction.rds', sep = '')))
      if(method == 'seurat'){
        prediction <- as.data.frame(prediction)
        colnames(prediction) <- 'seurat_pred'
      }
      if(method == 'singler'){
        colnames(prediction) <- 'singler_pred'
      }
      if(method == 'scmap'){
        colnames(prediction) <- 'scmap_pred'
      }
      if(method == 'garnett'){
        colnames(prediction) <- 'garnett_pred'
      }
      prediction <- prediction[order(row.names(prediction)), , drop = F]
      if(method != 'garnett'){
        orig_id <- cbind(orig_id, prediction)
      }else{
        cross_result <- merge(orig_id, prediction, by = 0, all.x = T, sort = F)
        write.csv(cross_result, file = file.path(RESULT_DIR, paste(dataset, 'result.csv', sep = '_')))
      }
    }
  }else{
    ###garnett doesn't produce results for tabula_full experiment
    method_list <- c('RPC', 'CP', 'singleCellNet', 'scID', 'SCINA', 'CHETAH', 'seurat', 'singler', 'scmap')
    for (method in method_list){
      prediction <- readRDS(file = file.path(DATA_DIR, paste('tabula', '_', method, '_', 'prediction.rds', sep = '')))
      if(method == 'seurat'){
        prediction <- as.data.frame(prediction)
        colnames(prediction) <- 'seurat_pred'
      }
      if(method == 'singler'){
        colnames(prediction) <- 'singler_pred'
      }
      if(method == 'scmap'){
        colnames(prediction) <- 'scmap_pred'
      }
      prediction <- prediction[order(row.names(prediction)), , drop = F]
      if(method != 'scmap'){
        orig_id <- cbind(orig_id, prediction)
      }else{
        cross_result <- merge(orig_id, prediction, by = 0, all.x = T, sort = F)
        write.csv(cross_result, file = file.path(RESULT_DIR, paste(dataset, 'result.csv', sep = '_')))
      }
    }
  }
}

#####Then move to metrics_calculation_example.ipynb#####