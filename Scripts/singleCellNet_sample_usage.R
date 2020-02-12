#######################################
##### singleCellNet sample usage ######
#######################################
library(singleCellNet)
library(dplyr)
set.seed(3456)

DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"


###For smaller dataset: cell type smaller than 4 cannot be split in singleCellNet and need to be exclude out#####
splitCommon<-function(sampTab, ncells, dLevel="cell_ontology_class"){
  cts<-unique(as.vector(sampTab[,dLevel]))
  trainingids<-vector()
  for(ct in cts){
    stX<-sampTab[sampTab[,dLevel]==ct, ,drop = F]
    if(nrow(stX) > 4){
      cat(ct, ": ")
      ccount<-nrow(stX)-3
      ccount<-min(ccount, ncells)
      cat(nrow(stX),"\n")
      trainingids<-append(trainingids, sample(rownames(stX), ccount))
    }
  }
  val_ids<-setdiff(rownames(sampTab), trainingids)
  list(train=sampTab[trainingids, , drop = F], val=sampTab[val_ids, , drop = F])
}
####################################################################################################################

#load query
query_count <- readRDS(file = file.path(DATA_DIR, 'query_count.rds'))
query_meta <- readRDS(file = file.path(DATA_DIR, 'query_meta.rds'))
gene_qeury <- rownames(query_count)

#load referecne/train
reference_count <- readRDS(file = file.path(DATA_DIR, 'ref_count.rds'))
reference_meta <- readRDS(file = file.path(DATA_DIR, 'ref_meta.rds'))
reference_meta <- droplevels(reference_meta)

#process ref
commonGenes <- intersect(rownames(reference_count), gene_qeury)
reference_count <- reference_count[commonGenes, ]

#Split for training and assessment, and transform training data
#dLevel need to set to the correct column name for meta data
reference_list <- splitCommon(reference_meta, ncells = 100, dLevel = "Group")
reference_train <- reference_list[[1]]
expTrain <- reference_count[, rownames(reference_train)]
tmpX <- weighted_down(expTrain, 1.5e3, dThresh = 0.25)
expTrain <- trans_prop(tmpX, 1e4)

#find best set of classifier genes
#replace Group by the correct column name for meta data
cgenes2 <- findClassyGenes(expTrain, reference_train, "Group", topX=10)
cgenesA<-cgenes2[['cgenes']]
grps<-cgenes2[['grps']]
length(cgenesA)

#find the best pairs
expT <- as.matrix(expTrain[cgenesA, ])
dim(expT)
xpairs<-ptGetTop(expT, grps, topX=25, sliceSize=5000)
length(xpairs)
#TSP transform the training data
pdTrain<-query_transform(expT[cgenesA, ], xpairs)
#train the classif
rf_tspAll<-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps, nRand=10, ntrees=1000)

#apply to query data
query_TransAll<-query_transform(query_count[cgenesA,], xpairs)
nqRand<- 10
query_classification <-rf_classPredict(rf_tspAll, query_TransAll, numRand=nqRand)
query_meta_tmp <-as.vector(query_meta$Group)
names(query_meta_tmp)<-rownames(query_meta)
grpRand<-rep("rand", nqRand)
names(grpRand)<-paste("rand_", 1:nqRand, sep='')
query_meta_tmp<-append(query_meta_tmp, grpRand)
query_meta_tmp <- assign_cate(query_classification, query_meta_tmp)
query_meta_tmp <- data.frame(query_meta_tmp)

#remove rand id from the data matrix
query_meta_tmp <- query_meta_tmp[-grep('^rand', row.names(query_meta_tmp)), ]
query_meta_tmp <- query_meta_tmp[, 'category', drop = F]
colnames(query_meta_tmp) <- 'singleCellNet_pred'

saveRDS(query_meta_tmp, file = file.path(RESULT_DIR, 'singleCellNet_prediction.rds'))
