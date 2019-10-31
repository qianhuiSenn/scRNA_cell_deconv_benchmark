#######################################
##### CP and RPC sample usage ######
#######################################
#Config########
library(EpiDISH)
library(matrixStats)
library(dplyr)
library(Seurat)

set.seed(3456)

DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"

##config##
rpcmax <- 100

##load query and ref matrix
query <- readRDS(file = file.path(DATA_DIR, 'query_norm.rds'))
ref <- readRDS(file = file.path(DATA_DIR, 'ref_psudo_bulk.rds'))

####Feature Selection top 2000######
reference_count <- readRDS(file = file.path(DATA_DIR, 'ref_count.rds'))
reference <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "reference")
reference <- NormalizeData(object = reference)
reference <- FindVariableFeatures(object = reference, selection.method = 'vst', nfeatures = 2000)
var_gene <- head(VariableFeatures(reference), 2000)

###further process ref psudo_bulk and query matrix using var gene ###
sharegene <- intersect(var_gene, row.names(ref))
ref <- ref[sharegene, ]
ref <- ref[is.finite(rowSums(ref)),]
refgenes <- row.names(ref)
query <- query[is.finite(rowSums(query)),]
scalegenes <- row.names(query)
sharedgenes <- intersect(refgenes, scalegenes)
ref <- ref[sharedgenes, ]
ref <- t(ref)
ref <- t(ref)
query <- query[sharedgenes,]


#######Cell deconvolution########
##cp###
cpresult <- epidish(beta.m = query, ref.m = ref, method = 'CP', 
                    constraint = c('inequality'))
cpcov <- cpresult$estF
cpcov <- t(cpcov)
cpcov <- t(cpcov)
cpcov <- as.data.frame(cpcov)
#cpcov$Undefined <- 0
cellsums <- rowSums(cpcov)
#cpcov$Undefined <- 1 - cellsums
cpcov <- t(cpcov)
cpcov <- t(cpcov)
prevalues <- rowMaxs(cpcov)
boolmatrix <- cpcov == prevalues
prelabels <- c()
candidates <- colnames(cpcov)
getprelabel <- function(line){
  prelabel <- candidates[line]
  prelabels <- c(prelabels, prelabel)
}
labelresult <- apply(boolmatrix, 1, getprelabel)
CP_pred <- as.data.frame(labelresult)
colnames(CP_pred) <- 'CP_pred'
saveRDS(CP_pred,  file = file.path(RESULT_DIR,'CP_prediction.rds'))

##RPC###
rpcresult <- epidish(beta.m = query, ref.m = ref, method = 'RPC', 
                     maxit = rpcmax)
rpccov <- rpcresult$estF
rpccov <- t(rpccov)
rpccov <- t(rpccov)
rpccov <- as.data.frame(rpccov)
rpccov$Undefined <- 0
cellsums <- rowSums(rpccov)
rpccov$Undefined <- 1 - cellsums
rpccov <- t(rpccov)
rpccov <- t(rpccov)
prevalues <- rowMaxs(rpccov)
boolmatrix <- rpccov == prevalues
prelabels <- c()
candidates <- colnames(rpccov)
getprelabel <- function(line){
  prelabel <- candidates[line]
  prelabels <- c(prelabels, prelabel)
}
labelresult <- apply(boolmatrix, 1, getprelabel)
RPC_pred <- as.data.frame(labelresult)
colnames(RPC_pred) <- 'RPC_pred'
saveRDS(RPC_pred,  file = file.path(RESULT_DIR,'RPC_prediction.rds'))