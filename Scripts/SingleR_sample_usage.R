################################
##### SingleR sample usage ######
################################

library(Matrix)
library(SingleR)
library(ggplot2)
library(knitr)
library(genefilter)
set.seed(3456)
DATA_DIR <- "/path/to/data"
RESULT_DIR <- "/path/to/result"
TMP_DIR <-  "/path/to/tmp"

# load psudo_bulk_reference
exp <- readRDS(file = file.path(DATA_DIR,'ref_psudo_bulk.rds'))
exp <- as.matrix(exp)
name = 'example_ref'
types <- as.character(colnames(exp))

example_ref = list(name = name, data = exp, types = types)
example_ref$de.genes = CreateVariableGeneSet(exp, types, 2000)

save(example_ref, file = file.path(RESULT_DIR, 'example_ref.RData'))

# load normalized qeury data
query <- readRDS(file = file.path(DATA_DIR, 'query_norm.rds'))

example_singler <- CreateSinglerObject(counts = query, cluster = NULL, annot = NULL, project.name = "example singler",
                                   min.genes = 200, technology = "",
                                   species = "", citation = "", ref.list = list(example_ref),
                                   normalize.gene.length = F, variable.genes = "de", fine.tune = T,
                                   reduce.file.size = T, do.signatures = F, do.main.types = F,
                                   temp.dir = TMP_DIR, numCores = SingleR.numCores)
saveRDS(example_singler$singler[[1]]$SingleR.single$labels, file = file.path(RESULT_DIR, 'singler_prediction.rds'))