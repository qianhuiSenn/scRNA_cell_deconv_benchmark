# Required Library
library(splatter)
RESULT_DIR <- 'path/to/result'

# 1. Simulation data 1 in manuscript
params.groups <- newSplatParams(batchCells = 2000, nGenes = 4000)

sim1 <- splatSimulateGroups(params.groups, group.prob = c(0.1, 0.3, 0.3, 0.1, 0.2), de.prob = c(0.3, 0.2, 0.2, 0.3, 0.4),
                            de.facScale = c(0.2, 0.5, 0.2, 0.5, 0.4),
                            dropout.shape=-0.5, dropout.mid=1,
                            dropout.type="experiment", verbose = F)

raw <- assays(sim1)[["counts"]]
truth <- assays(sim1)[["TrueCounts"]]

# Write metadata
write.csv(colData(sim1),"cell_info_sim1.csv")
write.csv(rowData(sim1),"gene_info_sim1.csv")

# Write data
write.csv(raw,"raw_assay1.csv")
write.csv(truth,"truth_assay1.csv")

# Save global R object
saveRDS(sim1,"splatterSim1.RDS")

# 2. simulation data 2 in manuscript (different de scales)
seed_list <- c(202001, 202002, 202003, 202004, 202005)
scale_list <- list(c(0.1, 0.3, 0.1, 0.3, 0.2), c(0.3, 0.5, 0.3, 0.5, 0.4), c(0.5, 0.7, 0.5, 0.7, 0.6), c(0.7, 0.9, 0.7, 0.9, 0.8))
name_list <- list('Low', 'Low_Moderate', ' Moderate', 'High')

sim_run <- function(scale, name, seed_opt){
  params.groups <- newSplatParams(batchCells = 2000, nGenes = 10000, seed = seed_opt)
  sim <- splatSimulateGroups(params.groups, group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2), de.prob = c(0.3, 0.3, 0.3, 0.3, 0.3),
                             de.facScale = scale,
                             dropout.shape=-0.5, dropout.mid=1,
                             dropout.type="experiment", verbose = F)
  #extract assay from each simulations
  truth <- assays(sim)[["TrueCounts"]]
  raw <- assays(sim)[["counts"]]
  cell_info <- colData(sim)
  cell_info <- data.frame(cell_info[, c(1, 3)], row.names = 1)
  # Write metadata
  write.csv(cell_info,paste("cell_info_sim2", name, seed_opt, '.csv', sep = ''))
  # Write data
  write.csv(raw, paste("raw_assay2", name, seed_opt, '.csv', sep = ''))
  write.csv(truth,paste("true_assay2", name, seed_opt, '.csv', sep = ''))
  # save simulation data
  saveRDS(sim, paste("splatterSim2", name, seed_opt, ".RDS", sep = ''))
}

for (i in c(1:4)){
  for (j in seed_list){
    sim_run(scale_list[[i]], name_list[[i]], j)
  }
}


# 3. simulation data 3 with increased classification labels in manuscript

ann_level <- list(10, 20, 30, 40, 50)
name_list <- list('10', '20', ' 30', '40', '50')

sim_run <- function(scale, name){
  set.seed(1234)
  params.groups <- newSplatParams(batchCells = 10000, nGenes = 20000)
  sim <- splatSimulateGroups(params.groups, group.prob = rep(1/scale, scale),
                             dropout.shape=-0.5, dropout.mid=1,
                             dropout.type="experiment", verbose = F)
  #extract assay from each simulations
  truth <- assays(sim)[["TrueCounts"]]
  raw <- assays(sim)[["counts"]]
  cell_info <- colData(sim)
  cell_info <- data.frame(cell_info[, c(1, 3)], row.names = 1)
  # Write metadata
  write.csv(cell_info, file = file.path(RESULT_DIR, paste("cell_info_sim3", name, '.csv', sep = '')))
  # Write data
  saveRDS(truth, file = file.path(RESULT_DIR, paste('true_simulation3_count',name,'group.rds', sep = '_')))
  saveRDS(raw, file = file.path(RESULT_DIR, paste('raw_simulation3_count',name,'group.rds', sep = '_')))
}

for (i in c(1:5)){
  sim_run(ann_level[[i]], name_list[[i]])
}

# 4. simulation data 4 with descending cell proportion for each cell group

seed_list <- c(202001, 202002, 202003, 202004, 202005, 202006, 202007, 202008, 202009, 202010)
for (seed_opt in seed_list){
  params.groups <- newSplatParams(batchCells = 2000, nGenes = 10000, seed = seed_opt)
  sim2 <- splatSimulateGroups(params.groups, group.prob = c(1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 10/1024, 4/1024, 2/1024),
                           dropout.shape=-0.5, dropout.mid=1,
                           dropout.type="experiment", verbose = F)
  
  raw <- assays(sim2)[["counts"]]
  truth <- assays(sim2)[["TrueCounts"]]
  cell_info <- colData(sim2)
  cell_info <- data.frame(cell_info[, c(1, 3)], row.names = 1)
  
  # Write metadata
  write.csv(cell_info, file = file.path(RESULT_DIR, paste("cell_info_rare_group_sim_", seed_opt, '.csv', sep = '')))
  # Write data
  saveRDS(truth, file = file.path(RESULT_DIR, paste('true_rare_group_sim_count',seed_opt,'.rds', sep = '_')))
  saveRDS(raw, file = file.path(RESULT_DIR, paste('raw_rare_group_sim_count',seed_opt,'.rds', sep = '_')))
}


# 5. Simulation data 5 in time and memory comparison

batchCells <- c(5000, 10000, 15000, 20000, 25000, 50000)
for (i in batchCells){
  params.groups <- newSplatParams(batchCells = i, nGenes = 20000)
  sim <- splatSimulateGroups(params.groups, group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2),
                             dropout.shape=-0.5, dropout.mid=1,
                             dropout.type="experiment", verbose = F)
  #extract assay from each simulations
  truth <- assays(sim)[["TrueCounts"]]
  raw <- assays(sim)[["counts"]]
  cell_info <- colData(sim)
  cell_info <- data.frame(cell_info[, c(1, 3)], row.names = 1)
  saveRDS(true, file = file.path(RESULT_DIR, paste('true_simulation5_count',i,'Cell.rds', sep = '_')))
  saveRDS(raw, file = file.path(RESULT_DIR, paste('raw_simulation5_count',i,'Cell.rds', sep = '_')))
  write.csv(cell_info,paste("cell_info_sim5", name, '.csv', sep = ''))
  saveRDS(sim,paste("splatterSim5_", i, ".RDS", sep = ''))
}


