####data path#####
DATA_DIR <- "path/to/data"
RESULT_DIR <- "path/to/result"
TMP_DIR <-  "path/to/tmp"

#####load data pancreas fluidigm data######
pancreas_fluidigm_norm <- readRDS(file = file.path(DATA_DIR, 'pancreas_fluidigmc1_normed.rds'))
pancreas_fluidigm_count <- readRDS(file = file.path(DATA_DIR, 'pancreas_fluidigmc1_count.rds'))

####calculate the gene counts#######
pancreas_flu_gene_count <- rowSums(pancreas_fluidigm_count)

####rule out genes with 0 TPM in all cells####
pancreas_flu_gene_count <- pancreas_flu_gene_count[pancreas_flu_gene_count!=0]

##start with 21563 genes######
####original gene distribution####
pdf(file = file.path(RESULT_DIR, 'Original gene count distribution in pancreas fluidigm.pdf'))
hist(log10(pancreas_flu_gene_count), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"), ylim = c(0, 2000), xlim = c(0, 10))
abline(v=log10(1), col="blue", lwd=2, lty=2)
dev.off()

####process the data to extract probability distribution for each gene####
pancreas_flu_gene_count <- as.factor(pancreas_flu_gene_count)
pancreas_flu_gene_count <- as.data.frame(pancreas_flu_gene_count)
pancreas_flu_gene_count_freq <- table(pancreas_flu_gene_count)
pancreas_flu_gene_count_prop <- prop.table(pancreas_flu_gene_count_freq)
pancreas_flu_gene_count_prop <- as.data.frame(pancreas_flu_gene_count_prop)

write.csv(pancreas_flu_gene_count, file = file.path(TMP_DIR, 'pancreas_flu_gene_count_tmp.csv'), row.names = T)
pancreas_flu_gene_count <- read.csv(file = file.path(TMP_DIR, 'pancreas_flu_gene_count_tmp.csv'))
colnames(pancreas_flu_gene_count) <- c('gene_name', 'counts')
pancreas_flu_gene_count$counts <- as.factor(pancreas_flu_gene_count$counts)

result_prop_table <- merge(pancreas_flu_gene_count_prop, pancreas_flu_gene_count, all.x = T, all.y = T, by.y = 'counts', by.x = 'pancreas_flu_gene_count')

####sample####
gene_list <- result_prop_table$gene_name
prob_list <- result_prop_table$Freq

#####function to sample the dataset#####
down_sample_gene <- function(data_count, data_norm, gene_list, prob_list, size){
  for (i in c(1:5)){
    set.seed(i)
    gene_sample <- sample(gene_list, size = size, replace = FALSE, prob = prob_list)
    data_count_sample <- data_count[gene_sample, ]
    data_norm_sample <- data_norm[gene_sample, ]
    saveRDS(data_count_sample, file = file.path(RESULT_DIR, paste('pancreas_fluidigm_downsample_gene_', size, '_sample_', i, '_count.rds', sep = '')))
    saveRDS(data_norm_sample, file = file.path(RESULT_DIR, paste('pancreas_fluidigm_downsample_gene_', size, '_sample_', i, '_norm.rds', sep = '')))
  }
  
}

####code to construct all the sampling#####
for (i in c(5000, 10000, 15000)){
  down_sample_gene(pancreas_fluidigm_count, pancreas_fluidigm_norm, gene_list, prob_list, i)
}


#####plot log count to check the distribution####
pdf(file = file.path(RESULT_DIR, 'Downsampled gene count distribution in pancreas fluidigm.pdf'))
for (size in c(5000, 10000, 15000)){
  for (i in c(1:5)){
    data <- readRDS(file = file.path(RESULT_DIR, paste('pancreas_fluidigm_downsample_gene_', size, '_sample_', i, '_count.rds', sep = '')))
    data <- rowSums(data)
    hist(log10(data), breaks=100, main="", col="grey80",
         xlab=expression(Log[10]~"average count"), ylim = c(0, 2000), xlim = c(0, 10))
    abline(v=log10(1), col="blue", lwd=2, lty=2)
  }
}
dev.off()
