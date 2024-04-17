# main_analysis.R

# Load the data
library(diffCircadian)
tt <- read.csv('Data/data_impute_zeit.csv')$Cage
data <- read.csv('Data/data_impute_zeit.csv')[, -1]
pheno <- read.csv('Data/pheno.csv')

# Run the analysis
for (phenotype in c('A', 'B', 'C')) {
  phenotype_samples <- pheno$Cage[which(pheno$Dose == phenotype)]
  phenotype_data <- data[, which(colnames(data) %in% phenotype_samples)]
  thisResult <- do.call(rbind, lapply(1:length(phenotype_samples), function(i, phenotype) {
    # data_total <- data[which(pheno$
    data_for_one_cage <- as.numeric(phenotype_data[,i])
    aLR <- diffCircadian::LR_rhythmicity(tt, data_for_one_cage)
    thisResult_row <- list(aLR$amp, aLR$phase, aLR$peakTime, aLR$offset, aLR$pvalue, aLR$R2)
    return(thisResult_row)
  }))
  
  colnames(thisResult) <- c('amp', 'phase', 'peakTime', 'offset', 'pvalue', 'R2')
  rownames(thisResult) <- phenotype_samples 
  dir.create('Results', showWarnings = FALSE)
  write.csv(file=paste('Results/phenotype_', phenotype, '_results.csv', sep=''), thisResult)
}
