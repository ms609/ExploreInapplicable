source('rPlot/functions.R')
source('rPlot/definitions.R')

props <- read.csv('matrixProperties.csv', header=TRUE)
countDims <- vapply(rownames(props), MatrixProperties, integer(6))

validReads <- props$nTax == countDims['nTax',]
props[validReads, c('nTokens', 'chars_na', 'tokens_na', 'tokens_amb')] <- 
  t(countDims[c('nTokens', 'nInapp', 'inappTokens', 'ambigTokens'), validReads])
props$prop_na <- props$tokens_na / props$nTokens
props$prop_chars_na <- props$chars_na / props$nChars
rownames(props) <- props$matrix

### Analyse precision:
### Numbers of MPTs
csvs <- list.files('islandCounts', pattern='.*\\.csv$');

mptCount <- vapply(csvs, function (csvFile) {
  content <- read.csv(paste0('islandCounts/', csvFile))
  table(content$X)
}, double(4))[-1, ]
colnames(mptCount) <- gsub('\\.csv', '', colnames(mptCount))
mpts <- t(mptCount)
props[, c('ambigMPTs', 'exstMPTs', 'inappMPTs')] <- NA
props[colnames(mptCount), c('ambigMPTs', 'exstMPTs', 'inappMPTs')] <- t(mptCount)




nexusFiles <- list.files('matrices', pattern='.*\\.nex$');# nexusFiles

for (nexusName in nexusFiles) {
  nexusRoot <- gsub('.nex', '', nexusName)
  if (!any(is.na(props[nexusRoot, c('consNodes_ambig', 'consNodes_all')]))) next
  cat("\nEvaluating", nexusRoot, "...\n")
 
  rTrees <- lapply(rDirectories, readRTrees, nexusName=nexusName)
  if (is.null(rTrees[[1]])) {
    cat(" ! R trees not found.\n")
    next
  }
  trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName), rTrees)
  cat(" - Trees loaded. Calculating consensus trees: ")
  consensi <- lapply(trees, consensus)
  nNodes <- vapply(consensi, function (tr) tr$Nnode, double(1))
  cat("done.\n")
  props[nexusRoot, c('consNodes_ambig', 'consNodes_exst', 'consNodes_inapp')] <- nNodes[2:4]
  
  props[nexusRoot, 'consNodes_ambig_exst'] <- consensus(consensi[[2]], consensi[[3]])$Nnode
  props[nexusRoot, 'consNodes_ambig_inapp']  <- consensus(consensi[[2]], consensi[[4]])$Nnode
  props[nexusRoot, 'consNodes_exst_inapp'] <- consensus(consensi[[4]], consensi[[3]])$Nnode
  props[nexusRoot, 'consNodes_all'] <- consensus(consensi[[4]], consensi[[3]], consensi[[2]])$Nnode
}

write.csv(props, 'matrixProperties.csv')

plot(I(consNodes_exst_inapp / consNodes_inapp) ~ prop_na, data=props)
plot(I(consNodes_inapp / nTax) ~ prop_na, data=props)
plot(I(consNodes_exst_inapp / consNodes_inapp) ~ prop_chars_na, data=props)
plot(I(consNodes_exst_inapp / consNodes_inapp) ~ nTokens, data=props)

plot(I(inappMPTs) ~ prop_na, data=props)
plot(I(ambigMPTs / ) ~ prop_na, data=props)


### Analyse accuracy:
