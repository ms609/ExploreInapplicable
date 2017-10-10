source('rPlot/functions.R')
source('rPlot/definitions.R')
library('VennDiagram')

props <- read.csv('matrixProperties.csv', header=TRUE, row.names=1)
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

props[, c('MDSarea_ambig', 'MDSarea_exst', 'MDSarea_inapp')]  <- NA
for (fileRoot in names(validReads)) {
  if (all(file.exists(paste0('treeSpaces/', fileRoot, '.qt.csv'), paste0('islandCounts/', fileRoot, '.csv')))) {
    treeDetails <- read.csv(paste0('islandCounts/', fileRoot, '.csv'))
    nTrees <- table(treeDetails[, 1])
    qtDists <- data.matrix(read.csv(paste0('treeSpaces/', fileRoot, '.qt.csv'), row.names=1))
    qtDists[qtDists == 0] <- 1e-9
    diag(qtDists) <- 0
    props[fileRoot, c('MDSarea_ambig', 'MDSarea_exst', 'MDSarea_inapp')] <- PlotKruskalTreeSpace3(qtDists, nTrees)
  }
}

vennTreeNames <- c('utrees_ambig', 'utrees_exst', 'utrees_inapp', 'trees_ambig_exst', 
          'trees_ambig_inapp', 'trees_exst_inapp', 'trees_all')
props[, vennTreeNames] <- NA
for (fileRoot in names(validReads)) {
  nexusName <- paste0(fileRoot, '.nex')
  trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName),
             lapply(rDirectories, readRTrees, nexusName=nexusName))
  nTrees <- vapply(trees, length, integer(1))
  vennTrees <- integer(7)
  vennTrees[1:3] <- nTrees[2:4]
  vennTrees[4:7] <- vapply(list(trees[c(2, 3)], trees[c(2, 4)], trees[c(3, 4)], trees[2:4]),
    function (x) length(unique(unlist(x, recursive=FALSE))), integer(1))
  vennTrees[7] <- sum(vennTrees[1:3]) - vennTrees[7]
  vennTrees[4:6] <- vapply(3:1, function (i) sum(vennTrees[1:3][-i]), integer(1)) - vennTrees[4:6]
  props[fileRoot, vennTreeNames] <- vennTrees
  
}



write.csv(props, 'matrixProperties.csv')

vennProps <- props[!is.na(props$consNodes_ambig), c('consNodes_ambig', 'consNodes_exst', 'consNodes_inapp', 
                 'consNodes_ambig_exst', 'consNodes_exst_inapp', 'consNodes_ambig_inapp', 'consNodes_all')]
vennCons <- colSums(vennProps)
names(vennCons) <- c('A', 'B', 'C', 'A&B', 'A&C', 'B&C', 'A&B&C')
vennCons[4:6] <- vennCons[4:6] - vennCons['A&B&C']
vennCons[1:3] <- c(
'A' =     vennCons['A'] - (vennCons['A&B'] + vennCons['A&C'] + vennCons['A&B&C']),
'B' =     vennCons['B'] - (vennCons['A&B'] + vennCons['B&C'] + vennCons['A&B&C']),
'C' =     vennCons['C'] - (vennCons['A&C'] + vennCons['B&C'] + vennCons['A&B&C']))
vennPlot <- venneuler(vennCons)
vennPlot
plot(vennPlot, col=treePalette[2:4], col.fn=function(x) x, border=treePalette[2:4], edges=1024, col.txt=NA, 
     main='Node presence by method: all datasets')
dev.copy(svg, file='vennTrees/_All_datasets.svg'); dev.off()
dev.copy(png, file='vennTrees/_All_datasets.png', width=800, height=800); dev.off()
     
if (!file.exists(paste0('vennTrees/', nexusRoot, '.png')) || OVERWRITE) {
  vennTrees <- integer(7)
  names(vennTrees) <- c('A', 'B', 'C', 'A&B', 'A&C', 'B&C', 'A&B&C')
  vennTrees[1:3] <- nTrees[2:4]
  vennTrees[4:7] <- vapply(list(trees[c(2, 3)], trees[c(2, 4)], trees[c(3, 4)], trees[2:4]),
    function (x) length(unique(unlist(x, recursive=FALSE))), integer(1))
  vennTrees[7] <- sum(vennTrees[1:3]) - vennTrees[7]
  vennTrees[4:6] <- vapply(3:1, function (i) sum(vennTrees[1:3][-i]), integer(1)) - vennTrees[4:6]
  vennPlot <- venneuler::venneuler(vennTrees)
  plot(vennPlot, col=treePalette[2:4], col.fn=function(x) x, border=treePalette[2:4],
      col.txt=NA, edges=1024, main=paste0('Shortest trees: ', nexusRoot))
  dev.copy(svg, file=paste0('vennTrees/', nexusRoot, '.svg', collapse='')); dev.off()
  dev.copy(png, file=paste0('vennTrees/', nexusRoot, '.png', collapse=''), width=800, height=800); dev.off()
  cat(" - Plotted Venn diagram of trees.\n")
}


