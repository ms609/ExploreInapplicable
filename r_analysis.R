setwd("C:/Research/ExploreInapplicable")

source('rPlot/functions.R')
source('rPlot/definitions.R')

props <- read.csv('matrixProperties.csv', header=TRUE, row.names=1)
countDims <- vapply(rownames(props), MatrixProperties, integer(6))

validReads <- props$nTax == countDims['nTax',]

#### Plot QUARTET TREESPACES ############
nValid <- sum(validReads)
#par(mfrow=c(6, 5), mar=rep(0.2, 4), bg='white')
par(bg='white')
quartPlots <- sort(names(validReads)[validReads])
for (fileRoot in quartPlots) {
  if (fileRoot %in% slowFiles) {cat("\n x ", fileRoot, "\n"); next}
  cat("\n - ", fileRoot, "\n")
  if (file.exists(paste0('quartetSpaces/', fileRoot, '.png'))) {
    cat("   > [File exists].\n")
    next
  }
  trees <- GetTrees(fileRoot)
  nTrees <- vapply(trees, length, integer(1))
  qtDistances3   <- GetQuartetDistances(fileRoot, trees[-1], forPlot=TRUE)
  PlotKruskalTreeSpace3(qtDistances3, nTrees[-1], legendPos=QuartetLegendPos(fileRoot), fileRoot, fill=TRUE)
  dev.copy(svg, file=paste0('quartetSpaces/', fileRoot, '.svg', collapse='')); dev.off()
  dev.copy(png, file=paste0('quartetSpaces/', fileRoot, '.png', collapse=''), width=1024, height=1024); dev.off()  
}






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





############ VENN DIAGRAMS  ~~~##################

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
for (fileRoot in names(validReads)[validReads]) {
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
for (fileRoot in names(validReads)[validReads]) {
  #if (!is.na(props[fileRoot, vennTreeNames[6]])) next
  cat("\n - ", fileRoot, "\n")
  vennTrees <- GetVennTrees(fileRoot)
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
dev.copy(svg, file='vennNodes/_All_datasets.svg'); dev.off()
dev.copy(png, file='vennNodes/_All_datasets.png', width=800, height=800); dev.off()
   
vennTrees <- colSums(props[validReads, vennTreeNames])
names(vennTrees) <- c('A', 'B', 'C', 'A&B', 'A&C', 'B&C', 'A&B&C')
vennPlot <- venneuler::venneuler(vennTrees)
plot(vennPlot, col=treePalette[2:4], col.fn=function(x) x, border=treePalette[2:4],
    col.txt=NA, edges=1024, main='Shortest trees: all datasets') # a-c clockwise from bottom
text(vennPlot$centers[, 1] + c(0, -0.05, 0.05), vennPlot$centers[, 2] + c(-0.08, 0.05, 0.05), paste0(englishName[2:4], ': ', vennTrees[1:3]))
text(mean(vennPlot$centers[, 1]) - 0.001, mean(vennPlot$centers[, 2]) + 0.03, vennTrees[7])
offX <- c(-0.005, 0.02, -0.02)
offY <- c(0.056, 0.00, 0.00)
for (i in 1:3) text(mean(vennPlot$centers[-i, 1]) + offX[i], mean(vennPlot$centers[-i, 2]) + offY[i], vennTrees[7-i])

dev.copy(svg, file='vennTrees/_All_datasets.svg'); dev.off()
dev.copy(png, file='vennTrees/_All_datasets.png', width=400, height=400); dev.off()



