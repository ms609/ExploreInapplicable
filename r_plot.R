source('rPlot/functions.R')
source('rPlot/definitions.R')
OVERWRITE <- FALSE

dev.off()
nexusFiles <- list.files('matrices', pattern='.*\\.nex$');# nexusFiles
#nexusName <- 'Loconte1991.nex'

for (nexusName in nexusFiles) {
  nexusRoot <- gsub('.nex', '', nexusName); cat("\nEvaluating", nexusRoot, "...\n")
  
  if (nexusRoot == 'Aguado2009') {
    cat (" ! Problem calculating quartet distaances: Error in 1:(zero.eig[1] - 1): NA/NaN argument")
    next
  }
  if (file.exists(paste0('treeSpaces/', nexusName, '.qt.png', collapse='')) && !OVERWRITE) {
    cat(" - Results already exist.\n")
    next
  }
  rTrees <- lapply(rDirectories, readRTrees, nexusName=nexusName)
  if (is.null(rTrees[[1]])) {
    cat(" ! R trees not found.\n")
    next
  }
  trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName), rTrees)

  nTrees <- vapply(trees, length, integer(1)); names(nTrees) <- allDirectories
  dirTrees <- TreeNumbers(nTrees)
  flatTrees <- unlist(trees, recursive=FALSE)
  treeSource <- rep(c(tntDirectories, rDirectories), nTrees)
  treeTitles <- paste(treeSource, unlist(sapply(nTrees, seq_len)), sep='_')
  treeCol <- paste(rep(treePalette, nTrees))
  treePCh <- rep(plotChars, nTrees)
  
  # Calculate tree scores
  scores <- vapply(allDirectories, function (dirPath) {
    TreeScorer <- if (dirPath %in% tntDirectories) phangorn::fitch else inapplicable::InapplicableFitch
    rawData <- read.nexus.data(paste0(dirPath, '/', nexusName, collapse=''))
    phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
    as.integer(vapply(flatTrees, TreeScorer, double(1), phyData, USE.NAMES=FALSE))
  }, integer(sum(nTrees)))
  rownames(scores) <- treeSource
  minScores <- apply(scores, 2, min)
  extraSteps <- scores - matrix(minScores, nrow(scores), ncol(scores), byrow=TRUE)
  
  # Plot tree scores
  oldPar <- par(mfrow=c(2,2), bg='white')
  for (dirPath in allDirectories) {
    dirScores <- extraSteps[dirTrees[[dirPath]], ]
    dirBreaks <- (min(dirScores) - 1):max(dirScores) + 0.5
    dirCol <- treePalette[dirPath]
    otherDirectories <- which(allDirectories != dirPath)
    
    yMax <- max(apply(dirScores[, otherDirectories], 2, function (x) max(table(x))))
    hist(0, breaks=dirBreaks, border='#ffffffff', ylim=c(0, yMax), axes=FALSE, font.main=1, cex.main=1,
         main=paste0(nexusRoot, " MPTs under ", dirPath), col.main=dirCol, xlab='Extra length') # Set up blank histogram
    axis(1, col=dirCol)
    axis(2, col=dirCol)
    
    for (i in otherDirectories) {
      hist(dirScores[, i], add=TRUE, breaks=dirBreaks,
           border=treePalette[i], col=paste0(treePalette[i], '99', collapse=''))
    }
  }
  dev.copy(svg, file=paste0('islandCounts/', nexusName, '.svg', collapse='')); dev.off()
  dev.copy(png, file=paste0('islandCounts/', nexusName, '.png', collapse='')); dev.off()
  par(oldPar)
  
  # Check that nothing beats the inapplicable trees on the inapplicable measure
  rawData <- read.nexus.data(paste0('inapplicable/', nexusName, collapse=''))
  phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
  scores <- vapply(flatTrees, InapplicableFitch, dataset=phyData, integer(1))
  names(scores) <- treeSource
  table(names(scores)[scores == min(scores)]); min(scores)
  charTypes <- vapply(readLines(paste0('charType/', nexusRoot, '.txt', collapse='')), substr, character(1), 1, 1, USE.NAMES=FALSE)
  charSummary <- paste(names(table(charTypes)), table(charTypes), sep=': ', collapse='; ')
  rfTitleText <- paste(nexusRoot, "R-F space\n",  length(rawData), 'taxa,', attr(phyData, 'nr'), 'chars -', charSummary)
  qtTitleText <- paste(nexusRoot, "Quartet space\n",  length(rawData), 'taxa,', attr(phyData, 'nr'), 'chars -', charSummary)

  # Crude analysis inspired by http://www.kmeverson.org/blog/visualizing-tree-space-in-r
  cat(" - Calculating RF disances...\n")
  rfDistances <- RFDistances(flatTrees);
  cat(" - Calculating quartet disances...\n")
  qtDistances <- QuartetDistances(flatTrees);

  rfSpace <- pcoa(rfDistances)
  PlotTreeSpace(rfSpace, nTrees, legendPos='bottomright', rfTitleText)
  dev.copy(svg, file=paste0('treeSpaces/', nexusName, '.rf.svg', collapse='')); dev.off()
  dev.copy(png, file=paste0('treeSpaces/', nexusName, '.rf.png', collapse='')); dev.off()
  cat(" - Printed RF SVG & PNG.\n")
  qtSpace <- pcoa(qtDistances)
  PlotTreeSpace(qtSpace, nTrees, legendPos='bottomleft', qtTitleText)
  dev.copy(svg, file=paste0('treeSpaces/', nexusName, '.qt.svg', collapse='')); dev.off()
  dev.copy(png, file=paste0('treeSpaces/', nexusName, '.qt.png', collapse='')); dev.off()
  cat("- Printed Quartet SVG & PNG.\n\n")
}

#PlotTreeSpace(pcSpace, nTrees, legendPos=qtLegendPos[[nexusName]], mainTitle=titleText)

# Really we want to something more sophisticated: e.g.
# HILLIS, D. M., HEATH, T. A., JOHN, K. St. and ANDERSON, F. 2005. Analysis and Visualization of Tree Space. Systematic Biology, 54, 471–482.
# or WILGENBUSCH, J. C., HUANG, W. and GALLIVAN, K. A. 2017. Visualizing phylogenetic tree landscapes. BMC Bioinformatics, 18, 85.
