if (!dir.exists('rPlot')) setwd('../') # Might need to move up one directory
source('rPlot/functions.R')
source('rPlot/definitions.R')
OVERWRITE <- FALSE

# nexusName <- 'Asher2005.nex'

props <- read.csv('matrixProperties.csv', header=TRUE, row.names=1)

par(mfrow=c(5, 6), bg='white', mar=rep(0.78, 4))
for (nexusName in sort(nexusFiles)) {
#for (nexusName in rep(nexusFiles[1], 30)) {
  nexusRoot <- gsub('.nex', '', nexusName)
  cat("\nEvaluating", nexusRoot, "...\n")

  trees <- GetTrees(nexusRoot)
  cat(" - Trees read OK.\n")
  nTrees <- vapply(trees, length, integer(1)); names(nTrees) <- allDirectories
  dirTrees <- TreeNumbers(nTrees)
  flatTrees <- unlist(trees, recursive=FALSE)
  
  treeSource <- rep(c(tntDirectories, rDirectories), nTrees)
  treeTitles <- paste(treeSource, unlist(sapply(nTrees, seq_len)), sep='_')
  treeCol <- paste(rep(treePalette, nTrees))
  treePCh <- rep(plotChars, nTrees)
 
  charTypes <- vapply(readLines(paste0('charType/', nexusRoot, '.txt', collapse='')), substr, character(1), 1, 1, USE.NAMES=FALSE)
  # This will misestimate in the matrix where chars without inapps were coded N or T instead of being left as X.
  charSummary <- paste(names(table(charTypes)), table(charTypes), sep=': ', collapse='; ')
  
  qtTitleText <- paste(nexusRoot, "Quartet space\n",  props[nexusRoot, 'nTax'], 'taxa,', props[nexusRoot, 'nChars'], 'chars -', charSummary)
  qtDistances <- GetQuartetDistances(nexusRoot, trees, forPlot=TRUE)
  cat(" - Got Quartet distances.\n")
  
  ambigTrees <- seq_len(nTrees[1])
  TreeSpacePanel(qtDistances[-ambigTrees, -ambigTrees], nTrees[-1], studyName[[nexusRoot]], 1.0)
  cat(" - Plotted treespaces.\n\n")
}
dev.copy(svg, file=paste0('publication_figures/Treespace_raw.svg', collapse='')); dev.off()
dev.copy(pdf, file=paste0('publication_figures/Treespace_raw.pdf', collapse=''), width=8, height=6); dev.off()
