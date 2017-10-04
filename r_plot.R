source('r_plot_functions.R')
treeLegendPos = list(
  'Asher2005.nex' = 'topright',
  'Wortley2006.nex' = 'bottomleft',
  'Geisler2001.nex' = 'bottomleft'
)

tntDirectories <- c('ambiguous', 'ambigAbsent', 'extraState')
rDirectories <- c('inapplicable')

nexusFiles <- list.files('matrices', pattern='.*\\.nex$'); nexusFiles
nexusName <- 'Wortley2006.nex'; nexusRoot <- gsub('.nex', '', nexusName); cat("Evaluating ", nexusRoot, "\n")

trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName), 
              lapply(rDirectories, readRTrees, nexusName=nexusName))

nTrees <- vapply(trees, length, integer(1))
flatTrees <- unlist(trees, recursive=FALSE)
treeSource <- rep(c(tntDirectories, rDirectories), nTrees)
treeTitles <- paste(treeSource, unlist(sapply(nTrees, seq_len)))
treeCol <- paste(rep(treePalette, nTrees))
treePCh <- rep(plotChars, nTrees)


# Check that nothing beats the inapplicable trees on the inapplicable measure
rawData <- read.nexus.data(paste0('inapplicable/', nexusName, collapse=''))
phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
scores <- vapply(flatTrees, InapplicableFitch, dataset=phyData, integer(1))
names(scores) <- treeSource
table(names(scores)[scores == min(scores)]); min(scores)
charTypes <- vapply(readLines(paste0('charType/', nexusRoot, '.txt', collapse='')), substr, character(1), 1, 1, USE.NAMES=FALSE)
charSummary <- paste(names(table(charTypes)), table(charTypes), sep=': ', collapse='; ')
titleText <- paste(nexusRoot, '-',  length(rawData), 'taxa,', attr(phyData, 'nr'), 'chars -', charSummary)

# Crude analysis inspired by http://www.kmeverson.org/blog/visualizing-tree-space-in-r
allDistances <- RFDistances(flatTrees);
pcSpace <- pcoa(allDistances)
PlotTreeSpace(pcSpace, nTrees, legendPos='bottomleft', titleText)

PlotTreeSpace(pcSpace, nTrees, legendPos=treeLegendPos[[nexusName]], mainTitle=titleText)
#dev.copy2pdf(file=paste0('treeSpaces/', nexusName, '.pdf', collapse=''))

# Really we want to something more sophisticated: e.g.
# HILLIS, D. M., HEATH, T. A., JOHN, K. St. and ANDERSON, F. 2005. Analysis and Visualization of Tree Space. Systematic Biology, 54, 471–482.
# or WILGENBUSCH, J. C., HUANG, W. and GALLIVAN, K. A. 2017. Visualizing phylogenetic tree landscapes. BMC Bioinformatics, 18, 85.
