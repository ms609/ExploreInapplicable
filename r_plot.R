source('r_plot_functions.R')
treeLegendPos = list(
  'Asher2005.nex' = 'topright'
)

tntDirectories <- c('ambiguous', 'ambigAbsent', 'extraState')
rDirectories <- c('inapplicable')

nexusFiles <- list.files('matrices', pattern='.*\\.nex$'); nexusFiles
nexusName <- nexusFiles[2]; cat("Evaluating ", nexusName, "\n")

trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName), 
              lapply(rDirectories, readRTrees, nexusName=nexusName))

nTrees <- vapply(trees, length, integer(1))
flatTrees <- unlist(trees, recursive=FALSE)
treeTitles <- paste(rep(c(tntDirectories, rDirectories), nTrees), unlist(sapply(nTrees, seq_len)))
treeCol <- paste(rep(treePalette, nTrees))
treePCh <- rep(plotChars, nTrees)

# Check that nothing beats the inapplicable trees on the inapplicable measure
rawData <- read.nexus.data(paste0('inapplicable/', nexusName, collapse=''))
phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
scores <- vapply(flatTrees, InapplicableFitch, dataset=phyData, integer(1))
names(scores) <- treeTitles
names(scores)[scores == min(scores)]

# Crude analysis inspired by http://www.kmeverson.org/blog/visualizing-tree-space-in-r
allDistances <- RFDistances(flatTrees);
PlotTreeSpace(pcoa(allDistances), nTrees, legendPos=treeLegendPos[[nexusName]])

#dev.copy2pdf(file=paste0('treeSpaces/', nexusName, '.pdf', collapse=''))

# Really we want to something more sophisticated: e.g.
# HILLIS, D. M., HEATH, T. A., JOHN, K. St. and ANDERSON, F. 2005. Analysis and Visualization of Tree Space. Systematic Biology, 54, 471–482.
# or WILGENBUSCH, J. C., HUANG, W. and GALLIVAN, K. A. 2017. Visualizing phylogenetic tree landscapes. BMC Bioinformatics, 18, 85.
