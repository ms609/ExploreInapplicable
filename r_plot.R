source('r_plot_functions.R')
treeLegendPos = list(
  'Asher2005.nex' = 'topright',
  'Aria2015' = 'topright',
  'DeAssis2011.nex' = 'bottomright',
  'Eklund2004.nex' = 'bottomright',
  'Griswold1999.nex' = 'bottomright',
  'OLeary1999.nex' = 'bottomright',
  'Rousset2004' = 'bottomleft',
  'Sano2011.nex' = 'bottomleft',
  'Sansom2010.nex' = 'bottomright',
  'Wortley2006.nex' = 'bottomleft',
  'Geisler2001.nex' = 'bottomleft'
)

tntDirectories <- c('ambiguous', 'ambigAbsent', 'extraState')
rDirectories <- c('inapplicable')

dev.off()
nexusFiles <- list.files('matrices', pattern='.*\\.nex$');# nexusFiles
nexusName <- 'Loconte1991.nex'; nexusRoot <- gsub('.nex', '', nexusName); cat("Evaluating ", nexusRoot, "\n")

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
rfTitleText <- paste(nexusRoot, "R-F space\n",  length(rawData), 'taxa,', attr(phyData, 'nr'), 'chars -', charSummary)
qtTitleText <- paste(nexusRoot, "Quartet space\n",  length(rawData), 'taxa,', attr(phyData, 'nr'), 'chars -', charSummary)

# Crude analysis inspired by http://www.kmeverson.org/blog/visualizing-tree-space-in-r
rfDistances <- RFDistances(flatTrees);
qtDistances <- QuartetDistances(flatTrees);

rfSpace <- pcoa(rfDistances)
PlotTreeSpace(rfSpace, nTrees, legendPos='bottomright', rfTitleText)
#dev.copy(svg, file=paste0('treeSpaces/', nexusName, '.rf.svg', collapse='')); dev.off()
qtSpace <- pcoa(qtDistances)
PlotTreeSpace(qtSpace, nTrees, legendPos='bottomright', qtTitleText)
#dev.copy(svg, file=paste0('treeSpaces/', nexusName, '.qt.svg', collapse='')); dev.off()


#PlotTreeSpace(pcSpace, nTrees, legendPos=treeLegendPos[[nexusName]], mainTitle=titleText)

# Really we want to something more sophisticated: e.g.
# HILLIS, D. M., HEATH, T. A., JOHN, K. St. and ANDERSON, F. 2005. Analysis and Visualization of Tree Space. Systematic Biology, 54, 471–482.
# or WILGENBUSCH, J. C., HUANG, W. and GALLIVAN, K. A. 2017. Visualizing phylogenetic tree landscapes. BMC Bioinformatics, 18, 85.
