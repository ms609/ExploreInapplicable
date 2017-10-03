install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
library(phangorn)
devtools::install_github('ms609/inapplicable')
require(inapplicable)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
treePalette <- cbPalette[c(6, 3, 2, 8)]
plotChars <- c(2, 1, 3, 4)

readTntTrees <- function (directory, nexusName) {
  unique(read.nexus(paste0(directory, '/', nexusName, '.nextrees', collapse='')))
}
readRTrees <- function (directory, nexusName) {
  allResults <- list.files(directory, paste0(nexusName, '.*\\-[[:digit:]]+.tre', collapse=''))
  resultScores <- vapply(allResults, function (string) {
    hits <- regexpr(pattern='\\-[[:digit:]]+', string)
    return(as.integer(substr(string, hits[1] + 1, hits[1] + attr(hits, 'match.length') - 1)))
  }, integer(1))
  unique(read.tree(paste0(directory, '/', allResults[which.min(resultScores)], collapse='')))
}
RFDistances <- function(treeList) {
  nTrees <- length(treeList)
  distances <- matrix(0, nTrees, nTrees)
  distTri <- unlist(sapply(1:(nTrees - 1), function (i) vapply((i+1):nTrees, function (j) {
    RF.dist(treeList[[i]], treeList[[j]])
  }, double(1))))
  distTri[distTri == 0] <- 1e-9
  distances[upper.tri(distances)] <- distTri
  distances[lower.tri(distances)] <- t(distances)[lower.tri(distances)] # Hat tip https://stackoverflow.com/questions/18165320/creating-a-symmetric-matrix-in-r
  distances
}


tntDirectories <- c('ambiguous', 'ambigAbsent', 'extraState')
rDirectories <- c('inapplicable')



nexusFiles <- list.files('matrices', pattern='.*\\.nex$')
nexusName <- nexusFiles[1]

trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName), 
              lapply(rDirectories, readRTrees, nexusName=nexusName))
nTrees <- vapply(trees, length, integer(1))
flatTrees <- unlist(trees, recursive=FALSE)
treeTitles <- paste(rep(c(tntDirectories, rDirectories), nTrees), unlist(sapply(nTrees, seq_len)))
treeCol <- paste(rep(treePalette, nTrees))
treePCh <- rep(plotChars, nTrees)

# Crude analysis inspired by http://www.kmeverson.org/blog/visualizing-tree-space-in-r

allDistances <- RFDistances(flatTrees)
allPcoa<-pcoa(allDistances)
plot(allPcoa$vectors[, 1], allPcoa$vectors[, 2], type = "p", xlab = "", ylab = "",
     axes = FALSE, main = paste(nexusName, "MDS tree space"), col=treeCol, pch=treePCh)
legend('top', legend=c(tntDirectories, rDirectories), cex = 1, pch = plotChars, col=treePalette)


# Really we want to something more sophisticated: e.g.
# HILLIS, D. M., HEATH, T. A., JOHN, K. St. and ANDERSON, F. 2005. Analysis and Visualization of Tree Space. Systematic Biology, 54, 471–482.
# or WILGENBUSCH, J. C., HUANG, W. and GALLIVAN, K. A. 2017. Visualizing phylogenetic tree landscapes. BMC Bioinformatics, 18, 85.
