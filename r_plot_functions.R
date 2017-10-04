#devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
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

PlotTreeSpace <- function (pcs, nTrees, legendPos = 'bottomleft') {
  x <- pcs$vectors[, 1]
  y <- pcs$vectors[, 2]
  plot(x, y, type = "p", xlab = "", ylab = "",
       axes = FALSE, main = paste(nexusName, "MDS tree space"), col=treeCol, pch=treePCh)
  
  # Plot convex hulls
  for (i in seq_along(nTrees)) {
    firstTree = if (i == 1) 1 else cumsum(nTrees)[i-1] + 1
    lastTree = cumsum(nTrees)[i]
    iTrees <- firstTree:lastTree
    convexHull <- chull(x[iTrees], y[iTrees])
    convexHull <- c(convexHull, convexHull[1])
    lines(x[iTrees][convexHull], y[iTrees][convexHull], col=treePalette[i])
  }
  
  legend(legendPos, bty='n', legend=paste0(c(tntDirectories, rDirectories), ' (', nTrees, ')'),
         cex = 0.75, pch = plotChars, col=treePalette)
}
