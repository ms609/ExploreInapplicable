if (!require('ape')) install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
#devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
library(phangorn)
devtools::install_github('ms609/inapplicable')
require(inapplicable)
if (!require(rtqdist)) install.packages('http://users-cs.au.dk/cstorm/software/tqdist/files/tqDist-1.0.0.tar.gz', repos=NULL, type='source') # You can download it from http://users-cs.au.dk/cstorm/software/tqdist/

readTntTrees <- function (directory, nexusName) {
  unique(read.nexus(paste0(directory, '/', nexusName, '.nextrees', collapse='')))
}

readRTrees <- function (directory, nexusName) {
  allResults <- list.files(directory, paste0(nexusName, '.*\\-[[:digit:]]+.tre', collapse=''))
  if (length(allResults) == 0) return (NULL)
  resultScores <- vapply(allResults, function (string) {
    hits <- regexpr(pattern='\\-[[:digit:]]+', string)
    return(as.integer(substr(string, hits[1] + 1, hits[1] + attr(hits, 'match.length') - 1)))
  }, integer(1))
  bestTreeFile <- paste0(directory, '/', allResults[which.min(resultScores)], collapse='')
  unique(read.tree(bestTreeFile))
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

QuartetDistances <- function (treeList) {
  if (class(treeList) == 'list') class(treeList) <- 'multiPhylo'
  write.tree(treeList, file='~temp.trees')
  rtqdist::allPairsQuartetDistance('~temp.trees')
}

TreeNumbers <- function (nTrees) {
  res <- lapply(seq_along(allDirectories), function (i) {
    firstTree <- if (i == 1) 1 else cumsum(nTrees)[i-1] + 1
    lastTree  <- cumsum(nTrees)[i]
    iTrees <- firstTree:lastTree
  })
  names(res) <- allDirectories
  res
}

PlotTreeSpace <- function (pcs, nTrees, legendPos = 'bottomleft', mainTitle) {
  x <- pcs$vectors[, 1]
  y <- pcs$vectors[, 2]
  plot(x, y, type = "p", xlab = "", ylab = "",
       axes = FALSE, col=treeCol, pch=treePCh)
  title(main = mainTitle, cex.main=0.81)
  # Plot convex hulls
  iTrees <- TreeNumbers(nTrees)
  for (i in seq_along(nTrees)) {
    convexHull <- chull(x[iTrees[[i]]], y[iTrees[[i]]])
    convexHull <- c(convexHull, convexHull[1])
    lines(x[iTrees[[i]]][convexHull], y[iTrees[[i]]][convexHull], col=treePalette[i])
  }
  
  legend(legendPos, bty='n', legend=paste0(c(tntDirectories, rDirectories), ' (', nTrees, ')'),
         cex = 0.75, pch = plotChars, col=treePalette)
}

PlotTreeSpace3 <- function (pcs, nTrees, legendPos = 'bottomleft', mainTitle) {
  x <- pcs$vectors[, 1]
  y <- pcs$vectors[, 2]
  ambigTrees <- seq_len(nTrees[[1]])
  plot(x, y, type = "p", xlab = "", ylab = "",
       axes = FALSE, col=treeCol[-ambigTrees], pch=treePCh[-ambigTrees])
  title(main = mainTitle, cex.main=0.81)
  # Plot convex hulls
  iTrees <- TreeNumbers(nTrees)
  for (i in seq_along(nTrees)[-1]) {
    convexHull <- chull(x[iTrees[[i]] - nTrees[1]], y[iTrees[[i]] - nTrees[1]])
    convexHull <- c(convexHull, convexHull[1])
    lines(x[iTrees[[i]] - nTrees[1]][convexHull], y[iTrees[[i]] - nTrees[1]][convexHull], col=treePalette[i])
  }
  
  legend(legendPos, bty='n', legend=paste0(c('Ambiguous', 'Extra state', 'Inapplicable'), ' (', nTrees[-1], ')'),
         cex = 0.75, pch = plotChars[-1], col=treePalette[-1])
}
