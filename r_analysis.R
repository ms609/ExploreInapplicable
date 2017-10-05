source('rPlot/functions.R')
source('rPlot/definitions.R')

### Analyse precision:
### Island sizes
csvs <- list.files('islandCounts', pattern='.*\\.csv$');# nexusFiles

islandSizes <- vapply(csvs, function (csvFile) {
  content <- read.csv(csvFile)
  table(content$X)
}, double(4))[-1, ]

ncol(islandSizes)
sum(islandSizes[3, ] == islandSizes[2, ])
sum(islandSizes[3, ] == islandSizes[1, ])
sum(islandSizes[2, ] == islandSizes[1, ])
sum(islandSizes[3, ] > islandSizes[2, ])
sum(islandSizes[3, ] > islandSizes[1, ])
sum(islandSizes[2, ] > islandSizes[1, ])

nexusFiles <- list.files('matrices', pattern='.*\\.nex$');# nexusFiles

precision <- matrix(NA, 4, length(nexusFiles))
colnames(precision) <- nexusFiles; rownames(precision) <- allDirectories
differences <- matrix(NA, 4, length(nexusFiles))
colnames(differences) <- nexusFiles; rownames(differences) <- c('diffAmbExt', 'diffAmbInap', 'diffExtrInap', 'allCons')
for (nexusName in nexusFiles) {
  nexusRoot <- gsub('.nex', '', nexusName); 
  cat("\nEvaluating", nexusRoot, "...\n")
 
  rTrees <- lapply(rDirectories, readRTrees, nexusName=nexusName)
  if (is.null(rTrees[[1]])) {
    cat(" ! R trees not found.\n")
    next
  }
  trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName), rTrees)
  consensi <- lapply(trees, consensus)
  
  cons23 <- consensus(consensi[[2]], consensi[[3]])$Nnode
  cons24 <- consensus(consensi[[2]], consensi[[4]])$Nnode
  cons34 <- consensus(consensi[[4]], consensi[[3]])$Nnode
  allCons <- consensus(consensi[2:4])$Nnode
  differences[, nexusName] <- c(cons23, cons24, cons34, allCons)
  
  nNodes <- vapply(consensi, function (tr) tr$Nnode, double(1))
  precision[, nexusName] <- nNodes
}

precisions <- precision[-1, !is.na(precision[1, ])]

ncol(precisions)
sum(precisions[3, ] == precisions[2, ])
sum(precisions[3, ] == precisions[1, ])
sum(precisions[2, ] == precisions[1, ])
sum(precisions[3, ] >  precisions[2, ])
sum(precisions[3, ] >  precisions[1, ])
sum(precisions[2, ] >  precisions[1, ])

agreements <- differences[, !is.na(differences[1, ])]
agreements

nodesLostByAddingAmbig <- (precisions['inapplicable', ] - agreements['diffAmbInap', ] ) /  precisions['inapplicable', ]
nodesLostByAddingExtSt <- (precisions['inapplicable', ] - agreements['diffExtrInap', ]) / precisions['inapplicable', ] 

1 - mean(nodesLostByAddingAmbig)
1 - mean(nodesLostByAddingExtSt)

ambigNodesListByAddingExtSt <- (precisions['ambigAbsent', ] - agreements['diffAmbExt', ]) / precisions['ambigAbsent', ] 
1-mean(ambigNodesListByAddingExtSt)
