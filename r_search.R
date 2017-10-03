install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
library(phangorn)
devtools::install_github('ms609/inapplicable')
require(inapplicable)

inappFiles <- list.files('inapplicable', pattern='.*\\.nex$')
filename <- 'Geisler2001.nex'#inappFiles[2]

rawData <- read.nexus.data(paste0('inapplicable/', filename, collapse=''))
phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
phyData
best <- ape::root(ape::nj(phangorn::dist.hamming(phyData)), names(rawData)[1], resolve.root=TRUE)
attr(best, 'pscore') <- 1e+7
bestScore <- 1e+7
bestHits <- 0


maxHitMenu <- c(2:16, 25, 40, 60, 100, 200)


for (i in 1:10000) {
  started <- Sys.time()
  cat ("\n > Setting maxHits =", maxHits <- sample(maxHitMenu, 1))
  best <- Ratchet(best, phyData, maxIt=1e+7, maxIter=3e+7, maxHits=maxHits, k = 1, verbosity=0)
  secsTaken <- as.numeric(difftime(Sys.time(), started, units='secs'))
  plot (best, cex=0.95, main=attr(best, 'pscore'))
  treeScore <- attr(best, 'pscore')
  if (treeScore < bestScore) {
    bestScore <- treeScore
    bestHits <- 0
  }
  bestHits <- bestHits + 1
  resultsFile <- paste0('inapplicable/', filename, '-', treeScore, '.tre')
  write.tree(best, resultsFile, append=file.exists(resultsFile))
  write.table(data.frame(maxHits, secsTaken, filename, treeScore), file="inapplicable/searchTimes.csv", sep=',', col.names=FALSE, row.names=FALSE, append=TRUE)
  if (bestHits >= 100) break;
  cat ("; found best score", treeScore, "in", secsTaken, 's; hit', bestHits)
}

