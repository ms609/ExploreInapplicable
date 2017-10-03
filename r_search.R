install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
library(phangorn)
devtools::install_github('ms609/inapplicable')
require(inapplicable)

filename <- 'Vinther2008.nex'
rawData <- read.nexus.data(paste0('matrices/', filename, collapse=''))
phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
if (!exists("best")) best <- ape::root(ape::nj(phangorn::dist.hamming(phyData)), names(rawData)[1], resolve.root=TRUE)
attr(best, 'pscore') <- 1e+7

for (i in 1:10000) {
  started <- Sys.time()
  cat ("\n > Setting maxHits =", maxHits <- sample(c(2:16, 25, 40, 60, 100, 200), 1))
  best <- Ratchet(best, phyData, maxIt=1e+7, maxIter=3e+7, maxHits=maxHits, k = 1, verbosity=0)
  secsTaken <- as.numeric(difftime(Sys.time(), started, units='secs'))
  plot (best, cex=0.95, main=attr(best, 'pscore'))
  treeScore <- attr(best, 'pscore')
  resultsFile <- paste0('inapplicable/', filename, '-', treeScore, '.tre')
  write.tree(best, resultsFile, append=file.exists(resultsFile))
  write.table(data.frame(maxHits, secsTaken, filename, treeScore), file="inapplicable/searchTimes.csv", sep=',', col.names=FALSE, row.names=FALSE, append=TRUE)
  cat ("; found best score", treeScore, "in", secsTaken, 's')
}

