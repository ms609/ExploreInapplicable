install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
library(phangorn)
devtools::install_github('ms609/inapplicable')
require(inapplicable)

filename <- 'Dikov2009.nex'#inappFiles[2]


inappFiles <- list.files('inapplicable', pattern='.*\\.nex$')
#filename <- 'OMeara2013.nex'

for (filename in inappFiles) {
  cat(" - loading from", filename, "\n")
  
  results <- list.files('inapplicable', pattern=paste0(filename, '.*\\.tre$'))
  if (any(results)) {
    cat ("Results exist for", filename, "; next file please.\n\n")
    next
  }
  rawData <- read.nexus.data(paste0('inapplicable/', filename, collapse=''))
  phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
    

  best <- ape::root(ape::nj(phangorn::dist.hamming(phyData)), names(rawData)[1], resolve.root=TRUE)
  attr(best, 'pscore') <- 1e+7
  bestScore <- 1e+7
  bestHits <- 0
  maxHitMenu <- c(2:16, 25, 40, 60, 100, 200)

  verbose <- 0
  for (i in 1:15000) {
    started <- Sys.time()
    cat ("\n > Setting maxHits =", maxHits <- sample(maxHitMenu, 1))
    best <- Ratchet(best, phyData, maxIt=1e+7, maxIter=5e+4, maxHits=maxHits, k = 1, verbosity=verbose) # Maxiter of 3e+7 was causing unhappy delays when there was a single best tree at a local optimum.
    secsTaken <- as.numeric(difftime(Sys.time(), started, units='secs'))
    #verbose <- if (secsTaken > 5) 3 else if (secsTaken > 3) 1 else if (secsTaken < 2) 0 else verbose
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
    cat ("; found best score", treeScore, "in", secsTaken, 's; hit', bestHits)
    if (bestHits >= 100) break;
  }
  cat("\n\n ----- \n\n > Completed" , filename, "\n")
}
