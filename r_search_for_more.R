install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
library(phangorn)
if (!require(TreeSearch)) devtools::install_github('ms609/TreeSearch')
require(TreeSearch)

TARGET <- 250

inappFiles <- list.files('inapplicable', pattern='.*\\.nex$')

for (filename in inappFiles) {
  cat(" - loading from", filename, "\n")

  results <- list.files('inapplicable', pattern=paste0(filename, '.*\\-[[:digit:]]+.tre$'))
  if (length(results) == 0) {
    cat (" > No results yet exist.\n")
    next
  }

  resultScores <- vapply(results, function (string) {
    hits <- regexpr(pattern='\\-[[:digit:]]+', string)
    return(as.integer(substr(string, hits[1] + 1, hits[1] + attr(hits, 'match.length') - 1)))
  }, integer(1))
  bestTreeFile <- paste0('inapplicable/', results[which.min(resultScores)], collapse='')
  bestHits <- length(readLines(bestTreeFile))
  if (bestHits >= TARGET) {
    cat (" > Already have TARGET hits.\n")
    next
  }

  rawData <- read.nexus.data(paste0('inapplicable/', filename, collapse=''))
  phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))

  best <- ape::read.tree(text=readLines(bestTreeFile)[[1]])
  bestScore <- attr(best, 'pscore') <- inapplicable::InapplicableFitch(best, phyData)

  maxHitMenu <- 2:20 # Small maxHits finds result faster, but chance of improving score per
                     # second spent increases with maxHits at least up to 150. (Data for higher not available.)

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
    #write.table(data.frame(maxHits, secsTaken, filename, treeScore), file="inapplicable/searchTimes.csv", sep=',', col.names=FALSE, row.names=FALSE, append=TRUE)
    cat ("; found best score", treeScore, "in", secsTaken, 's; hit', bestHits)
    if (bestHits >= TARGET) break;
  }
  cat("\n\n ----- \n\n > Completed" , filename, "\n")
}
