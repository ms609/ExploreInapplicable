install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
library(phangorn)
devtools::install_github('ms609/inapplicable')
require(inapplicable)

nexusName <- 'Vinther2008.nex' # TODO foreach (file in inapplicable directory)
tntDirectories <- c('ambiguous', 'ambigAbsent', 'extraState')
rDirectories <- c('inapplicable')

readTntTrees <- function (directory, nexusName) {
  read.nexus(paste0(directory, '/', nexusName, '.nextrees', collapse=''))
}
readRTrees <- function (directory, nexusName) {
  allResults <- list.files(directory, paste0('\\-[[:digit:]]+.tre', collapse=''))
  resultScores <- vapply(allResults, function (string) {
    hits <- regexpr(pattern='\\-[[:digit:]]+', string)
    return(as.integer(substr(string, hits[1] + 1, hits[1] + attr(hits, 'match.length') - 1)))
  }, integer(1))
  unique(read.tree(paste0(directory, '/', allResults[which.min(resultScores)], collapse='')))
}

trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName), 
              lapply(rDirectories, readRTrees, nexusName=nexusName))

## Do something with trees