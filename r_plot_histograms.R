source('rPlot/functions.R')
source('rPlot/definitions.R')
OVERWRITE <- FALSE

# nexusName <- 'Vinther2008.nex'
figuredDatasets <- c('DeAssis2011.nex', 'Asher2005.nex', 'Wetterer2000.nex', 'Vinther2008.nex')

for (nexusName in nexusFiles) {
  # par(mfrow=c(1, 1), bg='white')   #par(mfrow=c(3, 1), bg='white') for portrait, for example
  nexusRoot <- gsub('.nex', '', nexusName)
  cat("\nEvaluating", nexusRoot, "...\n")
  trees <- GetTrees(nexusRoot)
  cat(" - Trees read OK.\n")
  nTrees <- vapply(trees, length, integer(1)); names(nTrees) <- allDirectories
  dirTrees <- TreeNumbers(nTrees)
  flatTrees <- unlist(trees, recursive=FALSE)
  
  treeSource <- rep(c(tntDirectories, rDirectories), nTrees)
  treeTitles <- paste(treeSource, unlist(sapply(nTrees, seq_len)), sep='_')
  treeCol <- paste(rep(treePalette, nTrees))
  treePCh <- rep(plotChars, nTrees)
  
  if (file.exists(paste0('histograms/', nexusRoot, '-inapplicable-5in.pdf', collapse=''))  && !OVERWRITE) {
    cat(" - MPT histograms already exist.\n")    
  } else {
    cat(" - Generating MPT histograms.\n")    
    
    # Calculate tree scores
    scores <- GetTreeScores(nexusRoot, trees)
    minScores <- apply(scores, 2, min)
    maxScores <- apply(scores, 2, max)
    extraSteps <- scores - matrix(minScores, nrow(scores), ncol(scores), byrow=TRUE)
    
    # Plot newHists
    for (dirPath in allDirectories[-1]) {
      bestScore <- minScores[dirPath]
      dirScores <- extraSteps[dirTrees[[dirPath]], , drop=FALSE]
      
      if (max(maxScores - minScores) > 28) {
        dirBreaks <- ((-1:43) * 3) + 0.5 # Highest extra score is 129
      } else if (max(maxScores - minScores) > 8) {
        dirBreaks <- (-1:28) + 0.5
      } else {
        dirBreaks <- (-1:8) + 0.5
      }
      cat("   > ", dirPath, ": Max score", max(dirScores), "; max break", round(max(dirBreaks), 3), "\n")
      dirCol <- treePalette[dirPath]
      otherDirectories <- which(allDirectories != dirPath)[-1]
      
      yMax <- nrow(dirScores)# max(apply(dirScores[, otherDirectories, drop=FALSE], 2, function (x) max(table(x))))
      hist(0, breaks=dirBreaks, border='#ffffffff', ylim=c(0, yMax), axes=FALSE,            
           #main=paste0("Trees on ", englishName[dirPath], ' island (', nexusRoot, ")"), col.main=dirCol, cex.main=1.5, font.main=1, 
           main='', ann=FALSE,
           xlab="steps longer than method's best tree") # Set up blank histogram
      axis(1, col=dirCol, line=-0.1); mtext(side=1, text="Steps longer than method's best tree", line=2.35-0.1)
      axis(2, col=dirCol, line=-0.5); mtext(side=2, text='Frequency', line=2.35-0.5)
      text(max(dirBreaks), yMax * 0.8, paste('Best score =', minScores[dirPath]), cex=1.5, 
           pos=2, col=dirCol)
      
      for (i in otherDirectories) {
        hist(dirScores[, i, drop=FALSE], add=TRUE, breaks=dirBreaks,
             border=treePalette[i], col=paste0(treePalette[i], '99', collapse=''))
      }
      dev.copy(pdf, file=paste0('histograms/', nexusRoot, '-', dirPath, '-5in.pdf', collapse=''), width=5, height=5); dev.off()
    }
    
    #dev.copy(svg, file=paste0('histograms/portrait/', nexusRoot, '.svg', collapse='')); dev.off()
    #dev.copy(png, file=paste0('histograms/', nexusRoot, '-200.png', collapse=''), width=600, height=200); dev.off()
    #dev.copy(png, file=paste0('histograms/landscape/', nexusRoot, '-300.png', collapse=''), width=900, height=300); dev.off()
    #dev.copy(pdf, file=paste0('histograms/portrait/', nexusRoot, '-300.pdf', collapse=''), width=3, height=9); dev.off()
    #dev.copy(png, file=paste0('histograms/portrait/', nexusRoot, '-300.png', collapse=''), width=300, height=900); dev.off()
  }
}

