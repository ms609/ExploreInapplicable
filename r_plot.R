source('rPlot/functions.R')
source('rPlot/definitions.R')
OVERWRITE <- FALSE

# nexusName <- 'Asher2005.nex'

for (nexusName in nexusFiles) {
  par(mfrow = c(2, 2), bg = 'white')
  nexusRoot <- gsub('.nex', '', nexusName)
  cli::cli_h2(paste0("Evaluating ", nexusRoot, "..."))
  
  if (nexusRoot %in% c('avoidThisFile')) { # add slowFiles to this list if quick results wanted first
    cli::cli_alert_warning("Manual override")
    next
  }

  rTrees <- lapply(rDirectories, readRTrees, nexusName = nexusName)
  if (is.null(rTrees[[1]])) {
    cli::cli_alert_danger("R trees not found.")
    next
  }
  trees <- c(lapply(tntDirectories, readTntTrees, nexusName = nexusName), rTrees)
  nTrees <- vapply(trees, length, integer(1)); names(nTrees) <- allDirectories
  dirTrees <- TreeNumbers(nTrees)
  flatTrees <- unlist(trees, recursive=FALSE)
  
  treeSource <- rep(c(tntDirectories, rDirectories), nTrees)
  treeTitles <- paste(treeSource, unlist(sapply(nTrees, seq_len)), sep='_')
  treeCol <- paste(rep(treePalette, nTrees))
  treePCh <- rep(plotChars, nTrees)
  
  if (file.exists(paste0('histograms/landscape/', nexusRoot, '-300.png', collapse=''))  && !OVERWRITE) {
    cli::cli_alert_info("MPT histograms already exist.")
  } else {
    cli::cli_alert_info("Generating MPT histograms.")
    
    # Calculate tree scores
    scores <- GetTreeScores(nexusRoot, trees)
    minScores <- apply(scores, 2, min)
    extraSteps <- scores - matrix(minScores, nrow(scores), ncol(scores), byrow=TRUE)
    rm(scores)
    
    # Plot tree scores: 4x4
    cli::cli_alert("Plotting island scores...\n")
    for (dirPath in allDirectories) {
      dirScores <- extraSteps[dirTrees[[dirPath]], , drop=FALSE]
      dirBreaks <- -0.5:(max(dirScores) + 0.5)
      dirCol <- treePalette[dirPath]
      otherDirectories <- which(allDirectories != dirPath)
      
      yMax <- max(apply(dirScores[, otherDirectories, drop=FALSE], 2, function (x) max(table(x))))
      hist(0, breaks=dirBreaks, border='#ffffffff', ylim=c(0, yMax), axes=FALSE, font.main=1, cex.main=1,
           main=paste0("Trees on ", englishName[dirPath], ' island (', nexusRoot, ")"), col.main=dirCol, xlab="Steps longer than method's best tree") # Set up blank histogram
      axis(1, col=dirCol)
      axis(2, col=dirCol)
      
      for (i in otherDirectories) {
        hist(dirScores[, i, drop=FALSE], add=TRUE, breaks=dirBreaks,
             border=treePalette[i], col=paste0(treePalette[i], '99', collapse=''))
      }
    }
    dev.copy(svg, file=paste0('islandCounts4/', nexusRoot, '.svg', collapse='')); dev.off()
    dev.copy(png, file=paste0('islandCounts4/', nexusRoot, '.png', collapse=''), width=800, height=800); dev.off()
    
    
    # Plot tree scores: 3x3
    for (dirPath in allDirectories[-1]) {
      dirScores <- extraSteps[dirTrees[[dirPath]], , drop=FALSE]
      dirBreaks <- (-0.5):(max(dirScores) + 0.5)
      dirCol <- treePalette[dirPath]
      otherDirectories <- which(allDirectories != dirPath)[-1]
      
      yMax <- max(apply(dirScores[, otherDirectories, drop=FALSE], 2, function (x) max(table(x))))
      hist(0, breaks=dirBreaks, border='#ffffffff', ylim=c(0, yMax), axes=FALSE, font.main=1, cex.main=1,
           main=paste0("Trees on ", englishName[dirPath], ' island (', nexusRoot, ")"), col.main=dirCol, xlab="Steps longer than method's best tree") # Set up blank histogram
      axis(1, col=dirCol)
      axis(2, col=dirCol)
      
      for (i in otherDirectories) {
        hist(dirScores[, i, drop=FALSE], add=TRUE, breaks=dirBreaks,
             border=treePalette[i], col=paste0(treePalette[i], '99', collapse=''))
      }
    }
    rm (extraSteps)
    
    plot(-9, -9, axes=FALSE, , ylim=c(0, 100), xlim=c(0, 100), xlab='', ylab='')
    ySpace <- 12; yHeight = 7
    rect(10, ySpace * 1, 10 + yHeight - 1, yHeight + (ySpace * 1), col=paste0(treePalette[4], '99'), border=treePalette[4])
    rect(10, ySpace * 2, 10 + yHeight - 1, yHeight + (ySpace * 2), col=paste0(treePalette[3], '99'), border=treePalette[3])
    rect(10, ySpace * 3, 10 + yHeight - 1, yHeight + (ySpace * 3), col=paste0(treePalette[2], '99'), border=treePalette[2])
    text(20, ySpace * 1 + (yHeight / 2), englishName[4], pos=4)
    text(20, ySpace * 2 + (yHeight / 2), englishName[3], pos=4)
    text(20, ySpace * 3 + (yHeight / 2), englishName[2], pos=4)
    text(20, ySpace * 4 + (yHeight / 2), 'Method', pos=4, font=2)
    dev.copy(svg, file=paste0('islandCounts/', nexusRoot, '.svg', collapse='')); dev.off()
    dev.copy(png, file=paste0('islandCounts/', nexusRoot, '.png', collapse=''), width=800, height=800); dev.off()
  }
  
  par(mfrow=c(1, 1), bg='white')
  if (!file.exists(paste0('vennNodes/', nexusRoot, '.png')) || OVERWRITE) {
    cli::cli_alert("Calculating consensus trees")
    consensi <- lapply(trees, consensus)
    nNodes <- vapply(consensi, function (tr) tr$Nnode, double(1))
    cat("done.\n")
    vennCons <- integer(7)
    names(vennCons) <- c('A', 'B', 'C', 'A&B', 'A&C', 'B&C', 'A&B&C')
    vennCons[7] <- consensus(consensi[[4]], consensi[[3]], consensi[[2]])$Nnode
    vennCons[4] <- consensus(consensi[[2]], consensi[[3]])$Nnode - vennCons['A&B&C']
    vennCons[5] <- consensus(consensi[[2]], consensi[[4]])$Nnode - vennCons['A&B&C']
    vennCons[6] <- consensus(consensi[[4]], consensi[[3]])$Nnode - vennCons['A&B&C']
    vennCons[1:3] <- c(
    'A' =     nNodes[2] - (vennCons['A&B'] + vennCons['A&C'] + vennCons['A&B&C']),
    'B' =     nNodes[3] - (vennCons['A&B'] + vennCons['B&C'] + vennCons['A&B&C']),
    'C' =     nNodes[4] - (vennCons['A&C'] + vennCons['B&C'] + vennCons['A&B&C']))
    vennPlot <- venneuler::venneuler(vennCons)
    plot(vennPlot, col=treePalette[2:4], col.fn=function(x) x, border=treePalette[2:4], 
        col.txt=NA, edges=1024, main=paste0('Node presence by method: ', nexusRoot))
    dev.copy(svg, file=paste0('vennNodes/', nexusRoot, '.svg', collapse='')); dev.off()
    dev.copy(png, file=paste0('vennNodes/', nexusRoot, '.png', collapse=''), width=800, height=800); dev.off()
    cli::cli_alert_success("Plotted Venn diagram of nodes.")
  }
  
  if (!file.exists(paste0('vennTrees/', nexusRoot, '.png')) || OVERWRITE) {
    vennTrees <- integer(7)
    names(vennTrees) <- c('A', 'B', 'C', 'A&B', 'A&C', 'B&C', 'A&B&C')
    vennTrees[1:3] <- nTrees[2:4]
    vennTrees[4:7] <- vapply(list(trees[c(2, 3)], trees[c(2, 4)], trees[c(3, 4)], trees[2:4]),
      function (x) length(unique(unlist(x, recursive=FALSE))), integer(1))
    vennTrees[7] <- sum(vennTrees[1:3]) - vennTrees[7]
    vennTrees[4:6] <- vapply(3:1, function (i) sum(vennTrees[1:3][-i]), integer(1)) - vennTrees[4:6]
    vennPlot <- venneuler::venneuler(vennTrees)
    plot(vennPlot, col=treePalette[2:4], col.fn=function(x) x, border=treePalette[2:4],
        col.txt=NA, edges=1024, main=paste0('Shortest trees: ', nexusRoot))
    dev.copy(svg, file=paste0('vennTrees/', nexusRoot, '.svg', collapse='')); dev.off()
    dev.copy(png, file=paste0('vennTrees/', nexusRoot, '.png', collapse=''), width=800, height=800); dev.off()
    cli::cli_alert_success("Plotted Venn diagram of trees.")
  }
  par(mfrow=c(2, 2), bg='white')
  
  if (file.exists(paste0('treeSpaces/panels/', nexusRoot, '.svg', collapse='')) && !OVERWRITE) {
    cli::cli_alert_info("Treespace plots already exist.")
    next
  }
  
  charTypes <- vapply(readLines(paste0('charType/', nexusRoot, '.txt', collapse='')), substr, character(1), 1, 1, USE.NAMES=FALSE)
  # This will misestimate in the matrix where chars without inapps were coded N or T instead of being left as X.
  charSummary <- paste(names(table(charTypes)), table(charTypes), sep=': ', collapse='; ')
  
  
  rawData <- read.nexus.data(paste0('inapplicable/', nexusName, collapse = ''))
  qtTitleText <- paste(nexusRoot, "Quartet space\n",  length(rawData), 'taxa,', length(rawData[[1]]), 'chars -', charSummary)

  qtDistances <- GetQuartetDistances(nexusRoot, trees[-1], forPlot = TRUE)
  cli::cli_alert_success("Got Quartet distances.")
  
  ambigTrees <- seq_len(nTrees[1])

  par(mfrow = c(1, 1), mar = rep(0, 4))
  TreeSpacePanel(qtDistances, nTrees[-1], studyName[[nexusRoot]], 1.0)
  dev.copy(svg, file=paste0('treeSpaces/panels/', nexusRoot, '.svg', collapse=''), width=2, height=2); dev.off()
   
  
#  areaFile <- paste0('treeSpaces/', nexusRoot, '.hullAreas.csv')
#  write.csv(t(data.frame(rf.areas = rfAreas, qt.areas = qtAreas)), file=areaFile)
   
 
  par(mfrow=c(1, 1), bg='white')
  
#  PlotKruskalTreeSpace3(qtDistances[-ambigTrees, -ambigTrees], nTrees[-1], legendPos=QuartetLegendPos(nexusRoot), nexusRoot, fill=TRUE)
#  dev.copy(svg, file=paste0('quartetSpaces/', nexusRoot, '.svg', collapse='')); dev.off()
#  dev.copy(png, file=paste0('quartetSpaces/', nexusRoot, '.png', collapse=''), width=1024, height=1024); dev.off()
#  
  cat(" - Printed treespaces and saved to files.\nEvaluation complete.\n\n")
}
