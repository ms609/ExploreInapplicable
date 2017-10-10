source('rPlot/functions.R')
source('rPlot/definitions.R')
OVERWRITE <- FALSE

nexusFiles <- list.files('matrices', pattern='.*\\.nex$');# nexusFiles
#nexusName <- 'Eklund2004.nex'
slowFiles <- c('Aria2015', 'Aguado2009', 'Capa2011', 'Conrad2008', 'Dikow2009', 'Eklund2004', 
               'Geisler2001', 'Giles2015', 'OMeara2014', 'Wills2012', 'Schulze2007')

for (nexusName in nexusFiles) {
  par(mfrow=c(2, 2), bg='white')
  nexusRoot <- gsub('.nex', '', nexusName);
  cat("\nEvaluating", nexusRoot, "...\n")
  
  if (nexusRoot %in% c('avoidThisFile.nex')) { # add slowFiles to this list if quick results wanted first
    cat (" ! Manual override\n")
    next
  }

  
  rTrees <- lapply(rDirectories, readRTrees, nexusName=nexusName)
  if (is.null(rTrees[[1]])) {
    cat(" ! R trees not found.\n")
    next
  }
  trees <- c(lapply(tntDirectories, readTntTrees, nexusName=nexusName), rTrees)
  cat(" - Trees read OK.\n")
  nTrees <- vapply(trees, length, integer(1)); names(nTrees) <- allDirectories
  dirTrees <- TreeNumbers(nTrees)
  flatTrees <- unlist(trees, recursive=FALSE)
  
  treeSource <- rep(c(tntDirectories, rDirectories), nTrees)
  treeTitles <- paste(treeSource, unlist(sapply(nTrees, seq_len)), sep='_')
  treeCol <- paste(rep(treePalette, nTrees))
  treePCh <- rep(plotChars, nTrees)
  
  if (file.exists(paste0('islandCounts/', nexusRoot, '.png', collapse=''))  && !OVERWRITE) {
    cat(" - MPT histograms already exist.\n")    
  } else {
    
    # Calculate tree scores
    treeScoreFile <- paste0('islandCounts/', nexusRoot, '.csv')
    if (file.exists(treeScoreFile)) {
      scores <- data.matrix(read.csv(treeScoreFile))
    } else {      
      cat(" - Calculating tree scores...\n")
      scores <- vapply(allDirectories, function (dirPath) {
        TreeScorer <- if (dirPath %in% tntDirectories) phangorn::fitch else inapplicable::InapplicableFitch
        rawData <- read.nexus.data(paste0(dirPath, '/', nexusName, collapse=''))
        phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
        as.integer(vapply(flatTrees, TreeScorer, double(1), phyData, USE.NAMES=FALSE))
      }, integer(sum(nTrees)))
      rownames(scores) <- treeSource
      write.csv(scores, file=treeScoreFile)
    }
    minScores <- apply(scores, 2, min)
    extraSteps <- scores - matrix(minScores, nrow(scores), ncol(scores), byrow=TRUE)
    rm(scores)
    
    # Plot tree scores: 4x4
    cat(" - Plotting island scores...\n")
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
    dev.copy(svg, file=paste0('islandCounts/', nexusRoot, '.svg', collapse='')); dev.off()
    dev.copy(png, file=paste0('islandCounts/', nexusRoot, '.png', collapse=''), width=800, height=800); dev.off()
    
    
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
    dev.copy(svg, file=paste0('islandCounts/', nexusRoot, '-3.svg', collapse='')); dev.off()
    dev.copy(png, file=paste0('islandCounts/', nexusRoot, '-3.png', collapse=''), width=800, height=800); dev.off()
  }
  
  par(mfrow=c(1, 1), bg='white')
  if (!file.exists(paste0('vennNodes/', nexusRoot, '.png'))|| OVERWRITE) {
    cat(" - Calculating consensus trees: ")
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
    vennPlot <- venneuler(vennCons)
    plot(vennPlot, col=treePalette[2:4], col.fn=function(x) x, border=treePalette[2:4], col.txt=NA,
         main=paste0('Nodes recovered by method: ', nexusRoot))
    dev.copy(svg, file=paste0('vennNodes/', nexusRoot, '.svg', collapse='')); dev.off()
    dev.copy(png, file=paste0('vennNodes/', nexusRoot, '.png', collapse=''), width=512, height=512); dev.off()
    cat(" - Plotted Venn diagram of nodes.\n")
  }
  
  if (!file.exists(paste0('vennTrees/', nexusRoot, '.png')) || OVERWRITE) {
    vennTrees <- integer(7)
    names(vennTrees) <- c('A', 'B', 'C', 'A&B', 'A&C', 'B&C', 'A&B&C')
    vennTrees[1:3] <- nTrees[2:4]
    vennTrees[4:7] <- vapply(list(trees[c(2, 3)], trees[c(2, 4)], trees[c(3, 4)], trees[2:4]),
      function (x) length(unique(unlist(x, recursive=FALSE))), integer(1))
    vennTrees[7] <- sum(vennTrees[1:3]) - vennTrees[7]
    vennTrees[4:6] <- vapply(3:1, function (i) sum(vennTrees[1:3][-i]), integer(1)) - vennTrees[4:6]
    vennPlot <- venneuler(vennTrees)
    plot(vennPlot, col=treePalette[2:4], col.fn=function(x) x, border=treePalette[2:4], col.txt=NA,
         main=paste0('Shortest trees: ', nexusRoot))
    dev.copy(svg, file=paste0('vennTrees/', nexusRoot, '.svg', collapse='')); dev.off()
    dev.copy(png, file=paste0('vennTrees/', nexusRoot, '.png', collapse=''), width=512, height=512); dev.off()
    cat(" - Plotted Venn diagram of trees.\n")
  }
  par(mfrow=c(2, 2), bg='white')
  
  if (file.exists(paste0('treeSpaces/', nexusRoot, '.png', collapse='')) && !OVERWRITE) {
    cat(" - Treespace plots already exist.\n")
    next
  }
  charTypes <- vapply(readLines(paste0('charType/', nexusRoot, '.txt', collapse='')), substr, character(1), 1, 1, USE.NAMES=FALSE)
  charSummary <- paste(names(table(charTypes)), table(charTypes), sep=': ', collapse='; ')
  
  rawData <- read.nexus.data(paste0('inapplicable/', nexusName, collapse=''))
  rfTitleText <- paste(nexusRoot, "R-F space\n",  length(rawData), 'taxa,', length(rawData[[1]]), 'chars -', charSummary)
  qtTitleText <- paste(nexusRoot, "Quartet space\n",  length(rawData), 'taxa,', length(rawData[[1]]), 'chars -', charSummary)

  # Crude analysis inspired by http://www.kmeverson.org/blog/visualizing-tree-space-in-r
  rfFileName <- paste0('treeSpaces/', nexusRoot, '.rf.csv')
  qtFileName <- paste0('treeSpaces/', nexusRoot, '.qt.csv')
  if (file.exists(rfFileName)) {
    rfDistances <- data.matrix(read.csv(rfFileName, row.names=1))
  } else {
    cat(" - Calculating RF distances...\n")
    rfDistances <- RFDistances(flatTrees)
    write.csv(rfDistances, file=rfFileName)
  }
  if (file.exists(qtFileName)) {
    qtDistances <- data.matrix(read.csv(qtFileName, row.names=1))
  } else {
    cat(" - Calculating quartet distances...\n")
    qtDistances <- QuartetDistances(flatTrees)
    write.csv(qtDistances, file=qtFileName)
  }
  rm(flatTrees)
  ambiguousTrees <- 1:nTrees[1]
  
  PlotKruskalTreeSpace (rfDistances, nTrees, legendPos='bottomright', rfTitleText)
  rfAreas <- PlotKruskalTreeSpace3(rfDistances, nTrees, legendPos='bottomright', rfTitleText)
  cat(" - Printed RF treespace.\n")
  qtDistances[qtDistances == 0] <- 1e-9
  diag(qtDistances) <- 0
  PlotKruskalTreeSpace (qtDistances, nTrees, legendPos='bottomright', qtTitleText)
  qtAreas <- PlotKruskalTreeSpace3(qtDistances, nTrees, legendPos='bottomright', qtTitleText)
  
  areaFile <- paste0('treeSpaces/', nexusRoot, '.hullAreas.csv')
  write.csv(t(data.frame(rf.areas = rfAreas, qt.areas = qtAreas)), file=areaFile)
  
  
  #rfSpace <- modifiedPcoa(rfDistances, correction='none')
  #PlotTreeSpace(rfSpace, nTrees, legendPos='bottomright', rfTitleText)
  #rf3 <- modifiedPcoa(rfDistances[-ambiguousTrees, -ambiguousTrees], correction='none')
  #PlotTreeSpace3(rf3, nTrees, legendPos='bottomright', rfTitleText)
  # If you want to understand what's going on, try
  # PlotTreeSpace3D(rf3, nTrees, legendPos='bottomright', rfTitleText)

  # qtSpace <- modifiedPcoa(qtDistances, correction='none')
  # PlotTreeSpace(qtSpace, nTrees, legendPos='bottomleft', qtTitleText)
  # qt3 <- modifiedPcoa(qtDistances[-ambiguousTrees, -ambiguousTrees], correction='none')
  # PlotTreeSpace3(qt3, nTrees, legendPos='bottomleft', qtTitleText)
  # If you want to understand what's going on, try
  # PlotTreeSpace3D(qt3, nTrees, legendPos='bottomright', rfTitleText)
  
  dev.copy(svg, file=paste0('treeSpaces/', nexusRoot, '.svg', collapse='')); dev.off()
  dev.copy(png, file=paste0('treeSpaces/', nexusRoot, '.png', collapse=''), width=1024, height=1024); dev.off()
  cat(" - Printed Quartet treespace and saved to files.\nEvaluation complete.\n\n")
}
