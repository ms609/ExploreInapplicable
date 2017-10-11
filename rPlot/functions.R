if (!require('ape')) install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
#devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
require(phangorn)
devtools::install_github('ms609/inapplicable')
require(inapplicable)
library(MASS)
library(venneuler)

if (!require(rtqdist)) install.packages('http://users-cs.au.dk/cstorm/software/tqdist/files/tqDist-1.0.0.tar.gz', repos=NULL, type='source') # You can download it from http://users-cs.au.dk/cstorm/software/tqdist/

readTntTrees <- function (directory, nexusName) {
  cat ("   - Reading directory", directory, '... ')
  treeLines <- readLines(paste0(directory, '/', nexusName, '.nextrees', collapse=''), warn=FALSE)
  treeLinesLight <- gsub('tree\\d*|\\s+', '', treeLines)
  if (!all(treeLinesLight %in% unique(treeLinesLight))) stop ("Non-unique trees found")
  cat(" all trees unique. ")
  treeList <- read.nexus(paste0(directory, '/', nexusName, '.nextrees', collapse=''))
  if (class(treeList) == 'phylo') {
    treeList <- list(treeList)
    class(treeList) <- 'multiPhylo'
    cat("Single tree returned.\n")
    return(treeList)
  }
  cat ("Read", length(treeList), "trees.\n")
  treeList
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
  fileName <- paste0('~temp', substring(runif(1), 3), '.trees')
  write.tree(treeList, file=fileName)
  on.exit(file.remove(fileName))
  rtqdist::allPairsQuartetDistance(fileName)
}

TreeNumbers <- function (nTrees) {
  if (length(nTrees) == 3) nTrees <- c(0, nTrees)
  res <- lapply(seq_along(allDirectories), function (i) {
    firstTree <- if (i == 1) 1 else cumsum(nTrees)[i - 1] + 1
    lastTree  <- cumsum(nTrees)[i]
    iTrees <- firstTree:lastTree
  })
  names(res) <- allDirectories
  if (res[[1]][2] == 0) res[[1]] <- NULL
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
  
  legend(legendPos, bty='n', legend=paste0(englishName, ' (', nTrees, ')'),
         cex = 0.75, pch = plotChars, col=treePalette)
}

PlotTreeSpace3 <- function (pcs, nTrees, legendPos = 'bottomleft', mainTitle) {
  pts <- pcs$vectors
  x <- pts[, 1]
  y <- pts[, 2]
 
  ambigTrees <- seq_len(nTrees[[1]])
  plot(x, y, type = "p", xlab = "", ylab = "",
       axes = FALSE, col=treeCol[-ambigTrees], pch=treePCh[-ambigTrees])
  title(main = mainTitle, cex.main=0.81)
  
  iTrees <- TreeNumbers(nTrees)
  for (i in seq_along(nTrees)[-1]) {
    #xyz <- pts[TreeNumbers(nTrees)[[i]] - nTrees[1], 1:3, drop=FALSE]
    #hull <- geometry::convhulln(xyz, option='FA')
    #polygon(xyz[t(hull$hull), 1:2], border=paste0(treePalette[i], '33'), lty=1)
    convexHull <- chull(x[iTrees[[i]] - nTrees[1]], y[iTrees[[i]] - nTrees[1]])
    convexHull <- c(convexHull, convexHull[1])
    lines(x[iTrees[[i]] - nTrees[1]][convexHull], y[iTrees[[i]] - nTrees[1]][convexHull], col=treePalette[i])
  }
  
  legend(legendPos, bty='n', legend=paste0(c('Ambiguous', 'Extra state', 'Inapplicable'), ' (', nTrees[-1], ')'),
         cex = 0.75, pch = plotChars[-1], col=treePalette[-1])
}

PlotKruskalTreeSpace <- function (distances, nTrees, legendPos = 'bottomleft', mainTitle) {
  
  scaled <- MASS::isoMDS(distances, k = 2, trace=FALSE)$points # Kruskal's non-multimetric MDS
  x <- scaled[, 1]
  y <- scaled[, 2]
  plot(scaled, type = "p", xlab = "", ylab = "", axes = FALSE, col=treeCol, pch=treePCh)
  title(main = mainTitle, cex.main=0.81)
  # Plot convex hulls
  iTrees <- TreeNumbers(nTrees)
  for (i in seq_along(nTrees)) {
    convexHull <- chull(x[iTrees[[i]]], y[iTrees[[i]]])
    convexHull <- c(convexHull, convexHull[1])
    lines(x[iTrees[[i]]][convexHull], y[iTrees[[i]]][convexHull], col=treePalette[i])
  }
  
  legend(legendPos, bty='n', legend=paste0(englishName, ' (', nTrees, ')'),
         cex = 0.75, pch = plotChars, col=treePalette)
}

PlotKruskalTreeSpace3 <- function (distances, nTrees, legendPos = 'bottomleft', mainTitle=NULL, fill=FALSE, ...) {
  if (length(nTrees) == 4) {
    warning("You should only send three trees here.")
  } else {
    scaled <- MASS::isoMDS(distances, k = 2, trace=FALSE)$points # Kruskal's non-multimetric MDS
  }
  x <- scaled[, 1]
  y <- scaled[, 2]
  if (!is.null(mainTitle)) {
    plotCol <- rep(treePalette[2:4], nTrees)
    plotPCh <- rep(plotChars[2:4], nTrees)
    plot(x, y, type = "p", xlab = "", ylab = "",
         axes = FALSE, col=plotCol, pch=plotPCh, 
         main = mainTitle, cex.main=0.81,
         ...)
  }
  
  iTrees <- TreeNumbers(nTrees)

  hullArea <- double(3)
  names(hullArea) <- allDirectories[-1]
  for (i in seq_along(nTrees)) {
    convexHull <- chull(x[iTrees[[i]]], y[iTrees[[i]]])
    convexHull <- c(convexHull, convexHull[1])
    convX <- x[iTrees[[i]]][convexHull]
    convY <- y[iTrees[[i]]][convexHull]
    if (fill) {
      polygon(convX, convY, col=paste0(treePalette[i + 1], '4B'), border=treePalette[i + 1]) #4B = 30% alpha
    } else {
      if (!is.null(mainTitle)) lines(convX, convY, col=treePalette[i + 1])
    }
    hullArea[i] <- geometry::polyarea(convX, convY)
  }
  
  if (!is.null(mainTitle)) {
    legend(legendPos, bty='n', legend=paste0(c('Ambiguous', 'Extra state', 'Inapplicable'), ' (', nTrees, ')'),
           cex = 0.75, pch = plotChars[-1], col=treePalette[-1])
  }
  
  c(hullArea / hullArea['inapplicable'])
}


PlotTreeSpace3D <- function (pcs, nTrees, legendPos = 'bottomleft', mainTitle) {
  pts <- pcs$vectors
  
  plot3d(pts, col=treeCol[-(1:nTrees[1])], size=5, axes=FALSE, xlab='', ylab='', zlab='')
  
  #rgl.open()
  #rgl.bg(color='white')
  #rgl.points(pts, color=treeCol[-(1:nTrees[1])], size=5)
  #rgl.spheres(pts, color=treeCol[-(1:nTrees[1])], radius=0.2)
  limits <- function (x) c(-max(abs(x)), max(abs(x))) * 1.1
  rgl.lines(limits(pcs$vectors[, 1]), c(0, 0), c(0, 0), color = "blue")
  rgl.lines(c(0, 0), limits(pcs$vectors[, 2]), c(0, 0), color = "red")
  rgl.lines(c(0, 0), c(0, 0), limits(pcs$vectors[, 3]), color = "green")
  
  for (i in 2:4) {
    xyz <- pts[TreeNumbers(nTrees)[[i]] - nTrees[1], 1:3, drop=FALSE]
    hull <- geometry::convhulln(xyz, option='FA')
    triangles3d(xyz[t(hull$hull), ], col=treePalette[i], alpha=0.3)
  }
  
  aspect3d('iso')
  
}

MatrixProperties <- function (fileRoot) {
  fileName <- paste0('inapplicable', '/', fileRoot, '.nex', collapse='')
  if (!file.exists(fileName)) return (integer(6))
  rawData <- read.nexus.data(fileName)
  inapplicableTokens <- vapply(rawData, function (x) x == '-', logical(length(rawData[[1]])))
  inapplicableChars <- rowSums(inapplicableTokens) > 0
  ambiguousTokens <- vapply(rawData, function (x) x == '?', logical(length(rawData[[1]])))
  ret <-    c(nTax = length(rawData),
              nChar = length(rawData[[1]]),
              nInapp = sum(inapplicableChars),
              #whichInapp = inapplicableChars,
              nTokens = length(rawData) * length(rawData[[1]]),
              inappTokens = sum(inapplicableTokens),
              ambigTokens = sum(ambiguousTokens)
              )  
}

GetTrees <- function (fileRoot) c(lapply(tntDirectories, readTntTrees, nexusName=paste0(fileRoot, '.nex')),
                                     lapply(rDirectories, readRTrees, nexusName=paste0(fileRoot, '.nex')))
GetRTrees <- function (fileRoot) lapply(rDirectories, readRTrees, nexusName=paste0(fileRoot, '.nex'))

GetTreeScores <- function(fileRoot, trees = NULL) {
  fileRoot <- gsub('\\.csv$', '', fileRoot)
  scores <- NULL
  treeScoreFile <- paste0('islandCounts/', fileRoot, '.csv')
  if (file.exists(treeScoreFile)) {
    rawScores <- read.csv(treeScoreFile)
    scores <- as.matrix(rawScores[, -1])
    rownames(scores) <- rawScores[, 1]
  }
  if (is.null(nrow(scores)) || nrow(scores) != sum(nTrees)) {
    cat("   - Calculating tree scores...\n")
    if (is.null(trees)) trees <- GetTrees(fileRoot)
    nTrees <- vapply(trees, length, integer(1))
    flatTrees <- unlist(trees, recursive=FALSE)
    scores <- vapply(allDirectories, function (dirPath) {
      TreeScorer <- if (dirPath %in% tntDirectories) phangorn::fitch else inapplicable::InapplicableFitch
      rawData <- read.nexus.data(paste0(dirPath, '/', fileRoot, '.nex'))
      phyData <- phangorn::phyDat(rawData, type='USER', levels=c('-', 0:9))
      as.integer(vapply(flatTrees, TreeScorer, double(1), phyData, USE.NAMES=FALSE))
    }, integer(sum(nTrees)))
    rownames(scores) <- rep(c(tntDirectories, rDirectories), nTrees)
    write.csv(scores, file=treeScoreFile)
  }
  scores
}

GetVennTrees <- function (fileRoot, trees=GetTrees(fileRoot)) {
  nTrees <- vapply(trees, length, integer(1))
  treeDetails <- GetTreeScores(fileRoot, trees)
  
  trees <- trees[-1]
  treeDetails <- treeDetails[!(rownames(treeDetails) == 'ambiguous'), -1]
  extraSteps <- apply(treeDetails, 2, function (x) x - min(x))
 
  dimnames(extraSteps) <- NULL
  a_in_b_or_c <-  trees[[1]] %in% trees[[2]] | trees[[1]] %in% trees[[3]]
  b_in_c      <-  trees[[2]] %in% trees[[3]]
  if (any(c(a_in_b_or_c, b_in_c))) extraSteps <- extraSteps[-which(c(a_in_b_or_c, b_in_c)), ]
  
  treeIsOptimal <- extraSteps == 0
  colnames(treeIsOptimal) <- NULL
  
  vennTrees <- vapply(list(
    c(TRUE , FALSE, FALSE),
    c(FALSE, TRUE, FALSE),
    c(FALSE, FALSE, TRUE),
    c(TRUE , TRUE, FALSE),
    c(TRUE , FALSE, TRUE),
    c(FALSE , TRUE, TRUE),
    c(TRUE , TRUE , TRUE ),
    c(FALSE, FALSE, FALSE)), function (pattern) {
      sum(apply(treeIsOptimal, 1, identical, pattern))
    }, integer(1))
  
  if (sum(vennTrees) != nrow(extraSteps)) stop(fileRoot, ": Something's not right.")
  if (sum(vennTrees[-8]) != nrow(extraSteps)) warning(fileRoot, ": Suboptimal trees included in list.")
  vennTrees[-8]
}

Flatten <- function (trees) unlist(trees, recursive=FALSE)

GetRFDistances <- function (fileRoot, trees=GetTrees(fileRoot)) {
  flatTrees <- Flatten(trees)
  rfFileName <- paste0('treeSpaces/', fileRoot, '.rf.csv')
  if (file.exists(rfFileName)) {
    rfDistances <- data.matrix(read.csv(rfFileName, row.names=1))
  } else {
    cat(" - Calculating RF distances...\n")
    rfDistances <- RFDistances(flatTrees)
    write.csv(rfDistances, file=rfFileName)
  }
  rfDistances
}

GetQuartetDistances <- function (fileRoot, trees=GetTrees(fileRoot), forPlot=FALSE) {
  flatTrees <- Flatten(trees)
  qtFileName <- paste0('treeSpaces/', fileRoot, '.qt.csv')
  if (file.exists(qtFileName) && (file.mtime(qtFileName) > "2017-10-10 15:34:10 BST")) { # Don't read stale files, regenerate them
    qtDistances <- data.matrix(read.csv(qtFileName, row.names=1))
    if (length(trees) == 4) {
      if (sum(sapply(trees, length)) < ncol(qtDistances)) stop('Distances calculated from outdated trees.')
      ambigTrees <- seq_len(ncol(qtDistances) - sum(sapply(trees, length)))
      qtDistances <- qtDistances[-ambigTrees, -ambigTrees]
    }
  } else {
    cat(" - Calculating quartet distances for ", fileRoot, "...")
    if (length(trees) == 3) trees <- c(trees, GetRTrees(fileRoot))
    if (length(trees) != 4) trees <- GetTrees(fileRoot)
    qtDistances <- QuartetDistances(flatTrees)
    cat(" Done.\n")
    write.csv(qtDistances, file=qtFileName)
  }
  if (forPlot) {    
    qtDistances[qtDistances == 0] <- 1e-9
    diag(qtDistances) <- 0
  }
  qtDistances
}

modifiedPcoa <- function (D, correction = "none", rn = NULL) {
    centre <- function(D, n) {
        One <- matrix(1, n, n)
        mat <- diag(n) - One/n
        mat.cen <- mat %*% D %*% mat
    }
    bstick.def <- function(n, tot.var = 1, ...) {
        res <- rev(cumsum(tot.var/n:1)/n)
        names(res) <- paste("Stick", seq(len = n), sep = "")
        return(res)
    }
    D <- as.matrix(D)
    n <- nrow(D)
    epsilon <- sqrt(.Machine$double.eps)
    names <- rownames(D)
    
    CORRECTIONS <- c("none", "lingoes", "cailliez")
    correct <- pmatch(correction, CORRECTIONS)
    if (is.na(correct)) 
        stop("Invalid correction method")
    
    delta1 <- centre((-0.5 * D^2), n)
    trace <- sum(diag(delta1))
    D.eig <- eigen(delta1, TRUE) # We know we're symmetric
    min.eig <- min(D.eig$values)
    zero.eig <- which(abs(D.eig$values) < epsilon)
    D.eig$values[zero.eig] <- 0
    if (min.eig > -epsilon) {
        correct <- 1
        eig <- D.eig$values
        k <- length(which(eig > epsilon))
        rel.eig <- eig[1:k]/trace
        cum.eig <- cumsum(rel.eig)
        vectors <- sweep(D.eig$vectors[, 1:k], 2, sqrt(eig[1:k]), 
            FUN = "*")
        bs <- bstick.def(k)
        cum.bs <- cumsum(bs)
        res <- data.frame(eig[1:k], rel.eig, bs, cum.eig, cum.bs)
        colnames(res) <- c("Eigenvalues", "Relative_eig", "Broken_stick", 
            "Cumul_eig", "Cumul_br_stick")
        rownames(res) <- 1:nrow(res)
        rownames(vectors) <- names
        colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
            prefix = "Axis.")
        note <- paste("There were no negative eigenvalues. No correction was applied")
        out <- (list(correction = c(correction, correct), note = note, 
            values = res, vectors = vectors, trace = trace))
    } else {
        k <- n
        eig <- D.eig$values
        rel.eig <- eig/trace
        rel.eig.cor <- (eig - min.eig)/(trace - (n - 1) * min.eig)
        if (length(zero.eig) != 0) rel.eig.cor[zero.eig[1]] <- 0
        cum.eig.cor <- cumsum(rel.eig.cor)
        k2 <- length(which(eig > epsilon))
        k3 <- length(which(rel.eig.cor > epsilon))
        vectors <- sweep(D.eig$vectors[, 1:k2], 2, sqrt(eig[1:k2]), FUN = "*")
        if ((correct == 2) | (correct == 3)) {
            if (correct == 2) {
                c1 <- -min.eig
                note <- paste("Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 -", 
                  c1, ", except diagonal elements")
                D <- -0.5 * (D^2 + 2 * c1)
            }
            else if (correct == 3) {
                delta2 <- centre((-0.5 * D), n)
                upper <- cbind(matrix(0, n, n), 2 * delta1)
                lower <- cbind(-diag(n), -4 * delta2)
                sp.matrix <- rbind(upper, lower)
                c2 <- max(Re(eigen(sp.matrix, symmetric = FALSE, 
                  only.values = TRUE)$values))
                note <- paste("Cailliez correction applied to negative eigenvalues: D' = -0.5*(D +", 
                  c2, ")^2, except diagonal elements")
                D <- -0.5 * (D + c2)^2
            }
            diag(D) <- 0
            mat.cor <- centre(D, n)
            toto.cor <- eigen(mat.cor)
            trace.cor <- sum(diag(mat.cor))
            min.eig.cor <- min(toto.cor$values)
            zero.eig.cor <- which((toto.cor$values < epsilon) & 
                (toto.cor$values > -epsilon))
            toto.cor$values[zero.eig.cor] <- 0
            if (min.eig.cor > -epsilon) {
                eig.cor <- toto.cor$values
                rel.eig.cor <- eig.cor[1:k]/trace.cor
                cum.eig.cor <- cumsum(rel.eig.cor)
                k2 <- length(which(eig.cor > epsilon))
                vectors.cor <- sweep(toto.cor$vectors[, 1:k2], 
                  2, sqrt(eig.cor[1:k2]), FUN = "*")
                rownames(vectors.cor) <- names
                colnames(vectors.cor) <- colnames(vectors.cor, 
                  do.NULL = FALSE, prefix = "Axis.")
                bs <- bstick.def(k2)
                bs <- c(bs, rep(0, (k - k2)))
                cum.bs <- cumsum(bs)
            }
            else {
                if (correct == 2) 
                  cat("Problem! Negative eigenvalues are still present after Lingoes", 
                    "\\n")
                if (correct == 3) 
                  cat("Problem! Negative eigenvalues are still present after Cailliez", 
                    "\\n")
                rel.eig.cor <- cum.eig.cor <- bs <- cum.bs <- rep(NA, 
                  n)
                vectors.cor <- matrix(NA, n, 2)
                rownames(vectors.cor) <- names
                colnames(vectors.cor) <- colnames(vectors.cor, 
                  do.NULL = FALSE, prefix = "Axis.")
            }
            res <- data.frame(eig[1:k], eig.cor[1:k], rel.eig.cor, 
                bs, cum.eig.cor, cum.bs)
            colnames(res) <- c("Eigenvalues", "Corr_eig", "Rel_corr_eig", 
                "Broken_stick", "Cum_corr_eig", "Cum_br_stick")
            rownames(res) <- 1:nrow(res)
            rownames(vectors) <- names
            colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                prefix = "Axis.")
            out <- (list(correction = c(correction, correct), 
                note = note, values = res, vectors = vectors, 
                trace = trace, vectors.cor = vectors.cor, trace.cor = trace.cor))
        }
        else {
            note <- "No correction was applied to the negative eigenvalues"
            bs <- bstick.def(k3)
            bs <- c(bs, rep(0, (k - k3)))
            cum.bs <- cumsum(bs)
            res <- data.frame(eig[1:k], rel.eig, rel.eig.cor, 
                bs, cum.eig.cor, cum.bs)
            colnames(res) <- c("Eigenvalues", "Relative_eig", 
                "Rel_corr_eig", "Broken_stick", "Cum_corr_eig", 
                "Cumul_br_stick")
            rownames(res) <- 1:nrow(res)
            rownames(vectors) <- names
            colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                prefix = "Axis.")
            out <- (list(correction = c(correction, correct), 
                note = note, values = res, vectors = vectors, 
                trace = trace))
        }
    }
    class(out) <- "pcoa"
    out
}
