times <- as.data.frame(read.csv('searchTimes.csv', header=TRUE))
colnames(times) <- c('hits', 'time', 'dataset', 'length')

thisRow <- seq_len(nrow(times))[-1]
lastRow <- seq_len(nrow(times) - 1)
sameDataset <- times[thisRow, 'dataset'] == times[lastRow, 'dataset']
lowerScore <- times[thisRow, 'length'] < times[lastRow, 'length']
times[, 'newScore'] <- FALSE
times[thisRow, 'newScore'] <- sameDataset & lowerScore

hitStats <- vapply(sort(unique(times[, 'hits'])), function (nHits) {
  rows <- times[, 'hits'] == nHits
  nReps <- sum(rows)
  timeSpent <- sum(times[rows, 'time'])
  repsLucky <- sum(times[rows, 'newScore'])
  timeLucky <- sum(times[rows & times[, 'newScore'], 'time'])
  return (c(nHits, nReps, timeSpent, repsLucky, timeLucky))
}, double(5))
rownames(hitStats) <- c('nHits', 'nReps', 'timeSpent', 'repsLucky', 'timeLucky')
hitStats <- as.data.frame(t(hitStats))

plot(I(repsLucky / nReps) ~ nHits, data=hitStats)
plot(I(timeLucky / timeSpent) ~ nHits, data=hitStats)
plot(I(repsLucky / timeSpent) ~ nHits, data=hitStats)

hitGroups <- list(1:8, 9:15, 16:18, 19:21) # Warning: MANUALLY DEFINED
hitStats2 <- as.data.frame(t(vapply(hitGroups, function (x) colSums(hitStats[x, ]), double(5))))
paste0(c('2-9: ', '10-16: ', '25-60: ', '100-200: '), round(hitStats2[, 'repsLucky'] / hitStats2[, 'nReps'] * 100, 4), '%')
paste0(c('2-9: ', '10-16: ', '25-60: ', '100-200: '), round(hitStats2[, 'timeLucky'] / hitStats2[, 'timeSpent'] * 100, 4), '%')


