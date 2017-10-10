tntDirectories <- c('ambiguous', 'ambigAbsent', 'extraState')
rDirectories <- c('inapplicable')
allDirectories <- c(tntDirectories, rDirectories)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
treePalette <- cbPalette[c(6, 3, 2, 8)]
names(treePalette) <- allDirectories
plotChars <- c(2, 1, 3, 4)
names(plotChars) <- allDirectories
englishName <- c('Ambiguous (raw)', 'Ambiguous', 'Extra state', 'Inapplicable')
names(englishName) <- allDirectories

QuartetLegendPos <- function (fileRoot) {
  hardCoded <- qtLegendPos[[fileRoot]]
  return (if (is.null(hardCoded)) 'bottomright' else hardCoded)
}

rfLegendPos = list(
  'Asher2005' = 'topright',
  'Aria2015' = 'topright',
  'DeAssis2011' = 'bottomright',
  'Eklund2004' = 'bottomright',
  'Griswold1999' = 'bottomright',
  'Laconte1991' = 'topleft',
  'OLeary1999' = 'bottomright',
  'Rousset2004' = 'bottomleft',
  'Sano2011' = 'bottomleft',
  'Sansom2010' = 'bottomright',
  'Wortley2006' = 'bottomleft',
  'Geisler2001' = 'bottomleft'
)
qtLegendPos = list(
  'Asher2005' = 'topright',
  'Aria2015' = 'topright',
  'DeAssis2011' = 'bottomright',
  'Eklund2004' = 'bottomright',
  'Griswold1999' = 'topleft',
  'Laconte1991' = 'topleft',
  'OLeary1999' = 'bottomright',
  'Rousset2004' = 'bottomleft',
  'Sano2011' = 'bottomleft',
  'Sansom2010' = 'bottomright',
  'Wortley2006' = 'bottomleft',
  'Geisler2001' = 'bottomleft'
)


