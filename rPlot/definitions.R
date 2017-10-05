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


rfLegendPos = list(
  'Asher2005.nex' = 'topright',
  'Aria2015' = 'topright',
  'DeAssis2011.nex' = 'bottomright',
  'Eklund2004.nex' = 'bottomright',
  'Griswold1999.nex' = 'bottomright',
  'Laconte1991.nex' = 'topleft',
  'OLeary1999.nex' = 'bottomright',
  'Rousset2004' = 'bottomleft',
  'Sano2011.nex' = 'bottomleft',
  'Sansom2010.nex' = 'bottomright',
  'Wortley2006.nex' = 'bottomleft',
  'Geisler2001.nex' = 'bottomleft'
)
qtLegendPos = list(
  'Asher2005.nex' = 'topright',
  'Aria2015' = 'topright',
  'DeAssis2011.nex' = 'bottomright',
  'Eklund2004.nex' = 'bottomright',
  'Griswold1999.nex' = 'bottomright',
  'Laconte1991.nex' = 'topleft',
  'OLeary1999.nex' = 'bottomright',
  'Rousset2004' = 'bottomleft',
  'Sano2011.nex' = 'bottomleft',
  'Sansom2010.nex' = 'bottomright',
  'Wortley2006.nex' = 'bottomleft',
  'Geisler2001.nex' = 'bottomleft'
)


