tntDirectories <- c('ambiguous', 'ambigAbsent', 'extraState')
rDirectories <- c('inapplicable')
allDirectories <- c(tntDirectories, rDirectories)

nexusFiles <- list.files('matrices', pattern='.*\\.nex$');# nexusFiles
slowFiles <- c('Aria2015', 'Aguado2009', 'Capa2011', 'Conrad2008', 'Dikow2009', 'Eklund2004', 
               'Geisler2001', 'Giles2015', 'OMeara2014', 'Wills2012', 'Schulze2007')
slowNexus <- paste0(slowFiles, '.nex')
nexusFiles <- c(nexusFiles[!(nexusFiles %in% slowNexus)], slowNexus)
               

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

rfLegendPos <- list(
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
qtLegendPos <- list(
  'Asher2005' = 'topright',
  'Aria2015' = 'topright',
  'DeAssis2011' = 'bottomright',
  'Eklund2004' = 'bottomright',
  'Griswold1999' = 'topleft',
  'Geisler2001' = 'bottomleft',
  'Loconte1991' = 'topleft',
  'Longrich2010' = 'topleft',
  'OLeary1999' = 'bottomright',
  'Rougier2012' = 'bottomleft',
  'Rousset2004' = 'bottomleft',
  'Sano2011' = 'bottomright',
  'Sansom2010' = 'topright',
  'Shultz2007' = 'topleft',
  'Schulze2007' = 'bottomleft',
  'Vinther2008' = 'bottomleft',
  'Wilson2003' = 'bottomleft',
  'Wortley2006' = 'bottomleft'
)

studyName <- list(
  'Agnarsson2004' = 'Agnarsson 2004',           # AGNARSSON, I. 2004. Morphological phylogeny of cobweb spiders and their relatives (Araneae, Araneoidea, Theridiidae). Zoological Journal of the Linnean Society, 141, 447–626.
  'Capa2011'      = 'Capa et al. 2011',         # CAPA, M., HUTCHINGS, P., AGUADO, M. T. and BOTT, N. J. 2011. Phylogeny of Sabellidae (Annelida) and relationships with other taxa inferred from morphology and multiple genes. Cladistics, 27, 449–469.
  'DeAssis2011'   = "De Assis & Christoffersen\n2011",
                                                # DE ASSIS, J. E. and CHRISTOFFERSEN, M. L. 2011. Phylogenetic relationships within Maldanidae (Capitellida, Annelida), based on morphological characters. Systematics and Biodiversity, 9, 233–245.
  'OLeary1999'    = "O'Leary & Geisler 1999",   # O’LEARY, M. A. and GEISLER, J. H. 1999. The position of Cetacea within Mammalia: phylogenetic analysis of morphological data from extinct and extant taxa. Systematic Biology, 48, 455–490.
  'Rousset2004'   = 'Rousset et al. 2004',      # ROUSSET, V., ROUSE, G. W., SIDDALL, M. E., TILLIER, A. and PLEIJEL, F. 2004. The phylogenetic position of Siboglinidae (Annelida) inferred from 18S rRNA, 28S rRNA and morphological data. Cladistics, 20, 518–533.
  'Sano2011'      = 'Sano & Akimoto 2011',      # SANO, M. and AKIMOTO, S.-I. 2011. Morphological phylogeny of gall-forming aphids of the tribe Eriosomatini (Aphididae: Eriosomatinae). Systematic Entomology, 36, 607–627.
  'Sansom2010'    = 'Sansom et al. 2010',       # SANSOM, R. S., FREEDMAN, K., GABBOTT, S. E., ALDRIDGE, R. J. and PURNELL, M. A. 2010. Taphonomy and affinity of an enigmatic Silurian vertebrate, Jamoytius kerwoodi White. Palaeontology, 53, 1393–1409.
  'Schulze2007'   = 'Schulze et al. 2007',      # SCHULZE, A., CUTLER, E. B. and GIRIBET, G. 2007. Phylogeny of sipunculan worms: A combined analysis of four gene regions and morphology. Molecular Phylogenetics and Evolution, 42, 171–92.
  'Shultz2007'    = 'Shultz 2007',              # SHULTZ, J. W. 2007. A phylogenetic analysis of the arachnid orders based on morphological characters. Zoological Journal of the Linnean Society, 150, 221–265.
  'Wetterer2000'  = 'Wetterer et al. 2000',     # WETTERER, A. L., ROCKKMAN, M. V. and SIMMONS, N. B. 2000. Phylogeny of phyllostomid bats (Mammalia: Chiroptera): data from diverse morphological systems, sex chromosomes, and restriction sites. Bulletin of the American Museum of Natural History, 248, 1–200.
  'Wills2012'     = 'Wills et al. 2012',        # WILLS, M. A., GERBER, S., RUTA, M. and HUGHES, M. 2012. The disparity of priapulid, archaeopriapulid and palaeoscolecid worms in the light of new data. Journal of Evolutionary Biology, 25, 2056–2076.
  'Aguado2009'    = 'Aguado & San Martin 2009', # AGUADO, M. T. and SAN MARTÍN, G. 2009. Phylogeny of Syllidae (Polychaeta) based on morphological data. Zoologica Scripta, 38, 379–402.
  'Aria2015'      = 'Aria et al. 2015',         # ARIA, C., CARON, J. B. and GAINES, R. 2015. A large new leanchoiliid from the Burgess Shale and the influence of inapplicable states on stem arthropod phylogeny. Palaeontology, 58, 629–660.
  'Asher2005'     = 'Asher & Hofreiter 2006',   # ASHER, R. J. and HOFREITER, M. 2006. Tenrec phylogeny and the noninvasive extraction of nuclear DNA. Systematic biology, 55, 181–194.
  'Baker2009'     = 'Baker et al. 2009',        # BAKER, W. J., SAVOLAINEN, V., ASMUSSEN-LANGE, C. B., CHASE, M. W., DRANSFIELD, J., FOREST, F., HARLEY, M. M., UHL, N. W. and WILKINSON, M. 2009. Complete generic-level phylogenetic analyses of palms (Arecaceae) with comparisons of supertree and supermatrix approaches. Systematic Biology, 58, 240–256.
  'Bouchenak2010' = "Bouchenak-Khelladi\net al. 2010",
                                                # BOUCHENAK-KHELLADI, Y., VERBOOM, G. A., SAVOLAINEN, V. and HODKINSON, T. R. 2010. Biogeography of the grasses (Poaceae): a phylogenetic approach to reveal evolutionary history in geographical space and geological time. Botanical Journal of the Linnean Society, 162, 543–557.
  'Conrad2008'    = 'Conrad 2008',              # CONRAD, J. L. 2008. Phylogeny And Systematics Of Squamata (Reptilia) Based On Morphology. Bulletin of the American Museum of Natural History, 310, 1–182.
  'Dikow2009'     = 'Dikow 2009',               # DIKOW, T. 2009. A phylogenetic hypothesis for Asilidae based on a total evidence analysis of morphological and DNA sequence data (Insecta: Diptera: Brachycera: Asiloidea). Organisms Diversity and Evolution, 9, 165–188.
  'Eklund2004'    = 'Eklund et al. 2004',       # EKLUND, H., DOYLE, J. A. and HERENDEEN, P. S. 2004. Morphological phylogenetic analysis of living and fossil Chloranthaceae. International Journal of Plant Sciences, 165, 107–151.
  'Geisler2001'   = 'Geisler 2001',             # GEISLER, J. H. 2001. New morphological evidence for the phylogeny of Artiodactyla, Cetacea, and Mesonychidae. American Museum Novitates, 3344, 53.
  'Giles2015'     = 'Giles et al. 2015',        # GILES, S., FRIEDMAN, M. and BRAZEAU, M. D. 2015. Osteichthyan-like cranial conditions in an Early Devonian stem gnathostome. Nature, 520, 82–85.
  'Griswold1999'  = 'Griswold et al. 1999',     # GRISWOLD, C. E., CODDINGTON, J. A., PLATNICK, N. I. and FORSTER, R. R. 1999. Towards a phylogeny of entelegyne spiders (Araneae, Araneomorphae, Entelegynae). Journal of Arachnology, 27, 53–63.
  'Liljeblad2008' = 'Liljeblad et al. 2008',    # LILJEBLAD, J., RONQUIST, F., NIEVES-ALDREY, J. L., FONTAL-CAZALLA, F., ROS-FARRE, P., GAITROS, D. and PUJADE-VILLAR, J. 2008. A fully web-illustrated morphological phylogenetic study of relationships among oak gall wasps and their closest relatives (Hymenoptera: Cynipidae). .
  'Loconte1991'   = 'Loconte & Stevenson 1991', # LOCONTE, H. and STEVENSON, D. W. 1991. Cladistics of the Magnoliidae. Cladistics, 7, 267–296.
  'Longrich2010'  = 'Longrich et al. 2010',     # LONGRICH, N. R., SANKEY, J. and TANKE, D. 2010. ~Texacephale langstoni~, a new genus of pachycephalosaurid (Dinosauria: Ornithischia) from the upper Campanian Aguja Formation, southern Texas, USA. Cretaceous Research, 31, 274–284.
  'OMeara2014'    = "O'Meara &\nThompson 2014", # O’MEARA, R. N. and THOMPSON, R. S. 2014. Were There Miocene Meridiolestidans? Assessing the Phylogenetic Placement of Necrolestes patagonensis and the Presence of a 40 Million Year Meridiolestidan Ghost Lineage. Journal of Mammalian Evolution, 21, 271–284.
  'Rougier2012'   = 'Rougier et al. 2012',      # ROUGIER, G. W., WIBLE, J. R., BECK, R. M. D. and APESTEGUÍA, S. 2012. The Miocene mammal Necrolestes demonstrates the survival of a Mesozoic nontherian lineage into the late Cenozoic of South America. Proceedings of the National Academy of Sciences, 109, 20053–8.
  'Sharkey2011'   = 'Sharkey et al. 2011',      # SHARKEY, M. J., CARPENTER, J. M., VILHELMSEN, L., HERATY, J., LILJEBLAD, J., DOWLING, A. P. G., SCHULMEISTER, S., MURRAY, D., DEANS, A. R., RONQUIST, F., KROGMANN, L. and WHEELER, W. C. 2012. Phylogenetic relationships among superfamilies of Hymenoptera. Cladistics, 28, 80–112.
  'Sundue2010'    = 'Sundue et al. 2010',       # SUNDUE, M. A., ISLAM, M. B. and RANKER, T. A. 2010. Systematics of Grammitid Ferns (Polypodiaceae): Using Morphology and Plastid Sequence Data to Resolve the Circumscriptions of Melpomene and the Polyphyletic Genera Lellingeria and Terpsichore. Systematic Botany, 35, 701–715.
  'Vinther2008'   = 'Vinther et al. 2008',      # VINTHER, J., VAN ROY, P. and BRIGGS, D. E. G. 2008. Machaeridians are Palaeozoic armoured annelids. Nature, 451, 185–188.
  'Wilson2003'    = 'Wilson & Edgecombe 2003',  # WILSON, G. D. F. and EDGECOMBE, G. D. 2003. The Triassic isopod ~Protamphisopus wianamattensis~ (Chilton) and comparison by extant taxa (Crustacea, Phreatoicidea). Journal of Paleontology, 77, 454–470.
  'Wortley2006'   = 'Wortley & Scotland 2006',  # WORTLEY, A. H. and SCOTLAND, R. W. 2006. The effect of combining molecular and morphological data in published phylogenetic analyses. Systematic Biology, 55, 677–685.
  'Zanol2014'     = 'Zanol et al. 2014',        # ZANOL, J., HALANYCH, K. M. and FAUCHALD, K. 2014. Reconciling taxonomy and phylogeny in the bristleworm family Eunicidae (polychaete, Annelida). Zoologica Scripta, 43, 79–100.
  'Zhu2013'       = 'Zhu et al. 2013'           # ZHU, M., YU, X., AHLBERG, P. E., CHOO, B., LU, J., QIAO, T., QU, Q., ZHAO, W., JIA, L., BLOM, H. and ZHU, Y. 2013. A Silurian placoderm with osteichthyan-like marginal jaw bones. Nature, 502, 188–193.
)

