# Quick explanation
These plots show the MDS tree spaces, created using the Quartet distances between each 
pair of trees.  Tree spaces based on the Robinson Foulds distance are depicted in the folder 
`pcoaTreeSpaces`.

The different colours correspond to different optimizations:

- Ambiguous: all inapplicable tokens replaced with ?s and analysed in TNT
- AmbigAbsent: inapplicable tokens replaced with 'absent' in neomorphic characters,
               'ambiguous' in transformational characters; analysis in TNT
- extraState: inapplicable tokens replaced with 'absent' in neomorphic characters, 
              '9' (i.e. extra state) in transformational characters; analysis in TNT
- inapplicable: inapplicable tokens replaced with 'absent' in neomorphic characters, 
              '-' (i.e. inapplicable state) in transformational characters;
              analysis in R using inapplicable, which calls MorphyLib.
              
Searches were run until the best score had been independently hit 100 times (which means something
slightly different in the different implementations) and the number of unique MPTs is reported.