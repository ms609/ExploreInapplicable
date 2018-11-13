# Quick explanation
These plots show the MDS tree spaces of all optimal trees.
The top row defines tree distances using the Robinson Foulds distances between each pair of trees;
the bottom row uses the quartet distance.

The left column includes the trees that are optimal when the raw data matrices (`Ambiguous`) 
are used; the right column excludes, these trees, only comparing trees where inapplicable
tokens in neomorphic characters have been replaced with the 'absent' state.

The matrices in the right-hand column, therefore, only differ in the treatment of inapplicable
tokens in transformational series.

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