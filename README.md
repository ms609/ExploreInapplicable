# Exploring the Inapplicable Algorithm

## Obtaining matrices

Morphological matrices were selected from the literature with a view to obtaining a broad coverage
of both taxa and taxonomic levels.

Matrices with a trivial amount of inapplicable tokens (e.g. max 2 per transformation series)
were excluded, as were those with large numbers of taxa (>80), as an efficient implementation
of tree search algorithms that works with the inapplicable algorithm has not yet been constructed.

Comments and additional blocks were removed from the matrices to allow their easy parsing by
different software packages, and the matrices placed in the `matrices` directory.

Details of each matrix, including a full citation and interpretation details for each 
character, are provided in `r_dataset_details.Rmd`, [viewable here](https://cdn.rawgit.com/ms609/ExploreInapplicable/master/r_dataset_details.html).

## Annotating character types

Original character lists were obtained for each matrix.  For each matrix, a template annotation
file was generated using the Perl script `character_type_teplate.pl`.  This blank template
marked the characters that contained at least one inapplicable token.  Each of these characters
was then examined in the original character list and scored as neomorphic or transformational,
with the 'non-derived' character state flagged for neomorphic characters.

## Rewriting matrices

Each matrix was then copied to each of four subfolders, and modified for analysis:

- Folder `ambiguous`:
  *  all inapplicable tokens replaced with `?`

- Folder `ambigAbsent`: 
  * inapplicable tokens in transformational series characters with `?`
  * inapplicable tokens in neomorphic characters series replaced with non-derived token
      for that character
    
    
- Folder `extraState`: 
  * inapplicable tokens in transformational series characters with `9`
  * inapplicable tokens in neomorphic characters series replaced with non-derived token
      for that character
    
    
- Folder `inapplicable`: 
  * inapplicable tokens left as `-`
    
The Perl script `analyse_matrices.pl` performs these modifications.

## Identifying optimal trees

The Perl script `analyse_matrices.pl` also initiates parsimony analysis on most of the matrices.

Parsimony analysis is conducted in [TNT](http://www.lillo.org.ar/phylogeny/tnt/) using the 
script specified in `tnt_search.run` - in summary, using sectorial search,
the parsimony ratchet and tree drifting, to find the optimal tree length 100 times.  All
most parsimonious trees are saved to a file in the respective folder, and these files are 
converted from TNT's proprietary format to NEXUS format in a separate file ending `.nextrees`.
A strict consensus tree depicting the results of each analysis for are presented for each
matrix in PDF format in the `consTrees` folder.

Parsimony analysis using the new algorithm is conducted in R using the script `r_search.R`.

## Analysis of optimal trees

Once optimal trees have been collected under each method, the 'islands' of optimal trees were
compared using two methods, implemented in the script `r_plot.R`.

### Method 1. Length of trees on each island.

The MPTs on each island were analysed by the other methods.  The scores that were optimal under
one method are often suboptimal under another.  The number of extra steps associated with each
tree -- in a sense, how far the trees that are optimal under one method are from those
that are optimal under another -- is counted in `islandCounts` and plotted as a histogram
in the `histograms` directory.

### Method 2. Overlap of islands in tree space.

Distances between each pair of trees were calculated using the Robinson-Foulds distance and the
Quartet distances.  Principle components were generated from these distance matrices using the 
`ape` function `pcoa`, and this crude tree-space was plotted in two dimensions, with convex
hulls drawn around the trees derived from each method.  These treespaces are saved in the
`treeSpaces` directory.





