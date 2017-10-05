# Exploring the Inapplicable Algorithm

## Obtaining matrices

Morphological matrices were selected from the literature with a view to obtaining a broad coverage
of both taxa and taxonomic levels.

Matrices with a trivial amount of inapplicable tokens (e.g. max 2 per transformation series)
were excluded, as were those with large numbers of taxa (>80), as an efficient implementation
of tree search algorithms that works with the inapplicable algorithm has not yet been constructed.

Comments and additional blocks were removed from the matrices to allow their easy parsing by
different software packages, and the matrices placed in the `matrices` directory.

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

## Analysis