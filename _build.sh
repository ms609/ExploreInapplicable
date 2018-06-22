#!/bin/sh

Rscript -e "devtools::install_github('ropensci/rcrossref')" # Need version 0.8.1 or higher
Rscript -e "devtools::install_github('ms609/TreeSearch')" # Need version > 0.1.2
Rscript -e "rmarkdown::render('r_dataset_details.Rmd', 'bookdown::gitbook')"
