# clust: efficient spatial and spatio-temporal cluster detection using the LASSO
_cluST_


*clust* is an R package which can be used to efficiently detect spatial clusters and spatio-temporal clusters with count data on a lattice. Examples include rates of disease in a population broken down by a geographic unit. 

This package is still under development. Please report bugs or constructive
tips to issues [here](https://github.com/mkamenet3/cluST/issues)

Detailed examples of use can be found in the [vignette](https://github.com/mkamenet3/cluST/tree/renamefuncs/scripts/clust/vignettes)


The main functions in this package are ```clust``` and ```clust_sim``` which
both perform both spatial and spatio-temporal cluster detection using the lasso. The
former works with observed and expected data while the former is used for the
simulation of results. Diagnostic tools are also provided for further analysis
of background data. 

Additionally, we provide the flexibility for the user to further explore
risk ratios from the output as well as provide colormapping. These can then easily be mapped onto a user-generated map.


## Installation

*clust* was built on R version 3.1.2. 

Other package dependencies are:

- geosphere
- glmnet (>= 2.0.0)
- maps
- MASS  (>= 7.3.0)
- Matrix (>= 1.2-4)
- RColorBrewer
- Rcpp (>= 0.12.4)
- sp
- spdep


To download the latest version of *clust*:

```R
library("devtools")
devtools::install_github("mkamenet3/cluST")
```



