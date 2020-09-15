# clusso: Regularized Spatial and Spatio-Temporal Cluster Detection
_clusso_




*clusso* is an R package based on "Regularized Spatial and Spatio-Temporal Cluster Detection" (Kamenetsky, Lee, Zhu, Gangnon). *clusso* implements Poisson and quasi-Poisson regression with the Lasso (least absolute shrinkage and selection operator) penalty to identify a parsimonious set of disease clusters with elevated or depressed risk relative to the background rate in a study region. The number of clusters and tuning parameters are selected based on (quasi-)information criteria.

Code developed by and repository maintained by M.Kamenetsky.

This package is still under development. Please report bugs or constructive tips to issues [here](https://github.com/mkamenet3/clusso/issues)

Detailed examples can be found in the [vignette](github.com/mkamenet3/clusso/tree/master/vignettes)


The main function in this package is ```clusso()```, which performs both spatial and spatio-temporal cluster detection using the Lasso. Cluster detection can be performed across both space and space and time. Diagnostic tools are also provided for further analysis as well as plotting functions. 



## Installation

*clusso* was built on `R` version 4.0.2 "Taking Off Again" 

Package dependencies include:


- glmnet (>= 3.0.1)
- Matrix (>= 1.2-17)


Package imports include:

- ggplot2 (>= 3.2.1),
- dplyr (>= 0.8.3),
- tidyr (>= 1.0.0),
- magrittr (>= 1.5),
- data.table (>= 1.12.6),
- geosphere (>= 1.5-10),
- MASS (>= 7.3-51.4),
- RColorBrewer (>= 1.1-2),
- Rcpp (>= 1.0.3),
- SOAR (>= 0.99-11),
- stringi (>= 1.4.3)






To download the latest version of *clusso*:

```R
library("devtools")
devtools::install_github("mkamenet3/clusso")
```



