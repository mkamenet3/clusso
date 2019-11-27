# clusso: efficient spatial and spatio-temporal cluster detection using the LASSO
_clusso_




*clusso* is an R package is package that efficiently detects across space and space-time using the LASSO (least absolute shrinkage and selection operator). Current models include the Poisson and Binomial, as well as their quasi- counter parts, allowing for overdispersion. 

This package uses information critera ((Q)BIC, (Q)AIC, (Q)AICc) as well as cross-validation to determine the optimal path and thereby number of clusters in the study region. 

This package is still under development. Please report bugs or constructive tips to issues [here](https://github.com/mkamenet3/clusso/issues)

Detailed examples can be found in the [vignette](github.com/mkamenet3/clusso/tree/master/vignettes)


The main function in this package is ```clusso()```, which performs both spatial and spatio-temporal cluster detection using the LASSO. Cluster detection can be performed across both space and space and time. Diagnostic tools are also provided for further analysis as well as plotting functions. 



## Installation

*clusso* was built on R version 3.6.1. 

Package dependencies include:

- tidyverse (>= 1.3.0)
- glmnet (>= 3.0.1)
- Matrix (>= 1.2-17)
- data.table (>= 1.12.6)

Package imports include:

- geosphere (>= 1.5-10)
- MASS  (>= 7.3-51.4)
- RColorBrewer (>= 1.1-2)
- Rcpp (>= 1.0.3)
- SOAR (>= 0.99-11)
- stringi (>= 1.4.3)





To download the latest version of *clusso*:

```R
library("devtools")
devtools::install_github("mkamenet3/clusso")
```



