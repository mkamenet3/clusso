# clusso: efficient spatial and spatio-temporal cluster detection using the LASSO
_clusso_




*clusso* is an R package is package efficiently detects cluster in sparse space
and space-time matrices. Current options are to use the Poisson and quasi-Poisson models, the latter allowing for overdispersion. It uses AIC/AICc/BIC and their quasi-likelihood equivalents to determine the optimal path.

These functions parse out clusters from a background rate, smoothing over the background rate via the regularization. This package also provides extensions for simulations and related simulation diagnostics as well as color-coding for easy-to-make maps. 

This package is still under development. Please report bugs or constructive
tips to issues [here](https://github.com/mkamenet3/clusso/issues)

Detailed examples of use can be found in the [vignette](github.com/mkamenet3/clusso/tree/master/vignettes)


The main function in this package is ```clusso()```, which performs both spatial and spatio-temporal cluster detection using the LASSO. Cluster detection can be performed across both space and space and time. Diagnostic tools are also provided for further analysis as well as plotting functions. 



## Installation

*clusso* was built on R version 3.5.2. 

Package dependencies include:

- tidyverse (>= 1.2.1)
- glmnet (>= 2.0.0)
- Matrix (>= 1.2-4)
- data.table (>= 1.12.0)

Package imports include:

- geosphere
- MASS  (>= 7.3.0)
- RColorBrewer
- Rcpp (>= 1.0.0)
- SOAR (>= 0.99-11)
- stringi (>= 1.2.4)





To download the latest version of *clusso*:

```R
library("devtools")
devtools::install_github("mkamenet3/clusso")
```



