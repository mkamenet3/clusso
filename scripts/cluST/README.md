# cluST
_cluST_

*Space and Space-Time Cluster Detection Using the Lasso*

The main functions in this package are ```clust``` and ```clust.sim``` which
both perform both space and spacetime cluster detection using the lasso. The
former works with observed and expected data while the former is used for the
simulation of results. Diagnostic tools are also provided for further analysis
of background data. For example, percentage of simulations which accurately
detect cluster (percent true positives), percentage of simulations which
inaccurately put cluster in background rate (false negative), and percentage of
simulations which inaccurately determine a cluster (false positives).

Additionally, we provide the flexibility for the user to further explore
risk ratios from the output as well as provide colormapping using Rcolorbrewer.
These can then easily be mapped onto a user-generated map.

This package is still under development. Please report bugs or constructive
tips to issues [here] (https://github.com/mkamenet3/cluST/issues)
