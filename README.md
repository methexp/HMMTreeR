[![R build status](https://github.com/methexp/HMMTreeR/workflows/R-CMD-check/badge.svg)](https://github.com/methexp/HMMTreeR/actions)

# HMMTreeR

HMMTreeR is an R interface to the program HMMTree. It allows fitting latent-class multinomial processing tree models.

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)


Stahl, C., & Klauer, K.C. (2007).
HMMTree: A computer program for hierarchical multinomial processing tree models. *Behavior Research Methods*, *39*, 267-273.
doi: [10.3758/BF03193157](https://doi.org/10.3758/BF03193157)

Klauer, K.C. (2006).
Hierarchical multinomial processing tree models: A latent-class approach. *Psychometrika*, *71*, 7-31.
doi: [10.1007/s11336-004-1188-3](https://doi.org/10.1007/s11336-004-1188-3)

See the package vignette for a short introduction.
The model fitting can only be executed on a Windows machine.

To install `HMMTreeR`, use:

```
if(!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("methexp/HMMTreeR")
```
