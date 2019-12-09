
# rCOSA

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

rCOSA is an R package. The main output is a cluster happy dissimilarity matrix that can serve as input for proximity analysis methods.

## Installation

These are the commands to install and load `rCOSA`:

```r
install.packages('devtools');
devtools::install_github('mkampert/rCOSA')
```

## Example

A detailed overview on how to use the `rCOSA` package is given in the open-source article ["rCOSA: A Software Package for Clustering Objects on Subsets of Attributes"](https://link.springer.com/article/10.1007/s00357-017-9240-z) in the Journal of Classification (2017, Vol 34, issue 3, pp. 514 - 547).  

A quick basic example of code on how to use the `rCOSA` package is: 

```r
library(rCOSA)
data(ApoE3) # ?ApoE3
cosa_rslts <- cosa2(ApoE3)
hierclust(cosa_rslts$D) 
```

## Latest Developments

For the latest developments, check out ["Improved Strategies for Distance Based Clustering of Objects on Subsets of Attributes in High-Dimensional Data"](https://openaccess.leidenuniv.nl/handle/1887/74690)
