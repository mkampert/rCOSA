#' A Typical Synthetic Data Set from a COSA Model
#'
#' A typical synthetic data set on which rCOSA would do very good. To generate this data set yourself see the example.
#'  
#'   A Monte Carlo data set for 100 objects with 1,000 attributes which contains two groups, and background noise. The two groups share a subset of 15 attributes, of which each single observation of an object on its attributes comes from a normal distribution with mean -1.5 and standard deviation 0.2 Each group also has its own uniques subets of 15 attributes. Each observation of an object and attribute in these non-overlapping subsets are from a normal distribution with mean 1.5 and standard deviation 0.2. The background noise is i.i.d. and comes from a standard normal distribution.
#' 
#' 
#'
#' @docType data
#' @keywords datasets
#' @name Xtoy
#' @usage data(Xtoy)
#' @format A data frame with 100 rows and 1000 variables
NULL
