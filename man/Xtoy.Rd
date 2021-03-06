% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Xtoy.r
\docType{data}
\name{Xtoy}
\alias{Xtoy}
\title{A Typical Synthetic Data Set from a COSA Model}
\format{A data frame with 100 rows and 1000 variables}
\usage{
data(Xtoy)
}
\description{
A typical synthetic data set on which rCOSA would do very good. To generate this data set yourself see the example.
}
\details{
A Monte Carlo data set for 100 objects with 1,000 attributes which contains two groups, and background noise. The two groups share a subset of 15 attributes, of which each single observation of an object on its attributes comes from a normal distribution with mean -1.5 and standard deviation 0.2 Each group also has its own uniques subets of 15 attributes. Each observation of an object and attribute in these non-overlapping subsets are from a normal distribution with mean 1.5 and standard deviation 0.2. The background noise is i.i.d. and comes from a standard normal distribution.
}
\keyword{datasets}
