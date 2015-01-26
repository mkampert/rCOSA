\name{dendro}
\alias{dendro}
\title{Plot a dendrogram}
\usage{
  dendro(Dissim = CO1[[1]], main = "",
    leaflab = c("perpendicular", "textlike", "none")[3],
    classnrs = y, cex.main = NULL, denplot = T, sub = "",
    xlab = "", ylab = "", condpop)
}
\value{
  Returns a plot of a dendrogram of an object of class
  dist. When the groups are known and filled out, then each
  each object has its own group color.
}
\description{
  Plot a dendrogram
}

