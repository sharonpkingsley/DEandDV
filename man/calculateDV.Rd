\name{calculateDV}
\alias{calculateDV}
\title{Calculate variance of each gene}
\description{
Calcuates the variance for each gene, statistically correcting for dropouts
}
\usage{
calculateDV(x, plot=FALSE)
}
\arguments{
  \item{x}{Expression matrix for which to calculate variance for each gene}
  \item{plot}{Display plot of variance without dropout correction vs corrected variance}
}
\value{
Large Numeric of genes and corresponding variance values
}
\examples{
exprs <-as.matrix(read.table("matrix.txt", header=TRUE, sep="\t", row.names=1,as.is=TRUE))
variances<-calculateDV(exprs)
}

