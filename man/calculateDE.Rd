\name{calculateDE}
\alias{calculateDE}
\title{Calculate expected value of each gene}
\description{
Calcuates the expected expression value for each gene, statistically correcting for dropouts
}
\usage{
calculateDE(x, plot=FALSE)
}
\arguments{
  \item{x}{Expression matrix for which to calculate expected expression for each gene}
  \item{plot}{Display plot of means without dropout correction vs corrected means}
}
\value{
Large Numeric of genes and corresponding expected expression values
}
\examples{
exprs <-as.matrix(read.table("matrix.txt", header=TRUE, sep="\t", row.names=1,as.is=TRUE))
expectedExpression<-calculateDE(exprs)
}

