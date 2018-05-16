\name{calculateDEandDV}
\alias{calculateDEandDV}
\title{Conduct DE and DV analyses}
\description{
Calculate differentially expressed and differentially variable genes while statistically correcting for dropouts.
}
\usage{
calculateDEandDV(x, labels)
}
\arguments{
  \item{x}{Gene expression matrix of counts}
  \item{labels}{Vector of 0s and 1s for which group each sample belongs to}
}
\value{
Large Numeric of genes and corresponding Log Fold Change, t-test results and f-test results
}
\examples{
exprs <-as.matrix(read.table("matrix.txt", header=TRUE, sep="\t", row.names=1,as.is=TRUE))
DEDV<-calculateDEanDV(exprs, labels)
}
