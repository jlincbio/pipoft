\name{se.corrplot}
\alias{se.corrplot}
\title{Correlation Plot for logFC in < 5 Samples}
\usage{
se.corrplot(tables, sample = NULL, fill.both = FALSE)
}
\arguments{
\item{tables}{List of expression data; needs to be specified as table = list(...). Each table should contain 'logFC' and 'Symbol' as column names (as generated by prep.exparrays).}
\item{sample}{Vector of sample names; needs to be the same length as list of tables. Defaults to automatically generate from the original numerical order.}
\item{fill.both}{Should the correlation between i-th and j-th sets be plotted on both sides of the diagonal? Defaults to FALSE (only upper diagonal).}
}
\description{
This function plots logFC values from a list of fold-change tables, and calculates the corresponding Sperman's correlation coefficient.}
\examples{
mat <- se.corrplot(tables = list(fc1, fc2, fc3) # plots fc1 vs. fc2, fc1 vs. fc2 and fc2 vs. fc3.
}
\keyword{correlations}