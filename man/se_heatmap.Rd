\name{se.heatmap}
\alias{se.heatmap}
\title{Construction of Multiple-Sample Side Effect Heatmap}
\usage{
se.heatmap(se, samples = NULL,
	hm.colors = c("lightgray", "darkorange"), ...)
}
\arguments{
\item{se}{List of computed side effects}
\item{samples}{Vector of sample names; Defaults to labeling them numerically.}
\item{hm.colors}{Vector of colors to be shown, first indicating 'not common' and latter 'common'.}
\item{...}{Additional arguments passed to gplots::heatmap.2}
}
\description{
A wrapper function for 'heatmap.2' in gplots to construct a similar figure as Fig 3C in Lin et al., PLoS ONE 2019.
}
\examples{
se.heatmap(se)
}
\keyword{side effects, predict, heatmap}