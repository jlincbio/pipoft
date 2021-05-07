\name{se.predict}
\alias{se.predict}
\title{Construction of Multiple-Sample Side Effect Heatmap}
\usage{
se.predict(training, effects, expr, output = NULL, ...)
}
\arguments{
\item{training}{List of computed side effects}
\item{effects}{Vector of sample names; Defaults to labeling them numerically.}
\item{expr}{Vector of colors to be shown, first indicating 'not common' and latter 'common'.}
\item{output}{}
\item{...}{Additional arguments passed to gplots::heatmap.2}
}
\description{
A wrapper function for 'heatmap.2' in gplots to construct a similar figure as Fig 3C in Lin et al., PLoS ONE 2019.
}
\keyword{side effects, predict, heatmap}

