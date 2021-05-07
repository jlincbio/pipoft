\name{kegg.fcmap}
\alias{kegg.fcmap}
\title{Generates Color Coding for Attributes to be Used Alongside KEGG Pathway Data}
\usage{
kegg.fcmap(fc, 
	col.neg = "red", col.pos = "green", col.zero = "white",
	col.na = "purple", append.border = TRUE, grain = 0.5)
}
\arguments{
\item{fc}{A numerical vector of attributes (e.g. logFC) to be converted.}
\item{col.neg}{Color to denote negative values (fc < 0). Defaults to red.}
\item{col.pos}{Color to denote positive values (fc > 0). Defaults to green.}
\item{col.zero}{Color to denote values around 0 (fc ~ 0). Defaults to white.}
\item{col.na}{Color to denote missing data. Defaults to purple.}
\item{append.border}{Should borders be included in color codes? Defaults to TRUE.}
\item{grain}{How fine should the color gradients be defined? Defaults to 0.5 (color shift every half unit).}
}
\description{
'kegg.fcmap' is a convenience function that bins input data to heatmap-like color gradients to bed used for color-coding KEGG pathways. If TRUE, 'append.border' will set missing data to have white borders. The current version is limited to three colors with breaks set to be every half unit (0.5) of the numerical range. As 'kegg.pathway' also calls this directly, in most situations a direct call of this function is not necessary, unless in cases when color changes may be desirable.
}
\examples{
# this creates a plot of random points with color gradient indicating sign
test1 <- rnorm(100, mean = 0, sd = 3)
test1.col <- kegg.fcmap(test1, grain = 0.25, append.border = FALSE)
plot(test1, col = test1.col, pch = 19)
}
\keyword{kegg, off-targets}