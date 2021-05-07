\name{kegg.colors}
\alias{kegg.colors}
\title{Creates a KEGG Pathway Colored with Expression Data}
\usage{
kegg.colors(pathway, colors)
}
\arguments{
\item{pathway}{KEGG pathway accession ID (e.g. hsa05212 for pancreatic cancer).}
\item{colors}{Name of expression color maps corresponding to KEGG pathway specified in 'pathway'; use file prefixed as 'colorMap_' from 'kegg.pathway'.}
}
\description{
'kegg.colors' tries to create a color-coded KEGG pathway image via KEGG's open API. Results from 'kegg.pathway' can be directly routed to 'kegg.colors' by specifying the name of the color code output. If a graphical device can be opened, a preview of the generated map will be displayed. If the function has trouble saving the file, an API-generated URL will be displayed.
}
\examples{
kegg.colors(pathway = "hsa05212", colors = "colorMap_hsa05212.txt")
}
\keyword{kegg, off-targets}