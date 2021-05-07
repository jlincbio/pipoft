\name{kegg.pathway}
\alias{kegg.pathway}
\title{Generates Color Coding to Display Gene Expressions in a KEGG Pathway}
\usage{
kegg.pathway(x, pathway, ...)
}
\arguments{
\item{x}{Expression data table created by 'prep.exparrays'.}
\item{pathway}{Human KEGG pathway accession ID (hsa...).}
\item{...}{Arguments to be passed to 'kegg.fcmap'.}
}
\description{
'kegg.pathway' retrieves a list of gene identifiers from KEGG database for a specified 'pathway', then creates relevant resource files to be used with other KEGG functions such as 'kegg.colors' to layer color-coded heatmap-like expression data onto KEGG pathways. Color arguments such as 'color.pos' and 'color.neg' can be passed to 'kegg.fcmap' should color other than the default 'red-white-green' tricolor gradient be desired. Required fields in 'x' are 'logFC' and 'Symbol' (standard gene symbols). 
}
\examples{
kegg.pathway(x = "fcTable.txt", pathway = "hsa05212") # maps pancreatic cancer pathway data with expressions in fcTable.txt
}
\keyword{kegg, off-targets}