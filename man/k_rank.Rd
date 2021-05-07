\name{kegg.rank}
\alias{kegg.rank}
\title{Estimation of Off-Target Genes and KEGG Pathways by PI Polyamide}
\usage{
kegg.rank(expr, targets, kegg = NULL, sites = NULL,
	prefix = "keggRank", path = getwd(), alpha = 0.05)
}
\arguments{
\item{expr}{Expression table generated from 'prep.exparrays'.}
\item{targets}{Vector of genes assumed to be on-target.}
\item{kegg}{Generated KEGG gene database (from 'prep.kegg').}
\item{sites}{BED file of Chem-seq binding data, if available.}
\item{alpha}{Significance level (defaults to 0.05).}
\item{path}{Where should output be saved? Defaults to current directory.}
\item{prefix}{How should the files be named? Default with prefix 'pipoft'.}
}
\description{
This function performs the proposed off-target estimation procedure based on differential expressions and Chem-seq binding information (if available). The call returns a table considering name of the gene, the corresponding class, relative rank and relevant statistics and significance. The output is also saved in a tab-delimited table file.
}
\examples{
p.table <- kegg.rank(expr = "arrayResults.txt", targets = "KRAS")
# sets KRAS as primary target
}
\keyword{kegg, off-targets}