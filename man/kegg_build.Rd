\name{prep.kegg}
\alias{prep.kegg}
\title{Generation of a Local Human KEGG Genome Database}
\usage{
prep.kegg(prefix = "pipoft", threads = 1)
}
\arguments{
\item{prefix}{How should the files be named? Defaults with prefix 'pipoft'.}
\item{threads}{How many CPU threads should be used? Defaults to 1.}
}
\description{
'prep.kegg' accesses the KEGG database via its open API to create a local version of the human data for 'pipoft'. At this time only the human genome (identifier: hsa) is supported. Multithreading is supported when 'threads' is set to a valid value greater than 1.
}
\examples{
prep.kegg() # not much else is needed, really
}
\keyword{kegg, off-targets}