\name{prep.exparrays}
\alias{prep.exparrays}
\title{Analysis of Agilent 1-color Expression Arrays Data with LIMMA}
\usage{
prep.exparrays(targets, 
	annotations = NULL, use.ensembl = TRUE,
	plots = TRUE, norm.method = NULL, offset = NULL)
}
\arguments{
\item{targets}{Path to LIMMA-style targetfile. See 'limma::readTargets' for more details.}
\item{annotations}{Path to genome annotation file. Default option (NULL) is to use built-in annotations in the package (hg19).
}
\item{use.ensembl}{Should Ensembl annotations be used as well? Defaults to TRUE.}
\item{plots}{Should diagnostic MA-plots be created? Defaults to TRUE.}
\item{norm.method}{Array normalization method: choose between 'scale' or 'quantile' normalization. 'NULL' will try both and decide on a suitable method based on how variances are stabilized.}
\item{offset}{Amount of background offset to reduce variability in low intensity spots. Setting this to 'NULL' allows the function to determine one automatically.}
}
\description{
This is a semi-automated processing function to calculate differential gene expressions from 1-color expression arrays by limma. For the 'annotation' option, the input annotation file should be a 6-column BED-like file in hg19 coordinates, with the following columns: [chromosome start end accession 1/2 +/- symbol]. By setting 'use.noensembl' to FALSE, the function will only annotate data that have refFlat annotations (this is the same operating mode as Lin et. al, PLoS ONE 2016). 

For 'norm.method', the script will try to determine the optimal method based on the maximization of prior variances; it is recommended to leave this and 'offset' to the default settings for your first run, then set them in subsequent cross-sample or data re-analysis.
}
\examples{
# processes expression array data with scale normalization and an offset of 16.
prep.motif(targets = "limma_sample.txt", plots = FALSE, norm.method = "scale", offset = 16)
}
\keyword{polyamide, off-targets}