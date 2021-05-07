\name{prep.drugdata}
\alias{prep.drugdata}
\title{Training Set Generation for Drug-Side Effect Prediction}
\usage{
prep.drugdata(drugs, platforms, effects)
}
\arguments{
\item{drugs}{Summary table of reference drugs to be used (tab-delimited).}
\item{platforms}{Summary table of GEO microarray platforms (tab-delimited).}
\item{effects}{Drug side-effect data (RData), if available.}
}
\description{
To generate a properly formatted training set, two summary tables are required. The 'drugs' summary table should contain a field 'Compound' indicating the name of the drug, 'Platform' linking the type of expression array dataset to its proper GEO platform name, 'LogTransformed' indicating whether log transformation of the signal data is necessary (TRUE or FALSE), and a 'File' field containing the proper path to the preprocessed expression data. The 'platforms' summary table should contain a minimum of two fields, one 'Platform' indicating the name of the GEO platform type and 'File' indicating the path to the definition file. Both files will need to be tab-delimited flat files. 'prep.drugdata' then performs the proper ordering and expression transformation to generate a base training dataset. By specifying 'effects', an analysis-ready matrix is created for random forest learning.
}
\examples{
mat <- build.effects(
	sider.id = "sider_drugNames.tsv",
	sider.effects = "sider_meddra_all_se.tsv")
}
\keyword{drugbank, sider}