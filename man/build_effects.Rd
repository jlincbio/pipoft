\name{build.effects}
\alias{build.effects}
\title{Generates a Drug-Side Effect Data Matrix}
\usage{
build.effects(sider.id, sider.effects, kegg.filter = NULL,
	drugbank = NULL, drugbank.id = NULL, 
	drugbank.filter = (if (is.null(drugbank)) NULL else "inhibit"),
	cutoff = -1, path = getwd(), prefix = "pipoft")
}
\arguments{
\item{sider.id}{SIDER drug ID table (tab-delimited file).}
\item{sider.effects}{SIDER MedDRA side effects (table-delimited file).}
\item{kegg.filter}{Target gene list to limit search scope, if available.}
\item{drugbank}{DrugBank XML file for imposing keyword search restrictions.}
\item{drugbank.id}{List of DrugBank Uniprot ID's to limit search scope; required if kegg.filter is specified.}
\item{drugbank.filter}{What to search in DrugBank records? Defaults to inhibitory actions (rather, entries beginning with 'inhibit').}
\item{cutoff}{Minimum number of drugs for which a side effect is counted. Default: do not set a minimum number.}
\item{path}{Where should output be saved? Defaults to current directory.}
\item{prefix}{How should the files be named? Default with prefix 'pipoft'.}
}
\description{
This function generates a drug-side effect matrix to be used for model training. As neither DrugBank nor SIDER datasets are completely open-access, the user is expected to acquire use permissions and have the database stored locally. A 'm x n' matrix is returned, with 'm' drugs in the row and 'n' side effects stored in the column.
}
\examples{
mat <- build.effects(
	sider.id = "sider_drugNames.tsv",
	sider.effects = "sider_meddra_all_se.tsv",
	kegg.filter = "kegg_genesList.txt", 
	drugbank = "drugbank_all_full_database-5.1.1.xml.gz", 
	drugbank.id = "uniprot_links.csv.gz",
	drugbank.filter = "inhibit",
	cutoff = 20)
}
\keyword{drugbank, sider}