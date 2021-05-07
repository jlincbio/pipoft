\name{prep.setup}
\alias{prep.setup}
\title{Package Configuration for 'pipoft'}
\usage{
prep.setup(perl = NULL, cpan = NULL, silent = TRUE)
}
\arguments{
\item{perl}{Path to Perl interpreter; set to 'NULL' to use system defaults.}
\item{cpan}{Path to CPAN, to install missing packages. Currently unused.}
\item{silent}{Should success messages be suppressed? Defaults to 'TRUE'.}
}
\description{
This call is automatically invoked during the installation of 'pipoft', and will fail when missing Perl packages are found. There is no need to call this if the package is installed successfully. In situations where missing packages are present, please install those Perl packages.
}
\keyword{pipoft}