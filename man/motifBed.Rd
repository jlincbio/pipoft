\name{prep.motif}
\alias{prep.motif}
\title{Generates a BED File of PI Polyamide Genomic Binding Locations by Motif}
\usage{
prep.motif(fasta, motif, dest = NULL)
}
\arguments{
\item{fasta}{FASTA sequence file to be searched against.}
\item{motif}{PI polyamide binding motif.}
\item{dest}{Should the results be stored somewhere? Set to NULL if in system temp.}
}
\description{
'prep.motif' will search a given FASTA file to identify possible PI polyamide binding sites (provided as a character string consisted of ACGTW bases). The resultant BED file can be used to estimate probable binding sites for other off-target functions or for Chem-seq applications. Use 'dest' to specify an filename if you want to save the results to a 4-column BED file.
}
\examples{
# finds binding sites for a PI polyamide against TGWWGGCGW in hg19
prep.motif(fasta = "hg19.fa", motif = "TGWWGGCGW")

}
\keyword{polyamide, off-targets}