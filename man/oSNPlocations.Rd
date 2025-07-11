% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oSNPlocations.r
\name{oSNPlocations}
\alias{oSNPlocations}
\title{Function to extract genomic locations given a list of SNPs}
\usage{
oSNPlocations(
data,
GR.SNP = "dbSNP_Common",
verbose = TRUE,
placeholder = NULL,
guid = NULL
)
}
\arguments{
\item{data}{a input vector containing SNPs. SNPs should be provided as
dbSNP ID (ie starting with rs). Alternatively, they can be in the
format of 'chrN:xxx', where N is either 1-22 or X, xxx is genomic
positional number; for example, 'chr16:28525386'}

\item{GR.SNP}{the genomic regions of SNPs. By default, it is
'dbSNP_Common', that is, Common SNPs from dbSNP (version 151) plus GWAS
SNPs and their LD SNPs (hg19). Alternatively, the user can specify the
customised input. To do so, first save your RData file (containing an
GR object) into your local computer, and make sure the GR object
content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData
file name (with or without extension), plus specify your file RData
path in "placeholder". Note: you can also load your customised GR
object directly}

\item{verbose}{logical to indicate whether the messages will be
displayed in the screen. By default, it sets to true for display}

\item{placeholder}{the characters to tell the location of built-in
RData files. See \code{\link{oRDS}} for details}

\item{guid}{a valid (5-character) Global Unique IDentifier for an OSF
project. See \code{\link{oRDS}} for details}
}
\value{
an GR oject, with an additional metadata column called 'variant_id'
storing SNP location in the format of 'chrN:xxx', where N is either
1-22 or X, xxx is genomic positional number.
}
\description{
\code{oSNPlocations} is supposed to extract genomic locations given a
list of SNPs.
}
\note{
none
}
\examples{
\dontrun{
# a) provide the SNPs
## get lead SNPs reported in AS GWAS and their significance info (p-values)
gr <- ImmunoBase$AS$variant
data <- names(gr)

# b) find the location
snp_gr <- oSNPlocations(data=data, RData.location=RData.location)
}
}
\seealso{
\code{\link{oSNPlocations}}
}
