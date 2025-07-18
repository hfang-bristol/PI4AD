% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oSNP2nGenes.r
\name{oSNP2nGenes}
\alias{oSNP2nGenes}
\title{Function to define nearby genes given a list of SNPs}
\usage{
oSNP2nGenes(
data,
distance.max = 2e+05,
decay.kernel = c("rapid", "slow", "linear", "constant"),
decay.exponent = 2,
GR.SNP = "dbSNP_Common",
GR.Gene = "UCSC_knownGene",
include.TAD = c("none", "GM12878", "IMR90", "MSC", "TRO", "H1", "MES",
"NPC"),
verbose = TRUE,
placeholder = NULL,
guid = NULL
)
}
\arguments{
\item{data}{an input vector containing SNPs. SNPs should be provided as
dbSNP ID (ie starting with rs). Alternatively, they can be in the
format of 'chrN:xxx', where N is either 1-22 or X, xxx is genomic
positional number; for example, 'chr16:28525386'}

\item{distance.max}{the maximum distance between genes and SNPs. Only
those genes no far way from this distance will be considered as seed
genes. This parameter will influence the distance-component weights
calculated for nearby SNPs per gene}

\item{decay.kernel}{a character specifying a decay kernel function. It
can be one of 'slow' for slow decay, 'linear' for linear decay, and
'rapid' for rapid decay. If no distance weight is used, please select
'constant'}

\item{decay.exponent}{a numeric specifying a decay exponent. By
default, it sets to 2}

\item{GR.SNP}{the genomic regions of SNPs. By default, it is
'dbSNP_Common', that is, Common SNPs from dbSNP (version 151) plus GWAS
SNPs and their LD SNPs (hg19). Alternatively, the user can specify the
customised input. To do so, first save your RData file (containing an
GR object) into your local computer, and make sure the GR object
content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData
file name (with or without extension), plus specify your file RData
path in "placeholder". Note: you can also load your customised GR
object directly}

\item{GR.Gene}{the genomic regions of genes. By default, it is
'UCSC_knownGene', that is, UCSC known genes (together with genomic
locations) based on human genome assembly hg19. Alternatively, the user
can specify the customised input. To do so, first save your RData file
(containing an GR object) into your local computer, and make sure the
GR object content names refer to Gene Symbols. Then, tell "GR.Gene"
with your RData file name (with or without extension), plus specify
your file RData path in "RData.location". Note: you can also load your
customised GR object directly}

\item{include.TAD}{TAD boundary regions are also included. By default,
it is 'none' to disable this option. Otherwise, inclusion of a TAD
dataset to pre-filter SNP-nGene pairs (i.e. only those within a TAD
region will be kept). TAD datasets can be one of "GM12878"
(lymphoblast), "IMR90" (fibroblast), "MSC" (mesenchymal stem cell)
,"TRO" (trophoblasts-like cell), "H1" (embryonic stem cell), "MES"
(mesendoderm) and "NPC" (neural progenitor cell). Explanations can be
found at \doi{10.1016/j.celrep.2016.10.061}}

\item{verbose}{logical to indicate whether the messages will be
displayed in the screen. By default, it sets to true for display}

\item{placeholder}{the characters to tell the location of built-in
RData files. See \code{\link{oRDS}} for details}

\item{guid}{a valid (5-character) Global Unique IDentifier for an OSF
project. See \code{\link{oRDS}} for details}
}
\value{
a data frame with following columns:
\itemize{
\item{\code{Gene}: nearby genes}
\item{\code{SNP}: SNPs}
\item{\code{Dist}: the genomic distance between the gene and the SNP}
\item{\code{Weight}: the distance weight based on the genomic distance}
\item{\code{Gap}: the genomic gap between the gene and the SNP (in the
form of 'chr:start-end')}
\item{\code{TAD}: if applied, it can be 'Excluded' or the TAD boundary
region (in the form of 'chr:start-end') that the genomic interval falls
into. Also if SNP within the gene body, Gap and TAD will be SNP
location (in the form of 'chr:start-end')}
}
}
\description{
\code{oSNP2nGenes} is supposed to define nearby genes given a list of
SNPs within certain distance window. The distance weight is calcualted
as a decaying function of the gene-to-SNP distance.
}
\examples{
\dontrun{

# a) provide the seed SNPs with the significance info
## get lead SNPs reported in AS GWAS and their significance info (p-values)
gr <- ImmunoBase$AS$variant
data <- names(gr)

# b) define nearby genes
df_nGenes <- oSNP2nGenes(data=data, distance.max=200000,
decay.kernel="slow", decay.exponent=2, RData.location=RData.location)

# c) define nearby genes (considering TAD boundary regions in GM12878)
df_nGenes <- oSNP2nGenes(data=data, distance.max=200000,
decay.kernel="slow", decay.exponent=2, include.TAD='GM12878',
RData.location=RData.location)
}
}
\seealso{
\code{\link{oSNP2nGenes}}
}
