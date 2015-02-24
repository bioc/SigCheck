\name{knownSignatures}
\alias{knownSignatures}
\alias{resultsNKI}
\docType{data}
\title{
Previously identified gene signatures for use in \code{\link{sigCheckKnown}}
}
\description{
Previously identified gene signature sets. These include three signatures sets 
from Venet et. al.
}
\usage{data(knownSignatures)}
\format{
The data object \code{knownSignatures} is a list of sets of gene signatures. 
Each set is a list of gene signatures. 
Each signature is a vector of gene names. 

Gene signature sets include:

\itemize{

\item \code{"cancer"}: 48 gene signatures derived from cancer samples, 
from Venet et. al.

\item \code{"proliferation"}: 5 gene signatures comprising genes associated 
with cell proliferation, including a "super signature", from Venet et. al.

\item \code{"non.cancer"}: 3 gene signatures derived from non-cancer sources, 
from Venet et. al.
}
}

\details{
These data are taken directly from the supplemental material for Venet et. al 
"Most random gene expression signatures are significantly associated with 
breast cancer outcome".
}

\source{
http://www.ploscompbiol.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pcbi.1002240.s001
}

\references{
Venet, David, Jacques E. Dumont, and Vincent Detours.
"Most random gene expression signatures are significantly associated with 
breast cancer outcome." PLoS Computational Biology 7.10 (2011): e1002240.
}

\note{
Other signatures of interest can be downlaoded at http://www.broad.mit.edu/gsea/downloads.jsp#msigdb.
}

\examples{
data(knownSignatures)
names(knownSignatures)
names(knownSignatures$cancer)
knownSignatures$cancer$VANTVEER
}

\keyword{datasets}
