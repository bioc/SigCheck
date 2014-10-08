\name{knownSignatures}
\alias{knownSignatures}
\alias{resultsNKI}
\docType{data}
\title{
Previously identified gene signatures for use in \code{\link{sigCheckKnown}}
}
\description{
Previously identified gene signature sets. These include three signatures sets 
from Venet et. al. and two signature sets from the 
Molecular Signatures Database (MSigDB) (Subramanian et. al).
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

\item \code{"c2.curated"}: 4,722 curated gene signatures downloaded from the
Molecular Signatures Database (MSigDB) (Subramanian et. al).

\item \code{"c6.oncogenic"}: 189 gene signatures derived from oncogenic samples, 
downloaded from the Molecular Signatures Database (MSigDB) (Subramanian et. al).
}
}

\details{
These data are taken directly from the supplemental material for Venet et. al 
"Most random gene expression signatures are significantly associated with 
breast cancer outcome" and the Molecular Signatures Database.
}

\source{
http://www.ploscompbiol.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pcbi.1002240.s001
http://www.broadinstitute.org/gsea/msigdb/index.jsp
}

\references{
Venet, David, Jacques E. Dumont, and Vincent Detours.
"Most random gene expression signatures are significantly associated with 
breast cancer outcome." PLoS Computational Biology 7.10 (2011): e1002240.

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., 
Gillette, M. A., ... & Mesirov, J. P. (2005). Gene set enrichment analysis: 
a knowledge-based approach for interpreting genome-wide expression profiles. 
Proceedings of the National Academy of Sciences of the United States of America, 
102(43), 15545-15550.
}

\examples{
data(knownSignatures)
names(knownSignatures)
names(knownSignatures$cancer)
knownSignatures$cancer$VANTVEER
}

\keyword{datasets}
