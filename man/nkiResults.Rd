\name{nkiResults}
\alias{nkiResults}
\docType{data}
\title{
Precomputed list of results for a call to \code{\link{sigCheckAll}} using the 
\code{breastCancerNKI} dataset.
}
\description{
This object represents the results lists returned by a call to 
\code{\link{sigCheckAll}}. It is used by the vignette accompanying the 
\code{\link{SigCheck}} package as an example result. It was derived by 
running the code in the example below.
}
\usage{data(nkiResults)}

\examples{
\dontrun{
# This is how nkiResults is built
library(breastCancerNKI)
data(nki)
nki = nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
check <- sigCheck(nki, classes="e.dmfs", survival="t.dmfs",
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol",
                  validationSamples=which(nki$series=="NKI2"))
                  
nkiResults <- sigCheckAll(check,iterations=1000) 
}

# Example usage of nkiResults
data(nkiResults)
sigCheckPlot(nkiResults)
}
\keyword{datasets}
