\name{nkiResults}
\alias{nkiResults}
\docType{data}
\title{
Precomputed list of results for a call to \code{\link{sigCheck}} using the 
\code{breastCancerNKI} dataset.
}
\description{
This object represents the results lists returned by a call to 
\code{\link{sigCheck}}. It is used by the vignette accompanying the 
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
nkiResults = sigCheck(nki, classes="e.dmfs", annotation="HUGO.gene.symbol",
                   signature=knownSignatures$cancer$VANTVEER,
                   validationSamples=275:319, 
                   knownSignatures="cancer",nIterations=1000)
}

# Example usage of nkiResults
data(nkiResults)
nkiResults$checkClassifier
sigCheckPlot(nkiResults$checkRandom)
}
\keyword{datasets}
