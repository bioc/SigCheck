\name{classifyResults}
\alias{classifyResults}
\docType{data}
\title{
Precomputed list of results for a calls to \code{\link{sigCheckRandom}}
and \code{\link{sigCheckKnown}} using the 
\code{breastCancerNKI} dataset.
}
\description{
This object represents the results lists returned by calls to 
\code{\link{sigCheckRandom}} and \code{\link{sigCheckKnown}}. 
It is used by the vignette accompanying the 
\code{\link{SigCheck}} package as an example result. It was derived by 
running the code in the example below.
It loads to object, \code{slassifyRandom} and \code{classifyKnown}.
}
\usage{data(classifyResults)}

\examples{
\dontrun{
# This is how classifyResults is built
library(breastCancerNKI)
data(nki)
nki = nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
check <- sigCheck(nki, classes="e.dmfs",
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol",
                  validationSamples=101:ncol(nki),
                  scoreMethod="classifier")
                  
classifyRandom <- sigCheckRandom(check, iterations=1000)
classifyKnown <- sigCheckKnown(check)
}

# Example usage of classifyResults
data(classifyResults)
par(mfrow=c(1,2))
sigCheckPlot(classifyRandom, classifier=TRUE)
sigCheckPlot(classifyKnown, classifier=TRUE)

}
\keyword{datasets}
