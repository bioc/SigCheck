\name{sigCheckPlot}
\alias{sigCheckPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot results of a signature check
}
\description{
Plots a histogram of the classification performance scores for a check, showing 
how it compares to classification performance of the signature being checked, 
as well as to a hypothetical classifier that uses the mode of the training 
samples.
}
\usage{
sigCheckPlot(checkResults, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{checkResults}{
The list value returned by \code{\link{sigCheckRandom}}, 
\code{\link{sigCheckKnown}}, or \code{\link{sigCheckPermuted}} 
(or one the elements of the list returned by \code{\link{sigCheck}}).

Can also be the list returned by \code{\link{sigCheck}}, which which case 
each of the four results lists will be plotted in turn.
}
  \item{\dots}{
Additional arguments to be passed to the \code{\link{plot}} function.
}
}
\details{
Draws a line plot version of a histogram, with the x-axis representing the 
range of classification performance scores computed in the check, and the
y-axis representing how many times that score was obtained. In addition, 
vertical lines are plotted representing the classification performance of
the originally specified signature (solid red line) and the performance of 
a classifier that always predicts the mode value of the training samples
(dotted red line).
}
\value{
none
}

\note{
By default, \code{\link{sigCheck}} will call this function for all 
checks it runs.
}
\author{
Rory Stark with Justin Norden
}

\seealso{
\code{\link{sigCheck}},  \code{\link{sigCheckRandom}}, 
\code{\link{sigCheckKnown}}, \code{\link{sigCheckPermuted}}
}

\examples{
library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
results <- sigCheckRandom(nki, classes="e.dmfs", 
                          signature=knownSignatures$cancer$VANTVEER, 
                          annotation="HUGO.gene.symbol", 
                          validationSamples=275:319, nIterations=25)
sigCheckPlot(results)
}
