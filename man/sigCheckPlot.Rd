\name{sigCheckPlot}
\alias{sigCheckPlot}
\title{
Plot results of a signature check or set of checks
}
\description{
Plots results of a signature check, as returned by \code{\link{sigCheckRandom}},
\code{\link{sigCheckKnown}}, \code{\link{sigCheckPermuted}}, or 
\code{\link{sigCheckAll}}.
}

\usage{
sigCheckPlot(checkResults, classifier=FALSE, 
             title, nolegend=FALSE, ...)
}

\arguments{
  \item{checkResults}{
The list  returned by \code{\link{sigCheckRandom}}, 
\code{\link{sigCheckKnown}}, \code{\link{sigCheckPermuted}},
or \code{\link{sigCheck}}. 
}
\item{classifier}{
If a classifier was used in the original call to \code{\link{sigCheck}}, setting this
to \code{TRUE} will result in a plot showing how the primary signature compares
based on classification performance (rather than survival).
}
\item{title}{
Title string for plot. If missing, a default plot title will be generated.
}
\item{nolegend}{
IF \code{TRUE}, no legend will be included with the plot(s).
}

  \item{\dots}{
Additional arguments to be passed to the \code{\link{plot}} function.
}
}
\details{
For results based on survival analysis, the background distribution of p-values (in \code{1-log10()} format)
derived from the check (either random signatures, known signatures, or performance
using permuted data) is plotted. Up to two vertical red lines are also plotted: a solid red line representing the performance of the primary signature/data, 
and a dotted red line representing a p-value of 0.05. 
One or both of these may be missing if their performance falls outside the range
of the background distributions.

For results based on classification performance, the x-axis represents the 
range of classification performance scores computed in the check, and the
y-axis representing how many times that score was obtained. In addition, 
vertical lines are plotted representing the classification performance of
the originally specified signature (solid red line) and the performance of 
a classifier that always predicts the mode value of the training samples
(dotted red line).

If the results of \code{\link{sigCheckAll}} is passed in, all four results
plots are generated in a 2x2 grid.
}
\value{
none
}

\note{
Better formance of the signature being checked results in the solid red line 
being to the right of the background distribution.
For survival results, this indicates a lower-value.
For classification results,this indicates superior classification performance.
}
\author{
Rory Stark with Justin Norden
}

\seealso{
\code{\link{sigCheck}}, 
\code{\link{sigCheckAll}},
\code{\link{sigCheckRandom}}, 
\code{\link{sigCheckKnown}}, 
\code{\link{sigCheckPermuted}}
}

\examples{
#Disable parallel so Bioconductor build won't hang
library(BiocParallel)
register(SerialParam())

library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
ITERATIONS <- 5 # should be at least 1000 for real checks

## survival analysis with separate training and validation using SVM
check <- sigCheck(nki, classes="e.dmfs", survival="t.dmfs",
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol",
                  validationSamples=250:319,
                  scoreMethod="classifier", threshold=.33) 
                  
results <- sigCheckRandom(check,iterations=ITERATIONS)
par(mfrow=c(1,2))
sigCheckPlot(results)
sigCheckPlot(results, classifier=TRUE)
}
