\name{sigCheckRandom}
\alias{sigCheckRandom}
\title{
Check signature performance against signatures composed of randomly selected 
features
}
\description{
Performance of a signature is compared to performance of signatures composed of 
the same number of randomly-selected features.
}
\usage{
sigCheckRandom(check, iterations=10)
}

\arguments{
\item{check}{
A \code{\link{SigCheckObject}}, as returned by \code{\link{sigCheck}}.
}
  \item{iterations}{
The number of random signatures the primary signature will be compared to.
This should be at least 1,000 to compute a meaningful empirical p-value
for comparative performance.
}
}
\details{
\cite{sigCheckRandom} will select \code{iterations} signatures, each consisting
of the same number of features as are in the primary signature
provided in the call to \code{\link{sigCheck}} that created the 
\code{\link{SigCheckObject}} sampled at random from all available features.

Each random signature will be evaluated in the same manner as the primary
signature.
If survival data were supplied, a survival analysis will be carried out on the 
validation samples, and a p-value computed as a performance measure. 
If no survival data are available, the training samples will be used to train
a classifier, and the performance score will be percentage of validation samples
correctly classified. (If no validation samples are provided, leave-one-out
cross validation will be used to calculate the classification performance 
for each random signature).

An empirical p-value will be computed based on the percentile rank of the 
performance of the primary signature compared to a null distribution of
the performance of the random signatures.
}
\value{
A result list with the following elements:

\itemize{

\item \code{$checkType} is equal to \code{"Random"}.

\item \code{$tests} is the number of tests run (equal to \code{iterations}.)

\item \code{$rank} is the performance rank of the primary signature  
within the performance of the random signatures.

\item \code{$checkPval} is the empirical p-value computed using the scores
of the random signature as a null distribution. A value of zero indicates that
no random signatures performed as good or better than the primary signature.

\item \code{$survivalPval} represents the performance of the primary,
if survival data were provided.

\item \code{$survivalPvalsRandom} is a vector of performance scores (p-values) 
for each random signature on the validation samples, if survival data
were provided.

\item \code{$trainingPvalsRandom} is a vector of performance scores (p-values) 
for each random signature on the training samples, if survival data
and separate validation samples were provided.

\item \code{$sigPerformance} is the proportion of validation samples 
correctly classified by the primary signature if a classifier was used.

\item \code{$modePerformance} is the proportion of validation samples 
correctly classified using a mode classifier.

\item \code{$performanceRandom} is a vector of classification performance 
scores for each random signature, each indicating the proportion
of validation samples correctly classified if a classifier was used.
}
}

\author{
Rory Stark
}

\seealso{
\code{\link{sigCheck}}, \code{\link{sigCheckAll}}, 
\code{\link{sigCheckPermuted}}, \code{\link{sigCheckKnown}},
\code{\link{sigCheckPlot}}
}
\examples{
library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
ITERATIONS <- 5 # should be at least 20, 1000 for real checks

## survival analysis
check <- sigCheck(nki, classes="e.dmfs", survival="t.dmfs",
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol",
                  validationSamples=150:319)

randomResult <- sigCheckRandom(check, iterations=ITERATIONS)
randomResult$checkPval
sigCheckPlot(randomResult)
}

