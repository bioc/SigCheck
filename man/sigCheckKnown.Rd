\name{sigCheckKnown}
\alias{sigCheckKnown}
\title{
Check signature performance against a panel of known
signatures.
}
\description{
Performance of a  signature is compared to performance of
a panel of
known (previously identified) signature.
}
\usage{
sigCheckKnown(check, known="cancer")
}

\arguments{
\item{check}{
A \code{\link{SigCheckObject}}, as returned by \code{\link{sigCheck}}.
}
  \item{known}{
Either a character string specifying which set of signatures to use from the 
included sets in \code{\link{knownSignatures}}, 
or a list of previously identified signatures to compare performance against. 
Each element in the list should be a vector of feature labels. 
Default is to use the "cancer" signatures from the included k
\code{\link{knownSignatures}} data set, taken from Venet et. al.
}
}

\details{
Each specified known signature will be evaluated in the 
same manner as the primary signature.
If survival data were supplied, a survival analysis will be carried out on the 
validation samples, and a p-value computed as a performance measure. 
If no survival data are available, the training samples will be used to train
a classifier, and the performance score will be percentage of validation samples
correctly classified. (If no validation samples are provided, leave-one-out
cross validation will be used to calculate the classification performance 
for each known signature).

An empirical p-value will be computed based on the percentile rank of the 
performance of the primary signature compared to a null distribution of
the performance of the known signatures.
}

\value{
A result list with the following elements:

\itemize{

\item \code{$checkType} is equal to \code{"Known"}.

\item \code{$knownSigs} is the number of tests run (equal to the
number of known signatures indicated where at least one gene matches
a feature).

\item \code{$rank} is the performance rank of the primary signature  
within the performance of the known signatures.

\item \code{$checkPval} is the empirical p-value computed using the scores
of the known signature as a null distribution. A value of zero indicates that
no known signatures performed as good or better than the primary signature.

\item \code{$survivalPval} represents the performance of the primary signature,
if survival data were provided.

\item \code{$survivalPvalsKnown} is a vector of performance scores (p-values) 
for each known signature on the validation samples, if survival data
were provided.

\item \code{$trainingPvalsKnown} is a vector of performance scores (p-values) 
for each known signature on the training samples, if survival data
and separate validation samples were provided.

\item \code{$sigPerformance} is the proportion of validation samples 
correctly classified by the primary signature if a classifier was used.

\item \code{$modePerformance} is the proportion of validation samples 
correctly classified using a mode classifier.

\item \code{$performanceKnown} is a vector of classification performance
scores for each
known signature, each indicating the proportion
of validation samples correctly classified is a classifier was used.
}
}

\references{
Venet, David, Jacques E. Dumont, and Vincent Detours. 
"Most random gene expression signatures are significantly associated 
with breast cancer outcome." PLoS Computational Biology 7.10 (2011): e1002240.
}

\author{
Rory Stark
}

\seealso{
\code{\link{knownSignatures}}, \code{\link{sigCheck}}, 
\code{\link{sigCheckAll}}, \code{\link{sigCheckRandom}}, 
\code{\link{sigCheckPermuted}}, \code{\link{sigCheckPlot}}
}
\examples{
library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)

## survival analysis
check <- sigCheck(nki, classes="e.dmfs", survival="t.dmfs",
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol",
                  validationSamples=150:319)

knownResult <- sigCheckKnown(check)
knownResult$checkPval
knownResult$survivalPvalsKnown[knownResult$survivalPvalsKnown <
                               knownResult$checkPval]
sigCheckPlot(knownResult)
}

