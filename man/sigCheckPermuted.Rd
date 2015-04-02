\name{sigCheckPermuted}
\alias{sigCheckPermuted}
\title{
Check signature performance against performance on randomly permuted data.
}
\description{
Performance of a  signature on intact data is compared to
performance on permuted data and/or metadata. Data may be permuted by
by feature (expression values of each 
feature permuted across samples), samples (expression values of all features 
permuted within each sample). Metadata may be permuted by categories 
(permuted assignment of samples to 
classification categories) or survival (permuted assigment
of survival times to samples).
}
\usage{
sigCheckPermuted(check, toPermute="categories", iterations=10)
}


\arguments{
\item{check}{
A \code{\link{SigCheckObject}}, as returned by \code{\link{sigCheck}}.
}

  \item{toPermute}{
Character string or vector of strings indicating what should be permuted. 
Allowable values:
\itemize{
\item \code{"features"}: the expression values for each feature will be 
permuted (permutation by row).

\item \code{"samples"}: the expression values for each sample will be 
permuted (permutation by column).

\item \code{"survival"}: the values in \code{survival} will be permuted.

\item \code{"categories"}: the values in \code{classes} will be permuted.
}}

  \item{iterations}{
The number of permuted dataset the primary signature will be compared to.
This should be at least 1,000 to compute a meaningful empirical p-value
for comparative performance.
}
}

 
\details{
The primary signature  will be evaluated against each permuted dataset in the 
same manner as for the intact dataset.

If survival data were supplied, a survival analysis will be carried out on the 
validation samples, and a p-value computed as a performance measure. 
If no survival data are available, the training samples will be used to train
a classifier, and the performance score will be percentage of validation samples
correctly classified. (If no validation samples are provided, leave-one-out
cross validation will be used to calculate the classification performance 
for each permuted datdaset).

An empirical p-value will be computed based on the percentile rank of the 
performance of the  signature on the intact datdaset compared to a 
null distribution of
the performance of the signature on all the permuted datasets.
}

\value{
A result list with the following elements:

\itemize{

\item \code{$checkType} is equal to \code{"Permuted"}.

\item \code{$permute} is equal to the passed value of \code{toPermute}.

\item \code{$tests} is the number of tests run (equal to \code{iterations}.)

\item \code{$rank} is the performance rank of the signature on the 
intact dataset compared to its performance in the permuted datasets.

\item \code{$checkPval} is the empirical p-value computed using the performance
scores of the signature on permuted datdasets as a null distribution. 
A value of zero indicates that the signature did not perform better on 
any permuted datasets than it does using the intact data.

\item \code{$survivalPval} represents the performance of the primary signature
on the original dataste
if survival data were provided.

\item \code{$survivalPvalsPermuted} is a vector of performance scores (p-values) 
for each permuted datdaset, if survival data
were provided.

\item \code{$trainingPvalsPermuted} is a vector of performance scores (p-values) 
for each permuted datdaset, if survival data
and separate validation samples were provided.

\item \code{$sigPerformance} is the proportion of validation samples 
correctly classified using th eintact datdaset if a classifier was used.

\item \code{$modePerformance} is the proportion of validation samples 
correctly classified in the intact dataaset using a mode classifier.

\item \code{$performancePermuted} is a vector of classification 
performance scores for each
permuted datdaset, each indicating the proportion
of validation samples correctly classified if a classifier was used.
}
}


\author{
Rory Stark
}

\seealso{
\code{\link{sigCheck}}, \code{\link{sigCheckAll}}, 
\code{\link{sigCheckRandom}}, \code{\link{sigCheckKnown}}
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

par(mfrow=c(1,2))
permutedCategories <- sigCheckPermuted(check, toPermute="categories", 
                                       iterations=ITERATIONS)
permutedCategories$checkPval
sigCheckPlot(permutedCategories)
permutedSurvival <- sigCheckPermuted(check, toPermute="survival", 
                                     iterations=ITERATIONS)
permutedSurvival$checkPval
sigCheckPlot(permutedSurvival)
}

