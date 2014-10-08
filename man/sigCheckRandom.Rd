\name{sigCheckRandom}
\alias{sigCheckRandom}
\title{
Check classification performance of signatures composed of randomly selected 
features
}
\description{
Performance of a classification signature is compared to signatures composed of 
the same number of randomly-selected features.
}
\usage{
sigCheckRandom(expressionSet, classes, signature, 
               annotation, validationSamples, 
               classifierMethod = svmI, nIterations = 10, classifierScore)
}

\arguments{
  \item{expressionSet}{
An \code{\link{ExpressionSet}} object containing the data to be checked, 
including an expression matrix, feature labels, and samples.
}
  \item{classes}{
Specifies which label is to be used to determine the classification categories 
(must be one of \code{varLabels(expressionSet)}). There should be only two 
unique values in \code{expressionSet$classes}.
}
  \item{signature}{
A vector of feature labels specifying which features comprise the signature to 
be checked. These feature labels should match values as specified in the 
\code{annotation} parameter (default is row names in the expressionSet). 
Alternatively, this can be a integer vector of feature indexes.
}
  \item{annotation}{
Character string specifying which \code{\link{featureData}} field should be 
used as the annotation. If missing, the row names of the expressionSet are used
as the feature names.
}
  \item{validationSamples}{
Optional specification, as a vector of sample indices, of what samples  
should used for validation. If present, a classifier will be trained, using the 
specified signature and classification method, on the non-validation samples,
and it's performance evaluated by attempting to classify the validations samples. 
If missing, a leave-one-out (LOO) validation method will be used, where a 
separate classifier will be trained to classify each sample using the 
remaining samples.
}
  \item{classifierMethod}{
The MLInterfaces learnerSchema object indicating the machine learning method to 
use for classification. Default is \code{\link{svmI}} for linear 
Support Vector Machine classification.  See \code{\link{MLearn}} 
for available methods.
}
  \item{nIterations}{
The number of permutations to test and compare classification outcomes. 
}
\item{classifierScore}{
A performance measure of the baseline classifier. Generally the 
\code{classifierScore} element of the result list returned by 
\code{\link{sigCheckClassifier}}. If missing, \code{\link{sigCheckClassifier}} 
will be called to establish baseline performance.
}
}
\details{
First, the number of features in the passed signature that match features in 
the dataset is calculated. Next, \code{nIterations} signatures are generated 
and tested, each consisting of the same number of randomly selected features. 
Performance for each signature is determined by calling 
\code{\link{sigCheckClassifier}}.
}
\value{
A list with five elements:

\itemize{
\item \code{$sigPerformance} is the percentage of validationSamples correctly
classified (or, in the LOO case, the percentage of total samples correctly 
classified by classifiers trained using the remaining samples.)

\item \code{$modePerformance} is the percentage of validationSamples correctly 
classified by a "mode" classifier (or, in the LOO case, the percentage of total 
samples correctly classified by a "mode" classifier, which is equal the number 
of samples with the more-frequent category.) The "mode" classifier always 
predicts the category that appears most often in the training set. 
If the training set is balanced between categories, one category will 
always be predicted.

\item \code{$tests} is the number of tests run (equal to \code{nIterations}.)

\item \code{$rank} is the performance rank of the primary signature classifier 
amongst the performance of the random signatures.

\item \code{$performanceRandom} is a vector of performance scores (proportion 
of the validation set correctly predicted) for each random signature.
}
}

\author{
Justin Norden with Rory Stark
}

\seealso{
\code{\link{sigCheck}}, \code{\link{sigCheckClassifier}}, 
\code{\link{sigCheckPermuted}}, \code{\link{sigCheckKnown}}, 
\code{\link{MLearn}}
}
\examples{
library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
results <- sigCheckRandom(nki, classes="e.dmfs", 
                          signature=knownSignatures$cancer$VANTVEER, 
                          annotation="HUGO.gene.symbol", 
                          validationSamples=275:319)
}

