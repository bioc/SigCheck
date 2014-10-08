\name{sigCheckClassifier}
\alias{sigCheckClassifier}
\title{
Establish baseline classification performance for a signature
}
\description{
Compute classification performance of a signature by training one or more 
classifiers and testing their ability to predict validation samples.
}
\usage{
sigCheckClassifier(expressionSet, classes, signature, annotation, 
                   validationSamples, classifierMethod = svmI, ...)
}
%- maybe also 'usage' for other objects documented here.
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
A vector of feature labels specifying which features comprise the signature to be checked. These feature labels should match values as specified in the \code{annotation} parameter (default is row names in the expressionSet). Alternatively, this can be a integer vector of feature indexes.
}
  \item{annotation}{
Character string specifying which \code{\link{featureData}} field should be
used as the annotation. If missing, the row names of the expressionSet are used as the feature names.
}
  \item{validationSamples}{
Optional specification, as a vector of sample indices, of what samples in the 
\code{expressionSet} should used for validation. If present, a classifier will 
be trained, using the specified signature and classification method, on the 
non-validation samples, and it's performance evaluated by attempting to 
classify the validations samples. If missing, a leave-one-out (LOO) validation 
method will be used, where a separate classifier will be trained to classify 
each sample using the remaining samples.
}
  \item{classifierMethod}{
The MLInterfaces learnerSchema object indicating the machine learning method to 
use for classification. Default is \code{\link{svmI}} for 
linear Support Vector Machine classification.  See \code{\link{MLearn}} 
for available methods.}

\item{\dots}{
additional parameters to be passed to \code{\link{MLearn}} in support of the 
classification method specified in \code{classifierMethod}.
}
}
\details{
If \code{validationSamples} are specified, the \code{MLInterfaces} package is
used to train a classifier on the remaining samples. By default, a 
Support Vector Machine classifier is used, but any machine learning approach 
supported by \code{\link{MLearn}} can be specified. Baseline performance is
measured by the percentage of the validation samples classified correctly
(a confusion matrix of the results is also returned). If the validationSamples 
are not specified, a leave-one-out (LOO) approach is deployed, whereby each
sample in turn is used as the validation sample, resulting in as many 
classifiers being trained as there are samples.
}
\value{
A list with three elements:
\itemize{
\item \code{$sigPerformance} is the percentage of validationSamples correctly 
classified (or, in the LOO case, the percentage of total samples correctly 
classified by classifiers trained using the remaining samples.)

\item \code{$confusion} is a confusion matrix in the form of a table showing 
how many samples in each class were correctly or incorrectly classified, 
corresponding to True Positives, True Negative, False Positives, 
and False Negatives.

\item \code{$modePerformance} is the percentage of validationSamples correctly 
classified by a "mode" classifier (or, in the LOO case, the percentage of total 
samples correctly classified by a "mode" classifier, which is equal the number 
of samples with the more-frequent category.) The "mode" classifier always 
predicts the category that appears most often in the training set. 
If the training set is balanced between categories, one category will 
always be predicted.
}
}
\author{
Justin Norden with Rory Stark
}

\seealso{
\code{\link{sigCheck}}, \code{\link{sigCheckRandom}}, 
\code{\link{sigCheckPermuted}}, \code{\link{sigCheckKnown}}, 
\code{\link{MLearn}}
}
\examples{
library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
results <- sigCheckClassifier(nki, classes="e.dmfs", 
                              signature=knownSignatures$cancer$VANTVEER, 
                              annotation="HUGO.gene.symbol")
}

