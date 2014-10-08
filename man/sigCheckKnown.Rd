\name{sigCheckKnown}
\alias{sigCheckKnown}
\title{
Check classification performance of signature against a panel of 
known gene signatures
}
\description{
Compare the classification performance of a known panel of gene signatures to 
the signature being checked. By default, a panel of gene signatures from 
Venet et. al. is used.
}
\usage{
sigCheckKnown(expressionSet, classes, signature, annotation, validationSamples, 
              classifierMethod = svmI, classifierScore, knownSignatures="cancer")
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
used as the annotation. If missing, the row names of the expressionSet are 
used as the feature names.
}
  \item{validationSamples}{
Optional specification, as a vector of sample indices, of what samples in the  
should used for validation. If present, a classifier will be trained, using the 
specified signature and classification method, on the non-validation samples, 
and it's performance evaluated by attempting to classify the validations samples. 
If missing, a leave-one-out (LOO) validation method will be used, where a 
separate classifier will be trained to classify each sample using the remaining
samples.
}
  \item{classifierMethod}{
The MLInterfaces learnerSchema object indicating the machine learning method to
use for classification. Default is \code{\link{svmI}} for linear
Support Vector Machine classification.  See \code{\link{MLearn}} for
available methods.
}

\item{classifierScore}{
A performance measure of the baseline classifier. Generally the 
\code{classifierScore} element of the result list returned by 
\code{\link{sigCheckClassifier}}. 
If missing, \code{\link{sigCheckClassifier}} will be called to establish 
baseline performance.
}
 \item{knownSignatures}{
Either a character string specifying which set of signatures to use from the 
included sets in \code{\link{knownSignatures}}, or a list of previously 
identified signatures to compare performance against. Each element in the 
list should be a vector of feature labels. Default is to use the 
\code{"cancer"} signatures from the included \code{\link{knownSignatures}} 
data set, taken from Venet et. al.
}
}
\details{
\code{\link{sigCheckClassifier}} is called for each of the known signatures.
}
\value{
A list with six elements:
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

\item \code{$known} is a character string indicating which gene signature set 
was checked. Either one of the sets in \code{\link{knownSignatures}}, or the 
string \code{"user specified"}.

\item \code{$knownSigs} is the number of signatures evaluated (equal to 
\code{length(knownSignatures)}, minus any signatures with zero features 
matching the labels in \code{expressionSet}.)

\item \code{$rank} is the performance rank of the primary signature classifier 
on the original dataset amongst the performances of the known signatures 
on the same dataset.

\item \code{$performanceKnown} is a vector of performance scores (proportion 
of the validation set correctly predicted) for each known signature 
on the dataset.
}
}
\references{
Venet, David, Jacques E. Dumont, and Vincent Detours. 
"Most random gene expression signatures are significantly associated 
with breast cancer outcome." PLoS Computational Biology 7.10 (2011): e1002240.
}

\author{
Justin Norden with Rory Stark
}

\seealso{
\code{\link{knownSignatures}}, \code{\link{sigCheck}}, 
\code{\link{sigCheckClassifier}}, \code{\link{sigCheckRandom}}, 
\code{\link{sigCheckPermuted}}, \code{\link{MLearn}}
}
\examples{
library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
results <- sigCheckKnown(nki, classes="e.dmfs", 
                         signature=knownSignatures$cancer$VANTVEER, 
                         annotation="HUGO.gene.symbol", 
                         validationSamples=275:319)
}

