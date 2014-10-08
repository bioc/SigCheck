\name{sigCheck}
\alias{sigCheck}
\title{
Check classification potential of a gene signature against randomly selected 
gene signatures,  known gene signatures, and permuted expression sets.
}
\description{
High-level function for package \code{\link{SigCheck}} that runs all available 
checks against a classification signature.
}
\usage{
sigCheck(expressionSet, classes, signature, annotation, validationSamples, 
         classifierMethod = svmI, nIterations = 10, knownSignatures="cancer", 
         plotResults=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expressionSet}{
An \code{\link{ExpressionSet}} object containing the data to be checked, 
ncluding an expression matrix, feature labels, and samples.
}
  \item{classes}{
Specifies which label is to be used to determine the classification categories 
(must be one of varLabels(expressionSet)). There should be only 
two unique values in expressionSet$classes.
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
\code{expressionSet} should be used for validation. If present, a classifier 
will be trained, using the specified signature and classification method, on 
the non-validation samples, and it's performance evaluated by attempting to 
classify the validations samples. If missing, a leave-one-out (LOO) validation
method will be used, where a separate classifier will be trained to classify 
each sample using the reaming samples.
}
  \item{classifierMethod}{
The MLInterfaces learnerSchema object indicating the machine learning method 
to use for classification. Default is \code{\link{svmI}} for linear 
Support Vector Machine classification.  See \code{\link{MLearn}} for
available methods.
}
  \item{nIterations}{
For random gene and permutation tests, the number of iterations to run to 
compare classification outcomes. 
}
 \item{knownSignatures}{
Either a character string specifying which set of signatures to use from the 
included sets in \code{\link{knownSignatures}}, or a list of previously 
identified signatures to compare performance against. Each element in the 
list should be a vector of feature labels. Default is to use the 
\code{"cancer"} signatures from the included \code{\link{knownSignatures}} 
data set, taken from Venet et. al.
}
 \item{plotResults}{
if \code{TRUE}, will call \code{\link{sigCheckPlot}} four times to plot the 
results of all checks (laid out in a 2x2 plot matrix).
}
}
\details{
First, \code{sigCheck} calls \code{\link{sigCheckClassifier}} to establish the 
baseline performance of the signature being checked.

Next, it calls \code{\link{sigCheckRandom}} to check the performance of randomly 
selected signatures. 

This is followed by a call to \code{\link{sigCheckKnown}} to check the 
performance of the signature against a database of signatures previously 
identified to discriminate in other, generally unrelated domains.

Finally, two calls are made to \code{\link{sigCheckPermuted}} to check the 
performance of randomly permuted data; the first call permutes the rows
(\code{toPermute="features"}), while the second call permutes the categories 
(\code{toPermute="categories"}). 
}
\value{
A list containing five elements:
\itemize{
\item \code{$checkClassifier} is the result list returned by 
\code{\link{sigCheckClassifier}}.

\item \code{$checkRandom} is the result list returned by 
\code{\link{sigCheckRandom}}.

\item \code{$checkKnown} is the result list returned by 
\code{\link{sigCheckKnown}}.

\item \code{$checkPermutedFeatures} is the result list returned by 
\code{\link{sigCheckPermuted}}
with \code{toPermute="features"}.

\item \code{$checkPermutedCategories} is the result list returned by 
\code{\link{sigCheckPermuted}}
with \code{toPermute="categories"}.
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
\code{\link{sigCheckClassifier}}, \code{\link{sigCheckRandom}}, 
\code{\link{sigCheckPermuted}}, \code{\link{sigCheckKnown}}, 
\code{\link{MLearn}}
}
\examples{
library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)
results <- sigCheck(nki, classes="e.dmfs", 
                    annotation="HUGO.gene.symbol",
                    signature=knownSignatures$cancer$VANTVEER, 
                    validationSamples=275:319, nIterations=5)
}

