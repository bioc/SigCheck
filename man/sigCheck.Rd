\name{sigCheck}
\alias{sigCheck}

\title{
Create a \code{\link{SigCheckObject}} and establish baseline performance.
}

\description{
Main constructor for a \code{\link{SigCheckObject}}. Also establishes 
baseline survival analysis and/or classification performance.
}

\usage{
sigCheck(expressionSet, classes, survival, signature, 
         annotation, validationSamples, 
         scoreMethod="PCA1", threshold=median,
         classifierMethod=svmI, modeVal,
         survivalLabel, timeLabel,
         plotTrainingKM=TRUE, plotValidationKM=TRUE,
         impute=TRUE)
}


\arguments{

\item{expressionSet}{
An \code{\link{ExpressionSet}} object containing the data to be checked, 
including an expression matrix, feature labels, and samples.

\code{expressionSet} can also be an existing \code{\link{SigCheckObject}},
in which case everything will be inherited from the passed object 
except the values for any specified parameters,
}

\item{classes}{
Specifies which label is to be used to determine the prognostic categories 
(must be one of \code{varLabels(expressionSet))}. 
There should be only 
two unique values in expressionSet$classes.
}

\item{survival}{
Specifies which label is to be used to determine survival times.
(must be one of \code{varLabels(expressionSet))}.
This may be missing if only classification is being checked.
}

\item{signature}{
A vector of feature labels specifying which features comprise the signature to 
be checked. These feature labels should match values as specified in the 
\code{annotation} parameter. 
Alternatively, this can be a integer vector of feature indexes.
}

\item{annotation}{
Character string specifying which \code{\link{fvarLabels}} field should be 
used as the annotation. If missing, the row names of the \code{expressionSet} 
are used as the feature names.
}

\item{validationSamples}{
Optional specification, as a vector of sample indices, of what samples in the 
\code{expressionSet} should be considered validation samples.
If present, the main checks will be performed using only these samples.
If a the \code{scoreMethod} parameter is equal to \code{"classifier"},
the remaining samples will 
be used as a training set to construct a classifier that will be used to 
separate the training samples.
If a classifier is used, and \code{validationSamples} is not specified,
a leave-one-out (LOO) validation
method will be used, where a separate classifier will be trained to classify 
each sample using the all the remaining samples.
}

\item{scoreMethod}{
specification of how the samples should be split into groups for 
survival analysis. If a character sting, one of the following values:
\itemize{ 
\item \code{"PCA1"}: default scoring method for separating validation
samples into groups
by taking the value of the first principal component of the expression
values in the signature for each sample.
\item \code{"High"}: score used for separating validation samples into groups
for each sample is the mean value over 
all the expression values in the signature for each sample.
\item \code{"classifer"}: score used for separating validation 
samples into groups is determined by a classifier specified in the
\code{classifierMethod} parameter. If the \code{survival} parameter is 
specified, the classifier method must return a real-valued score for each
predicted sample.
}

\code{scoreMethod} can also be a user-defined function that computes a score.
The function should take a single parameter, an \code{\link{ExpressionSet}}, and
return a vector of score, one for each row.

if the \code{survival} parameter is missing, \code{scoreMethod} value must be 
\code{"classifier"}.
}
\item{threshold}{
specifies the threshold used for separating the validation samples into classed
based on the score derived using \code{scoreMethod}. 
Can be either a function, (with default \code{median}) or a number between 
zero and one indicating a percentile. Validation samples will be divided into
a group whose percentile scores are less than this value, and another 
group with percentile scores greater 
than or equal to this value.
\code{threshold} may also be a vector of two percentiles,
in which case samples
will be divided into High, Low, and Mid groups.
The survival p-value will be computed using only the high and low
groups, with the mid group samples excluded.
}

  \item{classifierMethod}{
  if the \code{scoreMethod} is equal to \code{"classifer"}, this specifies 
  what classifier to use. It is a \code{MLInterfaces} \code{learnerSchema}
  object indicating the machine learning method 
to use for classification. Default is \code{\link{svmI}} for linear 
Support Vector Machine classification.  See \code{\link{MLearn}} for
available methods.
}
  \item{modeVal}{
specifies which of the two category values (one of the values implied 
by the \code{classes} parameter) should be considered as the default
value when computing the performance of a "mode" classifier. Is missing,
the actual mode (most commonly occurring) value of the training set will be used.
}

\item{survivalLabel}{
String to use in the Y-axis of any Kaplan-Meier plots generated, this indicates
what aspect of survival is being predicted, such as time to recurrence or death.
}
\item{timeLabel}{
String to use in the X-axis of an Kaplan-Meier plots generated, this indicates
the units of time, such as days or months to outcome event.
}
\item{plotTrainingKM}{
if the \code{survival} and \code{validationSamples} parameters are provided,
a Kaplan-Meier plot can be plotted automatically for the training set samples if
this is \code{TRUE}. A value of \code{FALSE} will suppress the plot being
automatically generated.
}
\item{plotValidationKM}{
if the \code{survival} parameter is provided,
a Kaplan-Meier plot can be plotted automatically for the validation set samples if
this is \code{TRUE}. A value of \code{FALSE} will suppress the plot being
automatically generated.
Note that is the \code{validationSamples} parameter is missing, the resulting
plot will be over all samples.
}
\item{impute}{
if \code{TRUE}, missing data values in the \code{expressionSet} will be imputed.
If \code{FALSE}, any features with any missing values will be removed 
from the dataset.
}
}

\details{
This function constructs a new \code{\link{SigCheckObject}} and carried out
a baseline analysis, 
which will vary depending on which parameters are specified.

If the \code{survival} parameter is specified, a survival analysis
is carried out. 
If the \code{validationSamples} parameter is specified, this will be done 
separately on the validation samples and the remaining 
(training/discovery) samples. 
The main result is a p-value indicating the confidence that the samples are
separable into groups with distinct survival outcomes. This value is obtained
using the \code{\link{survdiff}} function in the \code{survival} package 
(and applying \code{\link{pchisq}} to the 
\code{$chisq} component of the result). The samples are separated into groups
using the \code{scoreMethod} and \code{threshold} parameters 
(and possibly the \code{classifierMethod} parameter).

If the \code{survival} parameter is not specified, then the \code{scoreMethod} 
parameter must be equal to \code{"classifier"}, and a pure classification 
analysis is completed (as was done in \code{SigCheck 1.0}). 
If the \code{validationSamples} parameter is specified, the remaining samples 
are used as a training set to construct a classifier that is used
to classify the validation samples. If \code{validationSamples} is not
specified, leave-one-out cross-validation is used whereby a separate
classifier is trained to predict each sample using all of the others.

}
\value{
If the baseline analysis can be completed, 
a \code{\link{SigCheckObject}} is returned.
}

\author{
Rory Stark with Justin Norden
}

\seealso{
\code{\link{sigCheckAll}}, \code{\link{sigCheckRandom}}, 
 \code{\link{sigCheckKnown}}, \code{\link{sigCheckPermuted}}.
}
\examples{
library(breastCancerNKI)
data(nki)
nki <- nki[,!is.na(nki$e.dmfs)]
data(knownSignatures)

## survival analysis 
check <- sigCheck(nki, classes="e.dmfs", survival="t.dmfs",
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol") 
check@survivalPval
check <- sigCheck(check, classes="e.dmfs", survival="t.dmfs",
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol",
                  scoreMethod="High", threshold=.33) 
check@survivalPval

## survival analysis with separate training and validation using SVM
check <- sigCheck(nki, classes="e.dmfs", survival="t.dmfs",
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol",
                  validationSamples=150:319,
                  scoreMethod="classifier") 
check
}

