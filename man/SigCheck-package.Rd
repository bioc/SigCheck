\name{SigCheck-package}
\alias{SigCheck-package}
\alias{SigCheck}
\docType{package}
\title{
Check a gene signature's survival and/or classification performance against 
random signatures, known signatures, and permuted data/metadata.
}
\description{
While gene signatures are frequently used to predict phenotypes
        (e.g. predict prognosis of cancer patients), it it not always
        clear how optimal or meaningful they are (cf David Venet,
        Jacques E. Dumont, and Vincent Detours' paper "Most Random Gene
        Expression Signatures Are Significantly Associated with Breast
        Cancer Outcome"). Based on suggestions in that paper,
        SigCheck accepts a data set (as an ExpressionSet) and a gene
        signature, and compares its performance on survival and/or 
        classification tasks against 
        a) random gene signatures of the same length; 
        b) known, related and unrelated gene signatures; 
        and c) permuted data and/or metadata.
}
\details{
\tabular{ll}{
Package: \tab SigCheck\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-05-18\cr
License: \tab Artistic-2.0\cr
}
To use \code{SigCheck}, first create anew SigCheck object using the 
function \code{\link{sigCheck}}. 
This will establish the baseline performance
of the signature. 
Next, either run specific checks, or use the high level function 
\code{\link{sigCheckAll}} to run all the core functions in turn.
The three core functions enable 
1) comparison of baseline performance against signatures 
composed of random genes (\code{\link{sigCheckRandom}});
2) comparison of baseline performance against known, and mostly unrelated, 
gene signatures (\code{\link{sigCheckKnown}}); and 
3) comparison of baseline performance against randomly
permuted data and/or metadata (\code{\link{sigCheckPermuted}}).

At a minimum, SigCheck requires a data set (as an \code{\link{ExpressionSet}}),
metadata indicating the membership of each sample in one of two classes,
and a signature (a subset of features in the ExpressionSet). 
If survival data are available, survival analyses are carried out.
Validation samples can be divided into two classes using one of the simple
default methods (based on overall expression value 
or their first principal component).
Alternatively, more sophisticated classification algorithms can be
deployed, using the
\code{\link{MLearn}} function from the \code{MLInterfaces} package 
to build a classifier (using \code{link{smvI}} by default). 
If no validation samples are 
specified, leave-one-out (LOO) cross-validation is utilized to build multiple 
classifiers, each predicting one sample.
If no survival data are provided, signatures are evaluated
based on classification performance.

Output of each check includes the distribution of random performance scores 
(either survival p-value or classification performance)
and the ranking of the passed signature within this distribution. 
An empirical p-value calculation based on 
this rank is also returned indicating confidence that the performance of 
the signature being checked has unique power.
}
\author{
First version written by Justin Norden with Rory Stark at the University of 
Cambridge, Cancer Resaerch UK Cambridge Institute.

Second version, including all survival analysis, written by Rory Stark at CRUK-CI.

Maintainer: Rory Stark <rory.stark@cruk.cam.ac.uk>

}
\references{
Venet, David, Jacques E. Dumont, and Vincent Detours. 
"Most random gene expression signatures are significantly associated with 
breast cancer outcome." PLoS Computational Biology 7.10 (2011): e1002240.
}

\keyword{package}

