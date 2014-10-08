\name{SigCheck-package}
\alias{SigCheck-package}
\alias{SigCheck}
\docType{package}
\title{
Check a gene signature's classification performance against random signatures, 
permuted data, and known signatures.
}
\description{
While gene signatures are frequently used to classify data (e.g. predict 
prognosis of cancer patients), 
it is not always clear how optimal or meaningful they are 
(cf David Venet, Jacques E. Dumont, Vincent Detours' paper 
"Most Random Gene Expression Signatures Are Significantly Associated with 
Breast Cancer Outcome"). 
Based on suggestions in that paper, SigCheck accepts a data set (as an 
ExpressionSet) and a gene signature, and compares its classification performance 
(using the MLInterfaces package) against a) random gene signatures of the same 
length; b) permuted data; and c) known, unrelated gene signatures.
}
\details{
\tabular{ll}{
Package: \tab SigCheck\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-26\cr
License: \tab Artistic-2.0\cr
}
SigCheck provides a high-level function, sigCheck, that runs all the core 
functions in turn. The four core functions enable 1) a genes signature's 
baseline classification performance to be established 
(\code{\link{sigCheckClassifier}}), 2) compares performance against signatures 
composed of random genes (\code{\link{sigCheckRandom}}), 3) 
compares performance against known, and mostly unrelated, gene signatures
(\code{\link{sigCheckKnown}}), and 4) compares performance against randomly
permuted data (\code{\link{sigCheckPermuted}}).

At a minimum, SigCheck requires a data set (as an \code{\link{ExpressionSet}})
and a signature (a subset of features in the ExpressionSet). It uses the 
\code{\link{MLearn}} funciton formt he \code{MLInterfaces} package to build a 
classifier (using \code{link{smvI}} by default) and measure its performance
against validation samples in the ExpressionSet; if no validation samples are 
specified, it uses leave-one-out (LOO) cross-validation to build multiple 
classifiers, each predicting one sample.

Output of each check includes the distribution of random performance scores 
(percentage of validation samples correctly classified) and the ranking of the 
passed signature in this distribution. A simple p-value calculation based on 
this rank is also returned.
}
\author{
Originally written by Justin Norden with Rory Stark at the University of 
Cambridge, Cancer Resaerch UK Cambridge Institute.

Maintainer: Rory Stark <rory.stark@cruk.cam.ac.uk>

}
\references{
Venet, David, Jacques E. Dumont, and Vincent Detours. 
"Most random gene expression signatures are significantly associated with 
breast cancer outcome." PLoS Computational Biology 7.10 (2011): e1002240.
}

\keyword{package}

