\name{sigCheckAll}
\alias{sigCheckAll}
\title{
Run a default set of checks on a gene signature.
}
\description{
High-level function for package \code{\link{SigCheck}} that runs a default
set of checks against a predictive signature.
}
\usage{
sigCheckAll(check,
            iterations=10, known="cancer", 
            plotResults=TRUE, ...)
}

\arguments{
  \item{check}{
A \code{\link{SigCheckObject}}, as returned by \code{\link{sigCheck}}.
}
  \item{iterations}{
Number of iterations to run to generate background distributions. 
This is how many random signatures the primary signature will be compared to, 
as well as how many of each type of permuted dataset will 
be generated for comparison.
}
  \item{known}{
Specification of a set of known (previously identified) signatures to 
compare to. See \code{\link{sigCheckKnown}} for more details.
}
  \item{plotResults}{
By default, plots of the results will be generated 
unless this is set or \code{FALSE}.
}
  \item{...}{
Extra parameters to pass through \code{\link{sigCheckPlot}}.
}
}

\details{
This high-level function will run four checks, plot the results, and
return a consolidated result set.

First, it calls \code{\link{sigCheckRandom}} to compare the performance 
of \code{interations} randomly selected signatures. 

Next, it calls \code{\link{sigCheckKnown}} to check the 
performance of the signature against a database of signatures previously 
identified to discriminate in other domains.

Finally, two calls are made to \code{\link{sigCheckPermuted}} to check the 
performance of randomly permuted metadata (and possibly expression data).
The first call permutes the category assignments (\code{toPermute="categories"}).
The second call permutes the survival data if they
are available (\code{toPermute="categories"}); 
otherwise it permutes the expression value for each gene (permuting
each row in the \code{\link{ExpressionSet}}, 
equivalent to \code{toPermute="features"}).

If \code{plotResults} is \code{TRUE}, the results are plotted. If a classifier
is involved, a set of four classification results is plotted in
a 2x2 grid, showing how the classification performance of the main
signature compares to that of a mode classifier and to the
distribution of performance values observed for
the random and known signature sets, as well as how it performs using the two 
type of permuted dataset. 
If survival data is available, another 2x2 grid is plotted showing how the 
baseline survival p-value compares to a p-value of 0.05 and to the
distribution of p-values observed for the random and known signatures, as well
as for the permuted data.
}
\value{
A list containing four elements, each containing the result of a check.
\itemize{

\item \code{$checkRandom} is the result list returned by 
\code{\link{sigCheckRandom}}.

\item \code{$checkKnown} is the result list returned by 
\code{\link{sigCheckKnown}}.

\item \code{$checkPermutedCategories} is the result list returned by 
\code{\link{sigCheckPermuted}}
with \code{toPermute="categories"}.

The final element of the result list will be one of the following:

\item \code{$checkPermutedSurvival} is the result list returned by 
\code{\link{sigCheckPermuted}}
with \code{toPermute="survival"}.

\item \code{$checkPermutedFeatures} is the result list returned by 
\code{\link{sigCheckPermuted}}
with \code{toPermute="features"}.
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
\code{\link{sigCheck}}, \code{\link{sigCheckRandom}}, 
\code{\link{sigCheckPermuted}}, \code{\link{sigCheckKnown}}, 
\code{\link{sigCheckPlot}}
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
                  
results <- sigCheckAll(check,iterations=20, 
                       known=knownSignatures$cancer[1:20])                  


## classification analysis
check <- sigCheck(nki, classes="e.dmfs", 
                  signature=knownSignatures$cancer$VANTVEER,
                  annotation="HUGO.gene.symbol",
                  validationSamples=275:319,
                  scoreMethod="classifier")
                  
results <- sigCheckAll(check,iterations=20, 
                       known=knownSignatures$cancer[1:20])  
}

