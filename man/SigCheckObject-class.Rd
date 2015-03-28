\name{SigCheckObject-class}
\Rdversion{1.1}
\docType{class}
\alias{SigCheckObject-class}
\alias{SigCheckObject}
%% \alias{annotation,SigCheckObject-method}
%% \alias{checkType,SigCheckObject-method}
%% \alias{classes,SigCheckObject-method}
%% \alias{classifier,SigCheckObject-method}
%% \alias{classifierMethod,SigCheckObject-method}
%% \alias{confusion,SigCheckObject-method}
%% \alias{modePerformance,SigCheckObject-method}
%% \alias{modeVal,SigCheckObject-method}
%% \alias{signature,SigCheckObject-method}
%% \alias{sigPerformance,SigCheckObject-method}
%% \alias{survival,SigCheckObject-method}
%% \alias{survivalClassificationScore,SigCheckObject-method}
%% \alias{survivalConfusionMatrix,SigCheckObject-method}
%% \alias{survivalMethod,SigCheckObject-method}
%% \alias{survivalPval,SigCheckObject-method}
%% \alias{survivalScores,SigCheckObject-method}
%% \alias{survivalTrainingClassificationScore,SigCheckObject-method}
%% \alias{survivalTrainingConfusionMatrix,SigCheckObject-method}
%% \alias{survivalTrainingPval,SigCheckObject-method}
%% \alias{survivalTrainingScores,SigCheckObject-method}
%% \alias{survivalType,SigCheckObject-method}
%% \alias{threshold,SigCheckObject-method}
%% \alias{timeUnits,SigCheckObject-method}
%% \alias{validationSamples,SigCheckObject-method}

\title{Class \code{"SigCheckObject"}}
\description{
The main object containing everything associated with an expression
dataset and a gene signatures.
Used for subsequent checks of the unique prognostic and/or classification performance of the signature.
Based on an \code{\link{ExpressionSet}} object.
}
\section{Objects from the Class}{
The preferred way to create a \code{\link{SigCheckObject}} is to use the 
\code{\link{sigCheck}} function.
}
\section{Slots}{
  \describe{
    \item{\code{checkType}:}{Object of class \code{"character"} ~~ }
    \item{\code{classes}:}{Object of class \code{"character"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{survival}:}{Object of class \code{"character"} ~~ }
    \item{\code{signature}:}{Object of class \code{"vector"} ~~ }
    \item{\code{validationSamples}:}{Object of class \code{"vector"} ~~ }
    \item{\code{survivalMethod}:}{Object of class \code{"character"} ~~ }
    \item{\code{threshold}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{survivalScores}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{survivalConfusionMatrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{survivalClassificationScore}:}{Object of class \code{"numeric"} ~~ }
    \item{
       \code{survivalPval}:}{Object of class \code{"numeric"}, representing
       performance of the signature in dividing the samples into sets with
       distinct survival prognosis.  
    }
    \item{\code{survivalTrainingScores}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{survivalTrainingConfusionMatrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{survivalTrainingClassificationScore}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{survivalTrainingPval}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{survivalLabel}:}{Object of class \code{"character"} ~~ }
    \item{\code{timeLabel}:}{Object of class \code{"character"} ~~ }
    \item{\code{classifierMethod}:}{Object of class \code{"learnerSchema"} ~~ }
    \item{\code{sigPerformance}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{confusion}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{modeVal}:}{Object of class \code{"character"} ~~ }
    \item{\code{modePerformance}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{classifier}:}{Object of class \code{"classifierOutput"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ExpressionSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "ExpressionSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "ExpressionSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
%%   \describe{
%%     \item{annotation}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{checkType}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{classes}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{classifier}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{classifierMethod}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{confusion}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{modePerformance}{\code{signature(object = "SigCheckObject")}: ... %% }
%%     \item{modeVal}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{signature}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{sigPerformance}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survival}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalClassificationScore}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalConfusionMatrix}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalMethod}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalPval}{
%%        \code{signature(object = "SigCheckObject")}: retrieves value of 
%%        \code{survivalPval} slot.
%%     }
%%     \item{survivalScores}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalTrainingClassificationScore}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalTrainingConfusionMatrix}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalTrainingPval}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalTrainingScores}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{survivalType}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{threshold}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{timeUnits}{\code{signature(object = "SigCheckObject")}: ... }
%%     \item{validationSamples}{\code{signature(object = "SigCheckObject")}: ... }
%% 	 }

No public methods, acess slots directly if required. 
}
%% \references{
%%  ~~put references to the literature/web site here~~
%% }
\author{
Rory Stark
}
\note{
More methods and documentation coming soon...
}



\seealso{
\code{\link{sigCheck}}
}
\examples{
showClass("SigCheckObject")
}
\keyword{classes}
