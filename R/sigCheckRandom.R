# Test Classification Potential of Gene Signature against Random Gene Signatures
# 
# Compares how an input gene signature classifies an expression set object. 
# Compares this classification using an SVM LOO CV analysis to taking random 
# gene signatures of equal length 
# to the input list and testing classification power.

sigCheckRandom <- function(expressionSet, classes, signature, annotation, 
                           validationSamples, 
                           classifierMethod=svmI, nIterations=10, 
                           classifierScore){
    
    expressionSet <- .sigCheckNA(expressionSet)
    
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    featNames <- .sigCheckSignature(expressionSet, signature, annotation,
                                    bReturnFeatures=TRUE)
    featNames <- featNames[!is.na(featNames)]
    
    if(missing(validationSamples)) {
        bParallel=FALSE
    } else if(class(validationSamples) == 'xvalSpec') {
        bParallel=FALSE
    } else bParallel=TRUE
    
    #Test the Input Gene List Indicies if classifierScore not passed in.
    if (missing(classifierScore)==TRUE){
        #Get the Classifier Percentage Accuracy from  sigCheckClassifier 
        classifierScore <- 
            sigCheckClassifier(expressionSet=expressionSet, classes=classes,
                               signature=signature,
                               validationSamples=validationSamples,
                               classifierMethod=classifierMethod)$sigPerformance
    }
    
    #Random signature Testing
    sigLength <- length(signature)
    if(bParallel) {
        randomCorrectList <- bplapply(1:nIterations, .sigCheckClassifier,
                                      featNames, sigLength, expressionSet,
                                      classes, validationSamples,
                                      classifierMethod)
        randomCorrectList <- unlist(randomCorrectList)
    } else {
        randomCorrectList=NULL
        for (i in 1:nIterations){
            result <- 
                .sigCheckClassifier(i,featNames, sigLength, expressionSet, 
                                    classes, validationSamples, 
                                    classifierMethod)
            randomCorrectList <- c(randomCorrectList, result)
        }
    }
    
    #Decide on the significance level of the given gene set by permutation
    allScoreValues <- c(classifierScore, randomCorrectList)
    sortScores <- sort(allScoreValues, decreasing=TRUE)
    # Max as if there are multiple matches -- mean instead
    geneListRank <- mean(which(sortScores == classifierScore))
    nullPerf <- 
        .sigCheckClassifierNull(expressionSet, classes, validationSamples)
    output <- list(sigPerformance=classifierScore, modePerformance=nullPerf,
                   tests=nIterations, rank=geneListRank,
                   performanceRandom=randomCorrectList)
    return(output)
}  

.sigCheckClassifier <- function(i, featNames, sigLength, expressionSet, classes,
                                validationSamples, classifierMethod) {
    randSig <- sample(nrow(expressionSet), sigLength)
    result <- 
        sigCheckClassifier(expressionSet=expressionSet, classes=classes,
                           signature=randSig, 
                           validationSamples=validationSamples,
                           classifierMethod=classifierMethod)$sigPerformance
    gc()
    return(result)
}
