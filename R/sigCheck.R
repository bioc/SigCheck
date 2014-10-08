# Test Classification Potential of Gene Signature as compared to randomly 
# selected gene signatures, permuted expression sets
# and known gene signatures.
# 
# Compares how an input gene signature classifies an expression set object. 
# As default performs a leave
# one out analysis, if a validation set is not provided.

sigCheck <- function(expressionSet, classes, signature, annotation, 
                     validationSamples, 
                     classifierMethod=svmI, nIterations=10, 
                     knownSignatures="cancer", plotResults=TRUE){
    
    expressionSet <- .sigCheckNA(expressionSet)
    
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    
    print('Check classifier for baseline performance...')
    classifierScores <- 
        sigCheckClassifier(expressionSet=expressionSet, classes=classes, 
                           signature=signature,annotation=annotation,
                           classifierMethod=classifierMethod, 
                           validationSamples=validationSamples)
    geneSigPercentCorrect <- classifierScores$sigPerformance
    geneSigConfusionMatrix <- classifierScores$confusion
    
    print('Check random signatures...')
    randomGeneOutput  <-
        sigCheckRandom(expressionSet=expressionSet, classes=classes, 
                       signature=signature,annotation=annotation,
                       classifierMethod=classifierMethod, 
                       validationSamples=validationSamples,
                       nIterations=nIterations,
                       classifierScore=geneSigPercentCorrect)
    print('Check known signatures...')
    knownGenesOutput  <- 
        sigCheckKnown(expressionSet=expressionSet, classes=classes, 
                      signature=signature,annotation=annotation,
                      classifierMethod=classifierMethod, 
                      validationSamples=validationSamples,
                      classifierScore=geneSigPercentCorrect,
                      knownSignatures=knownSignatures)
    print('Check permuted features...')
    permuteRowsOutput <- 
        sigCheckPermuted(expressionSet=expressionSet, classes=classes, 
                         signature=signature, annotation=annotation,
                         classifierMethod=classifierMethod, 
                         validationSamples=validationSamples,
                         nIterations=nIterations,
                         classifierScore=geneSigPercentCorrect,
                         toPermute="features")
    
    print('Check permuted categories...')
    permuteCategoriesOutput <- 
        sigCheckPermuted(expressionSet=expressionSet, classes=classes, 
                         signature=signature, annotation=annotation,
                         classifierMethod=classifierMethod,
                         validationSamples=validationSamples,
                         nIterations=nIterations,
                         classifierScore=geneSigPercentCorrect,
                         toPermute="categories")    
    
    output <- list(checkClassifier=classifierScores,
                   checkRandom=randomGeneOutput,
                   checkKnown=knownGenesOutput,
                   checkPermutedFeatures=permuteRowsOutput,
                   checkPermutedCategories=permuteCategoriesOutput
    )
    
    if(plotResults) {
        par(mfrow=c(2,2))
        sigCheckPlot(output)
    }
    return(output)
}
