# Test Classification Potential of Gene Signature against 
#  Expression Set Permutations
# 
# Compares how an input gene signature classifies an expression set object. 
# Compares this classification using an SVM LOO CV analysis to taking random 
# gene signatures of equal length 
# to the input list and testing classification power.

sigCheckPermuted <- function(expressionSet, classes, signature, annotation, 
                             validationSamples, 
                             classifierMethod=svmI, nIterations=10, 
                             classifierScore, toPermute="features"){
    
    expressionSet <- .sigCheckNA(expressionSet)
    
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    
    #Test the Input Gene List Indicies if classifierScore not passed in.
    if (missing(classifierScore)==TRUE){
        #Get the Classifier Percentage Accuracy from sigCheckClassifier 
        classifierScore <- 
            sigCheckClassifier(expressionSet=expressionSet, 
                               classes=classes, 
                               signature=signature, 
                               validationSamples=validationSamples,
                               classifierMethod=classifierMethod)$sigPerformance
    }
    if(missing(validationSamples)) {
        bParallel=FALSE
    } else if(class(validationSamples) == 'xvalSpec') {
        bParallel=FALSE
    } else bParallel=TRUE
    
    #Permute Testing
    permuteScoreList <- c()
    
    if(bParallel) {
        permuteScoreList <- 
            bplapply(1:nIterations, .sigCheckPermuted,
                     expressionSet, toPermute, classes, signature,
                     validationSamples, classifierMethod)
        permuteScoreList <- unlist(permuteScoreList)
    } else {
        for (i in 1:nIterations){            
            res <- .sigCheckPermuted(i, expressionSet, toPermute,
                                     classes, signature,
                                     validationSamples, classifierMethod)        
            permuteScoreList <- c(permuteScoreList, res)
        }
    }
    
    #Compare to Permute Gene List
    #Decide on the significance level of the given gene set by permutation
    allScoreValues <- c(classifierScore, permuteScoreList)
    sortScores <- sort(allScoreValues, decreasing=TRUE)
    # Max as if there are multiple matches -- mean instead
    geneListRank <- mean(which(sortScores == classifierScore))
    nullPerf <- 
        .sigCheckClassifierNull(expressionSet, classes, validationSamples)
    output <- list(sigPerformance=classifierScore, modePerformance=nullPerf,
                   permute=toPermute, tests=nIterations,
                   rank=geneListRank ,performancePermuted=permuteScoreList)
    return(output)
}

.sigCheckPermuted <- function(iter, permuteEset, toPermute, classes, signature,
                              validationSamples, classifierMethod) {    
    
    exMatrix <- exprs(permuteEset)
    #permuteEset <- expressionSet
    
    if ("features" %in% toPermute){
        # Choose Random Rows of the ESET
        mix <- t(apply(exMatrix, 1, sample))
        colnames(mix) <- colnames(exMatrix)
        #Replaces Old Matrix with The new One
        exprs(permuteEset) <- mix
    }
    
    if ("samples" %in% toPermute){
        # Permute the Columns of the ESET
        mix <- apply(exMatrix, 2, sample)
        rownames(mix) <- rownames(exMatrix)
        #Replaces Old Matrix with The new One
        exprs(permuteEset) <- mix
    }
    
    if ("categories" %in% toPermute){
        category <- which(varLabels(permuteEset) %in% classes)
        if(class(validationSamples) == 'xvalSpec') {
            permuteEset[[category]] <- sample(permuteEset[[category]])
        } else {
            permuteEset[[category]][validationSamples] <- 
                sample(permuteEset[[category]][validationSamples])
            permuteEset[[category]][-validationSamples] <- 
                sample(permuteEset[[category]][-validationSamples])
        }
    }
    
    #Testing on Permuted ESET
    permuteScore <- 
        sigCheckClassifier(expressionSet=permuteEset, 
                           classes=classes, signature=signature, 
                           validationSamples=validationSamples,
                           classifierMethod=classifierMethod)$sigPerformance
    gc()
    return(permuteScore)
}
