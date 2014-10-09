# Test Classification Potential of Gene Signature against Known Gene Signatures
# 
# Compares how an input gene signature classifies an expression set object. 
# Compares this classification using an SVM LOO CV analysis to known gene 
# signatures using the same classification scheme.
# Known gene signatures can either be passed in, or if not the package will 
# test against 48 known signatures from the v'ant Veer paper.

sigCheckKnown <- function(expressionSet, classes, signature, annotation, 
                          validationSamples, 
                          classifierMethod=svmI, classifierScore, 
                          knownSignatures="cancer"){
    
    if(length(knownSignatures)==1) {
        sigName <- knownSignatures
        data(knownSignatures,envir=environment())
        sig <- names(knownSignatures) %in% sigName
        if(sum(sig)==0) {
            stop('Invalid known signature set: ',sigName)
        }
        knownSignatures <- knownSignatures[[which(sig)]]
    } else {
        sigName <- "user specified"
    }
    
    expressionSet <- .sigCheckNA(expressionSet)
    
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    
    #Test the Input Gene List Indicies if classifierScore not passed in.
    if (missing(classifierScore)==TRUE){
        #Get the Classifier Percentage Accuracy from sigCheckClassifier 
        classifierScore <-
            sigCheckClassifier(expressionSet, classes=classes,
                               signature=signature,
                               validationSamples=validationSamples,
                               classifierMethod=classifierMethod)$sigPerformance
    }
    
    #Loop through Known Gene Signatures for Classification
    knownGeneSignatureScoreVector <- c(rep(0,length(knownSignatures)))
    names(knownGeneSignatureScoreVector) <- names(knownSignatures)
    totalSignatureGenesInEsets <- 0
    totalGenesInSignature <- 0
    
    if(missing(validationSamples)) {
        bParallel=FALSE
    } else if(class(validationSamples) == 'xvalSpec') {
        bParallel=FALSE
    } else bParallel=TRUE
    
    if(bParallel) {
        knownGeneSignatureScoreVector <- 
            bplapply(knownSignatures, .sigCheckKnown,
                     expressionSet, classes, validationSamples,
                     classifierMethod, annotation)
        knownGeneSignatureScoreVector <- unlist(knownGeneSignatureScoreVector)
    } else {
        for (i in 1:length(knownSignatures)){
            knownGeneSignatureScoreVector[i] <- 
                .sigCheckKnown(knownSignatures[[i]],
                               expressionSet, classes, validationSamples,
                               classifierMethod, annotation)
        }
    }
    
    
    #  Decide on the significance level of the given gene set by permutation
    allScoreValues <- c(classifierScore, knownGeneSignatureScoreVector)
    sortScores <- sort(allScoreValues, decreasing=TRUE)
    #   Max as if there are multiple matches -- mean instead
    geneListRank <- mean(which(sortScores == classifierScore))
    
    #check matching features
    sigfeatures = NULL
    for(sig in knownSignatures) {
        sigfeatures = unique(c(sigfeatures,sig))
    }
    matches <- sigfeatures %in% .sigCheckSignature(expressionSet, signature, 
                                                annotation,bReturnFeatures=TRUE)
    if(sum(matches) < length(sigfeatures)) {
        warning(
            sprintf("NOTE: %d unmatched features in known signatures (out of %d)",
                        length(sigfeatures)-sum(matches),length(sigfeatures)),
            call.=FALSE)   
    }
    
    #return
    nullPerf <- 
        .sigCheckClassifierNull(expressionSet, classes, validationSamples)
    output <- list(sigPerformance=classifierScore, modePerformance=nullPerf,
                   known=sigName, 
                   knownSigs=length(knownGeneSignatureScoreVector),
                   rank=geneListRank,
                   performanceKnown=knownGeneSignatureScoreVector)
    return(output)
}

.sigCheckKnown <- function(geneSig, expressionSet, classes, validationSamples,
                           classifierMethod, annotation) {
    
    signature <- .sigCheckSignature(expressionSet,toupper(geneSig), annotation)
    #if no overlap
    if(length(signature) == 0){
        return(0)
    }
    tempGeneSigScore <- 
        sigCheckClassifier(expressionSet, classes=classes, signature=signature,
                           validationSamples=validationSamples,
                           classifierMethod=classifierMethod)$sigPerformance
    gc()
    return(tempGeneSigScore)
}


