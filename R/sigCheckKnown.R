# Test Classification Potential of Gene Signature against Known Gene Signatures
# 
# Compares how an input gene signature classifies an expression set object. 
# Compares this classification using an SVM LOO CV analysis to known gene 
# signatures using the same classification scheme.
# Known gene signatures can either be passed in, or if not the package will 
# test against 48 known signatures from the v'ant Veer paper.

sigCheckKnown <- function(check, known="cancer"){
    
    bParallel <- TRUE
    
    
    expressionSet     <- check
    classes           <- check@classes
    signature         <- check@signature
    annotation        <- check@annotation 
    validationSamples <- check@validationSamples 
    survival          <- check@survival
    threshold         <- check@threshold
    nosecond          <- FALSE
    
    if(check@checkType=="Classifier") {
        classifier         <- check
        modeVal            <- check@modeVal
        Method             <- check@classifierMethod
        classifierScore    <- check@sigPerformance
        classifierSurvival <- check@survivalPval
        doMethod           <- .sigCheckKnownClassifier
        if(!length(validationSamples)) {
            bParallel <- FALSE
            nosecond  <- TRUE
            validationSamples <- xvalSpec("LOO")
        } 
        if(survival=="") {
            nosecond <- TRUE
        }
    } else if(check@checkType == "Survival") {
        Method           <- check@survivalMethod
        survivalLabel     <- check@survivalLabel
        timeLabel        <- check@timeLabel
        survivalPval     <- check@survivalPval
        doMethod         <- .sigCheckKnownSurvival
        if(!length(validationSamples)) {
            nosecond <- TRUE
        }
    } else {
        stop("Invalid SigCheck object.")
    }
    
    if(length(known)==1) {
        sigName <- known
        data(knownSignatures,envir=environment())
        sig <- names(knownSignatures) %in% sigName
        if(sum(sig)==0) {
            stop('Invalid known signature set: ',sigName)
        }
        known <- knownSignatures[[which(sig)]]
    } else {
        sigName <- "user specified"
    }
    
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    
    #Loop through Known Gene Signatures for Classification
    knownResult1 <- c(rep(0,length(known)))
    names(knownResult1) <- names(known)
    knownResult2 <- knownResult1
    totalSignatureGenesInEsets <- 0
    totalGenesInSignature <- 0
    
    if(bParallel) {
        knownList <- 
            bplapply(known, doMethod,
                     expressionSet, classes, validationSamples,
                     Method, annotation,
                     survival, threshold)
        if(nosecond){
            knownResult1 <- unlist(knownList)
        } else {
            knownResult1 <- sapply(knownList,function(x){x$result1})
            knownResult2 <- sapply(knownList,function(x){x$result2})
        }
    } else {
        knownResult1  <- NULL
        knownResult2  <- NULL
        for (i in 1:length(known)){
            checkit <- 
                doMethod(known[[i]],
                         expressionSet, classes, validationSamples,
                         Method, annotation,
                         survival, threshold)
            if(nosecond) {
                knownResult1 <- c(knownResult1, checkit)
            } else {
                knownResult1 <- c(knownResult1, checkit$result1)
                knownResult2 <- c(knownResult2, checkit$result2)
            }
        }
    }
    if(nosecond){
        names(knownResult1) = names(known)
    }
    output <- NULL
    if(check@checkType=="Classifier") {
        allScoreValues <- c(classifierScore, knownResult1)
        sortScores     <- sort(allScoreValues, decreasing=TRUE)
        geneListRank   <- mean(which(sortScores == classifierScore))
        checkPval      <- .sigCheckPval(classifierScore,knownResult1,lt=FALSE)
        nullPerf <- 
            .sigCheckClassifierNull(expressionSet, classes, check@validationSamples,
                                    modeVal=modeVal)
        output <- list(checkType="Known",
                       sigPerformance=classifierScore, 
                       modePerformance=nullPerf,
                       checkPval=checkPval,
                       performanceKnown=knownResult1)
    } else {
        allScoreValues       <- c(survivalPval,knownResult1)
        sortScores           <- sort(allScoreValues, decreasing=FALSE)
        geneListRank         <- mean(which(sortScores == survivalPval))
        geneListRank         <- mean(which(sortScores == survivalPval))
        checkPval            <- .sigCheckPval(survivalPval,knownResult1,lt=TRUE)
        output$checkType <- "Known"
        output$survivalPvalsKnown <- knownResult1
        output$survivalPval  <- survivalPval
        output$checkPval <- checkPval
    }
    
    #check matching features
    sigfeatures = NULL
    for(sig in known) {
        sigfeatures = unique(c(sigfeatures,sig))
    }
    matches <- sigfeatures %in% .sigCheckSignature(expressionSet, signature, 
                                                   annotation,
                                                   bReturnFeatures=TRUE)
    if(sum(matches) < length(sigfeatures)) {
        warning(
            sprintf("NOTE: %d unmatched features in known signatures (out of %d)",
                    length(sigfeatures)-sum(matches),length(sigfeatures)),
            call.=FALSE)   
    }
    
    #return
    output$known     <- sigName 
    output$knownSigs <- length(knownResult1)
    output$rank      <- geneListRank
    
    if(check@checkType=="Classifier") {
        if(!nosecond) {
            output$survivalPval  <- classifierSurvival
            output$survivalPvalsKnown <- knownResult2
        }
    } else {
        if(!nosecond) {
            output$trainingPvalsKnown <- knownResult2
        }
    }
    
    return(output)
}

.sigCheckKnownClassifier <- function(geneSig, expressionSet, 
                                     classes, validationSamples,
                                     Method, annotation, survival, threshold) {
    
    signature <- .sigCheckSignature(expressionSet,toupper(geneSig), annotation)
    #if no overlap
    if(length(signature) == 0){
        return(list(result1=NA,result2=NA)) 
    }
    expressionSet@signature <- signature
    classifier <- 
        sigCheckClassifier(expressionSet, 
                           plotTrainingKM=FALSE, plotValidationKM=FALSE)
    gc()
    if(survival=="") {
        return(result1=classifier@sigPerformance)        
    } else {
        return(list(result1=classifier@sigPerformance,
                    result2=classifier@survivalPval)) 
    }
    
}

.sigCheckKnownSurvival <- function(geneSig, expressionSet, 
                                   classes, validationSamples,
                                   Method, annotation, survival, threshold) {
    
    signature <- .sigCheckSignature(expressionSet,toupper(geneSig), annotation)
    #if no overlap
    if(length(signature) == 0){
        return(list(result1=NA,result2=NA)) 
    }
    expressionSet@signature <- signature
    check <- 
        suppressWarnings(
            sigCheckSurvival(expressionSet,
                             plotTrainingKM=FALSE, plotValidationKM=FALSE)
        )
    gc()
    if(is.null(validationSamples)) {
        return(result1=check@survivalPval)        
    } else {
        return(list(result1=check@survivalPval,
                    result2=check@survivalTrainingPval)) 
    }
    
}


