# Test Classification Potential of Gene Signature against Random Gene Signatures
# 
# Compares how an input gene signature classifies an expression set object. 
# Compares this classification using an SVM LOO CV analysis to taking random 
# gene signatures of equal length 
# to the input list and testing classification power.

sigCheckRandom <- function(check, iterations=10){
    
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
        doMethod           <- .sigCheckClassifierRandom
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
        doMethod         <- .sigCheckSurvivalRandom
        if(!length(validationSamples)) {
            nosecond <- TRUE
        }
    } else {
        stop("Invalid SigCheck object.")
    }
    
    #expressionSet <- .sigCheckNA(expressionSet)
    
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    featNames <- .sigCheckSignature(expressionSet, signature, annotation,
                                    bReturnFeatures=TRUE)
    featNames <- featNames[!is.na(featNames)]
    
    #Random signature Testing
    sigLength <- length(signature)
    if(bParallel) {
        randomList <- bplapply(1:iterations, doMethod,
                               featNames, sigLength, expressionSet,
                               classes, validationSamples,
                               Method, survival, threshold)
        if(nosecond){
            randomResult  <- unlist(randomList)
        } else {
            randomResult  <- sapply(randomList,function(x){x$result1})
            randomResult2 <- sapply(randomList,function(x){x$result2})
        }
    } else {
        randomResult   <- NULL
        randomResult2 <- NULL
        for (i in 1:iterations){
            checkit <- 
                doMethod(i,featNames, sigLength, expressionSet, 
                         classes, validationSamples, 
                         Method, survival, threshold)
            if(nosecond) {
                randomResult <- c(randomResult, checkit)
            } else {
                randomResult  <- c(randomResult,  checkit$result1)
                randomResult2 <- c(randomResult2, checkit$result2)
            }
        }
    }
    
    output <- NULL
    if(check@checkType=="Classifier") {
        allScoreValues <- c(classifierScore, randomResult)
        sortScores     <- sort(allScoreValues, decreasing=TRUE)
        geneListRank   <- mean(which(sortScores == classifierScore))
        checkPval      <- .sigCheckPval(classifierScore,randomResult,lt=FALSE)
        nullPerf <- 
            .sigCheckClassifierNull(expressionSet, classes, check@validationSamples,
                                    modeVal=modeVal)
        output <- list(checkType="Random",
                       sigPerformance=classifierScore, 
                       modePerformance=nullPerf,
                       checkPval=checkPval,
                       performanceRandom=randomResult)
    } else {
        allScoreValues <- c(survivalPval,randomResult)
        sortScores     <- sort(allScoreValues, decreasing=FALSE)
        geneListRank   <- mean(which(sortScores == survivalPval))
        checkPval      <- .sigCheckPval(survivalPval,randomResult,lt=TRUE)
        output$checkType <- "Random"
        output$survivalPvalsRandom <- randomResult
        output$survivalPval <- survivalPval
        output$checkPval <- checkPval
    }
    
    output$tests <- iterations
    output$rank <- geneListRank
    
    if(check@checkType=="Classifier") {
        if(!nosecond) {
            output$survivalPval <- classifierSurvival
            output$survivalPvalsRandom <- randomResult2
        }
    } else {
        if(!nosecond) {
            output$trainingPvalsRandom <- randomResult2
        }
    }
    
    return(output)
}  

.sigCheckClassifierRandom <- function(i, featNames, sigLength, 
                                      expressionSet, classes,
                                      validationSamples, classifierMethod,
                                      survival=NULL, threshold) {
    expressionSet@signature <- sample(nrow(expressionSet), sigLength)
    classifier <- 
        suppressWarnings(
            sigCheckClassifier(expressionSet, 
                               plotTrainingKM=FALSE, plotValidationKM=FALSE)
        )
    gc()
    if(survival=="") {
        return(result1=classifier@sigPerformance)        
    } else {
        return(list(result1=classifier@sigPerformance,
                    result2=classifier@survivalPval)) 
    }
}

.sigCheckSurvivalRandom <- function(i, featNames, sigLength, 
                                    expressionSet, classes,
                                    validationSamples, survivalMethod,
                                    survival, threshold) {
    expressionSet@signature <- sample(nrow(expressionSet), sigLength)
    check <- 
        suppressWarnings(
            sigCheckSurvival(expressionSet, 
                             plotTrainingKM=FALSE, plotValidationKM=FALSE)
        )
    gc()
    if(validationSamples=="") {
        return(result1=check@survivalPval)        
    } else {
        return(list(result1=check@survivalPval,
                    result2=check@survivalTrainingPval)) 
    }
}

.sigCheckPval = function(performance,scores,lt=T) {
    if(lt) {
        pval <- 1 - (sum(performance < scores)/length(scores))
    } else {
        pval <- 1 - (sum(performance > scores)/length(scores))
    }
    accuracy <- ceiling(log10(length(scores)))
    pval = signif(pval,accuracy)
    return(pval)
}
