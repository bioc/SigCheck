# Test Classification Potential of Gene Signature against 
#  Expression Set Permutations
# 
# Compares how an input gene signature classifies an expression set object. 
# Compares this classification using an SVM LOO CV analysis to taking random 
# gene signatures of equal length 
# to the input list and testing classification power.

sigCheckPermuted <- function(check, toPermute="categories", iterations=10){
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
        if(!length(validationSamples)) {
            nosecond <- TRUE
        }
    } else {
        stop("Invalid SigCheck object.")
    }  
    
    #expressionSet <- .sigCheckNA(expressionSet)
    
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    
    #Permute Testing
    permuteResult1 <- c()
    permuteResult2 <- c()
    if(bParallel) {
        checkit <- 
            bplapply(1:iterations, .sigCheckPermuted, 
                     expressionSet, check@checkType, 
                     toPermute, classes, signature,
                     validationSamples, Method,
                     survival, threshold, nosecond)
        if(nosecond){
            permuteResult1  <- unlist(checkit)
        } else {
            permuteResult1 <- sapply(checkit,function(x){x$result1})
            permuteResult2 <- sapply(checkit,function(x){x$result2})
        }
    } else {
        permuteResult1  <- NULL
        permuteResult2 <- NULL
        for (i in 1:iterations){    
            checkit <- .sigCheckPermuted(i, expressionSet, check@checkType,
                                         toPermute, classes, signature,
                                         validationSamples, Method,
                                         survival, threshold, nosecond)        
            if(nosecond) {
                permuteResult1 <- c(permuteResult1, checkit)
            } else {
                permuteResult1 <- c(permuteResult1, checkit$result1)
                permuteResult2 <- c(permuteResult2, checkit$result2)
            }
        }
    }
    
    output <- NULL
    if(check@checkType=="Classifier") {
        allScoreValues <- c(classifierScore, permuteResult1)
        sortScores     <- sort(allScoreValues, decreasing=TRUE)
        geneListRank   <- mean(which(sortScores == classifierScore))
        checkPval      <- .sigCheckPval(classifierScore,permuteResult1,lt=FALSE)
        nullPerf <- 
            .sigCheckClassifierNull(expressionSet, classes, check@validationSamples,
                                    modeVal=modeVal)
        output <- list(checkType="Permuted",
                       sigPerformance=classifierScore, 
                       modePerformance=nullPerf,
                       checkPval=checkPval,
                       performancePermuted=permuteResult1)
    } else {
        allScoreValues <- c(survivalPval,permuteResult1)
        sortScores     <- sort(allScoreValues, decreasing=FALSE)
        geneListRank   <- mean(which(sortScores == survivalPval))
        checkPval      <- .sigCheckPval(survivalPval,permuteResult1,lt=TRUE)
        output$checkType             <- "Permuted"
        output$survivalPvalsPermuted <- permuteResult1
        output$survivalPval          <- survivalPval
        output$checkPval             <- checkPval
    }
    
    output$permute <- toPermute
    output$tests   <- iterations
    output$rank   <- geneListRank
    
    if(check@checkType=="Classifier") {
        if(!nosecond) {
            output$survivalPval <- classifierSurvival
            output$survivalPvalsPermuted <- permuteResult2
        }
    } else {
        if(!nosecond) {
            output$trainingPvalsPermuted <- permuteResult2
        }
    }
    
    return(output)
    
}

.sigCheckPermuted <- function(iter, permuteEset, checkType,
                              toPermute, classes, signature,
                              validationSamples, Method,
                              survival=NULL, threshold, nosecond) {    
    
    exMatrix <- exprs(permuteEset)
    #permuteEset <- expressionSet
    
    if ("features" %in% toPermute){
        mix <- t(apply(exMatrix, 1, sample))
        colnames(mix) <- colnames(exMatrix)
        exprs(permuteEset) <- mix
    }
    
    if ("samples" %in% toPermute){
        mix <- apply(exMatrix, 2, sample)
        rownames(mix) <- rownames(exMatrix)
        exprs(permuteEset) <- mix
    }
    
    if ("categories" %in% toPermute){
        category <- which(varLabels(permuteEset) %in% classes)
        if(is.null(validationSamples)) {
            permuteEset[[category]] <- sample(permuteEset[[category]])
        } else {
            permuteEset[[category]][validationSamples] <- 
                sample(permuteEset[[category]][validationSamples])
            permuteEset[[category]][-validationSamples] <- 
                sample(permuteEset[[category]][-validationSamples])
        }
    }
    
    if ("survival" %in% toPermute){
        if(is.null(survival)) {
            stop("No survival times specified.")
        }
        category <- which(varLabels(permuteEset) %in% classes)
        if(is.null(validationSamples)) {
            permuteEset[[category]] <- sample(permuteEset[[category]])
        } else {
            permuteEset[[category]][validationSamples] <- 
                sample(permuteEset[[category]][validationSamples])
            permuteEset[[category]][-validationSamples] <- 
                sample(permuteEset[[category]][-validationSamples])
        }
    }
    permuteEset@signature <- signature
    if(checkType=="Classifier") {
        classifier <- 
            sigCheckClassifier(permuteEset, 
                               plotTrainingKM=FALSE, plotValidationKM=FALSE)
        gc()
        if(nosecond) {
            return(result1=classifier@sigPerformance)        
        } else {
            return(list(result1=classifier@sigPerformance,
                        result2=classifier@survivalPval)) 
        }
    } else {
        check <- 
            suppressWarnings(
                sigCheckSurvival(permuteEset,
                                 plotTrainingKM=FALSE, plotValidationKM=FALSE)
            )
        gc()
        if(nosecond) {
            return(result1=check@survivalPval)        
        } else {
            return(list(result1=check@survivalPval,
                        result2=check@survivalTrainingPval)) 
        }   
    }    
}
