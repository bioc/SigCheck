# Test Classification Potential of Gene Signature
# 
# Compares how an input gene signature classifies an expression set object. 
# As default performs a leave
# one out analysis, if a validation set is not provided.

sigCheckClassifier <- function(SigCheck,  
                               plotTrainingKM=TRUE, plotValidationKM=TRUE,
                               impute=FALSE, ...){
    
    classifyFormula <- formula(sprintf("%s~.",SigCheck@classes))
    category <- which(varLabels(SigCheck) %in% SigCheck@classes)
    if(length(unique(SigCheck[[category]])) != 2) {
        stop('ERROR: classes can specify only two categories for classification.')
    }
    if(class(SigCheck[[category]])!="factor") {
        SigCheck[[category]] <- as.factor(SigCheck[[category]])
    }
    
    if (!length(SigCheck@validationSamples)){
        noValidation <- TRUE
        validationSamples <- xvalSpec("LOO")
        #Check to see if validation indices is already an xvalSpec object.
    } else {
        noValidation <- FALSE
        validationSamples <- SigCheck@validationSamples
    }
    
    if(class(validationSamples) == 'xvalSpec'){
        trainingSet <- validationSamples
    } else {
        totalNumSamples <- length(colnames(SigCheck))
        trainingSet <- 1:totalNumSamples
        trainingSet <- trainingSet[-validationSamples]
    }
    
    SigCheck <- .sigCheckNA(SigCheck,impute=impute)
    
    #Subset Whole Matrix Based Off of GeneSigIndices
    
    saveSig   <- SigCheck@signature
    signature <- .sigCheckSignature(SigCheck, SigCheck@signature, 
                                    SigCheck@annotation)
    
    subsetExpressionSet <- SigCheck[signature,]
    
    classifier <- 
        MLearn(classifyFormula, subsetExpressionSet, SigCheck@classifierMethod,
               trainInd=trainingSet,...)  
    classifierConfusionMatrix <- confuMat(classifier)
    classifierScore <- .sigCheckClassifyScore(classifier)
    
    if((SigCheck@survival!="")) {
        if(noValidation) {
            validationSamples=NULL
        }
        classifierSurvival <- 
            .sigCheckClassifierSurvival(expressionSet=SigCheck, 
                                        classes=SigCheck@classes,
                                        signature=signature,
                                        annotation=SigCheck@annotation,
                                        validationSamples=validationSamples,
                                        classifier=classifier,
                                        survival=SigCheck@survival,
                                        threshold=SigCheck@threshold,
                                        plotTrainingKM=plotTrainingKM,
                                        plotValidationKM=plotValidationKM,
                                        survivalLabel=SigCheck@survivalLabel, 
                                        timeLabel=SigCheck@timeLabel)
    }
    
    if(noValidation) {
        nullPerf <- 
            .sigCheckClassifierNull(SigCheck, SigCheck@classes, NULL,
                                    modeVal=SigCheck@modeVal)        
    } else {
        nullPerf <- 
            .sigCheckClassifierNull(SigCheck, SigCheck@classes, 
                                    validationSamples,modeVal=SigCheck@modeVal)
    }  
    
    SigCheck@checkType <- "Classifier"
    SigCheck@sigPerformance <- classifierScore
    SigCheck@confusion <- classifierConfusionMatrix
    SigCheck@modePerformance <- nullPerf
    SigCheck@classifier <- classifier
    SigCheck@signature <- saveSig
    
    if(SigCheck@survival!="") {
        SigCheck@survivalPval <- classifierSurvival$validationPval
        if(!missing(validationSamples)) {
            SigCheck@survivalTrainingPval <- classifierSurvival$trainPval
        }
    }
    return(SigCheck)
}

.sigCheckClassifyScore <- function(classOut) {
    given <- as.character(classOut@testOutcomes)
    predict <- as.character(classOut@testPredictions)
    correct <- sum(predict == given)
    score <- correct / length(given)
    return(score)    
}

.sigCheckClassifierNull <- function(expressionSet, classes, validationSamples,
                                    modeVal){
    if(modeVal=="") {
        modeVal=NULL
    }
    category <- which(varLabels(expressionSet) %in% classes)
    catVals <- unique(expressionSet[[category]])
    if(!length(validationSamples)){
        train <- expressionSet[[category]]
        valid <- expressionSet[[category]]
    } else {
        train <- expressionSet[[category]][-validationSamples]
        valid <- expressionSet[[category]][validationSamples]
    }
    if(is.null(modeVal)) {
        num1 <- sum(train %in% catVals[1])
        num2 <- sum(train %in% catVals[2])
        pick <- which.max(c(num1,num2))
        pickval <- catVals[pick]
    } else {
        if(modeVal %in% catVals) {
            pickval <- modeVal
        } else {
            stop("Mode value must be one of",catVals[1], "or", catVals[2])
        }
    }
    predict <- sum(valid %in% pickval)
    score <- predict/length(valid)
    return(score)
}

.sigCheckSignature <- function(expressionSet, signature, annotation,
                               bReturnFeatures=FALSE) {
    
    if(annotation=="") {
        annotation=NULL
    }
    if(is.null(annotation)) {
        features <- rownames(exprs(expressionSet))
    } else if(length(annotation)>1) {
        if(length(annotation != nrow(expressionSet))) {
            stop("Annotation doesn't have same number of labels as there are features",
                 call.=FALSE)
        } else {
            features <- annotation
        }
    } else {
        anno <- rownames(varMetadata(featureData(expressionSet))) %in% annotation
        if(sum(anno) == 0) {
            stop('Invalid annotation: ',annotation, call.=FALSE)
        }
        features <- featureData(expressionSet)[[which(anno)]]
    }
    if(bReturnFeatures) {
        return(features)
    }
    
    if(class(signature)=="integer") {
        if(max(signature)>length(features)) {
            stop('Invalid signature: feature number greater than number of features.',
                 call.=FALSE)
        }
    } else {
        matches <- signature %in% features
        if(sum(matches) < length(signature)) {
            warning(sprintf("NOTE: %d unmatched features in signature",
                            length(signature)-sum(matches)),call.=FALSE)       
        }
        signature <- which(features %in% signature)
    }
    return(signature)
}

.sigCheckNA <- function(expressionSet, impute=FALSE) {
   if(impute) {
      nas <- apply(exprs(expressionSet),1,function(x)sum(is.na(x)))
      if(sum(nas)) {
         message(sprintf("NOTE: %d values in %d features imputed (%2.4f%%).",
                      sum(nas),sum(nas>0),
                      (sum(nas)/(nrow(expressionSet)*ncol(expressionSet)))*100))      
         exprs(expressionSet) <- t(impute(t(exprs(expressionSet))))
      }
   } else {
      nas <- apply(exprs(expressionSet),1,function(x)sum(is.na(x))>0)
      if(sum(nas)) {
         warning(sprintf("NOTE: %d features with NA values removed",
                         sum(nas)),call.=FALSE)
      }
      expressionSet <- expressionSet[!nas,]
   }
   return(expressionSet)
}


.sigCheckClassifierSurvival <- function(expressionSet, classes, 
                                        signature, annotation, 
                                        validationSamples, 
                                        classifier, 
                                        survival, threshold=median, 
                                        plotTrainingKM, plotValidationKM,
                                        survivalLabel, timeLabel="Time"){
    
    expressionSet <- .sigCheckNA(expressionSet)
    
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    
    if(survivalLabel=="") {
        survivalLabel=survival
    }
    
    classes  <- .sigCheckPheno(expressionSet,classes,bAsFactor=TRUE)
    survival <- .sigCheckPheno(expressionSet,survival,bAsFactor=FALSE)
    
    #Get classification scores    
    savemfrow <- par("mfrow")
    if(length(validationSamples)) {
        par(mfrow=c(1,2))
    }
    
    if(!is.null(validationSamples)) {
        scores    <- trainScores(classifier)     
        trainPval <- .sigCheckPlotClassifierKM(expressionSet,classes,
                                               survival,scores,threshold,
                                               plotKM=plotTrainingKM,
                                               survivalLabel,timeLabel,
                                               main="Survival: Training Set")   
    } else {
        trainPval <- 1
    }
    
    scores <- testScores(classifier)
    validationPval <- .sigCheckPlotClassifierKM(expressionSet,classes,
                                                survival,scores,threshold,
                                                plotKM=plotValidationKM,
                                                survivalLabel,timeLabel,
                                                main="Survival: Validation Set")    
    
    par(mfrow=savemfrow)
    return(list(trainPval=trainPval,validationPval=validationPval))
    
}

.sigCheckPheno <- function(expressionSet,classes,bAsFactor=TRUE) {
    classnum <- which(varLabels(expressionSet) %in% classes)
    classvec <- expressionSet[[classnum]]
    if(bAsFactor) {
        if(class(classvec)!="factor") {
            classvec <- as.factor(classvec)
        }
    }
    return(classvec)
}

.sigCheckPlotClassifierKM <- function(expressionSet,classes,survival,
                                      scores,threshold,
                                      plotKM,survivalLabel,timeLabel,
                                      ...) {
    scorelabels <- rownames(scores)
    scores      <- scores[,1]
    
    validnames <- which(sampleNames(expressionSet) %in% scorelabels)
    classes <- classes[validnames]
    survival <- survival[validnames]
    
    if(is.function(threshold)) {
        threshold <- threshold(scores)    
    } else {
        threshold <- as.numeric(quantile(scores,threshold))
    }
    
    categories <- scores > threshold
    if(categories[1]) {
        l1 <- ">"; l2 <- "<"
    } else {
        l1 <- "<"; l2 <- ">"
    }
    if(length(unique(categories))==1){
        if(categories[1]) {
            warning("All scores above threshold.")
        } else {
            warning("All scores below threshold.")
        }
        return(1)
    }
    categories <- factor(categories,
                         labels=c(sprintf("%s%01.02f",l1,threshold), 
                                  sprintf("%s%01.02f",l2,threshold)))   
    
    pval <- .sigCheckPlotKM(survival, classes, categories, pvalOnly=!plotKM, 
                            xlab=timeLabel, 
                            ylab=survivalLabel, 
                            ...)
    return(pval)
}

