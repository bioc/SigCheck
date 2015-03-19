# Test Classification Potential of Gene Signature
# 
# Compares how an input gene signature classifies an expression set object. 
# As default performs a leave
# one out analysis, if a validation set is not provided.

sigCheckClassifier <- function(expressionSet, classes, signature, 
                               annotation, validationSamples,
                               classifierMethod=svmI,...){
    
    
    classifyFormula <- formula(sprintf("%s~.",classes))
    category <- which(varLabels(expressionSet) %in% classes)
    if(length(unique(expressionSet[[category]])) != 2) {
        stop('ERROR: classes can specify only two categories for classification.')
    }
    if(class(expressionSet[[category]])!="factor") {
        expressionSet[[category]] <- as.factor(expressionSet[[category]])
    }
    
    if (missing(validationSamples)){
        validationSamples <- xvalSpec("LOO")
        #Check to see if validation indices is already an xvalSpec object.
    }
    if(class(validationSamples) == 'xvalSpec'){
        trainingSet <- validationSamples
    } else {
        totalNumSamples <- length(colnames(expressionSet))
        trainingSet <- 1:totalNumSamples
        trainingSet <- trainingSet[-validationSamples]
    } 
    
    expressionSet <- .sigCheckNA(expressionSet)
    
    #Subset Whole Matrix Based Off of GeneSigIndices
    expressionMatrix <- exprs(expressionSet)
    signature <- .sigCheckSignature(expressionSet, signature, annotation)
    
    subsetExpressionSet <- expressionSet[signature,]
    
    classificationGeneList <- 
        MLearn(classifyFormula, subsetExpressionSet, classifierMethod,
               trainInd=trainingSet,...)  
    classifierConfusionMatrix <- confuMat(classificationGeneList)
    classifierScore <- .sigCheckClassifyScore(classificationGeneList)
    nullPerf <- 
        .sigCheckClassifierNull(expressionSet, classes, validationSamples)
    output <- list(sigPerformance=classifierScore,
                   confusion=classifierConfusionMatrix,
                   modePerformance=nullPerf)
    return(output)
}

.sigCheckClassifyScore <- function(classOut) {
    given <- as.character(classOut@testOutcomes)
    predict <- as.character(classOut@testPredictions)
    correct <- sum(predict == given)
    score <- correct / length(given)
    return(score)    
}

.sigCheckClassifierNull <- function(expressionSet, classes, validationSamples){
    if(missing(validationSamples)) {
        validationSamples <- xvalSpec("LOO")
    }
    category <- which(varLabels(expressionSet) %in% classes)
    catVals <- unique(expressionSet[[category]])
    if(class(validationSamples) == 'xvalSpec'){
        train <- expressionSet[[category]]
        valid <- expressionSet[[category]]
    } else {
        train <- expressionSet[[category]][-validationSamples]
        valid <- expressionSet[[category]][validationSamples]
    }
    num1 <- sum(train %in% catVals[1])
    num2 <- sum(train %in% catVals[2])
    pick <- which.max(c(num1,num2))
    pickval <- catVals[pick]
    predict <- sum(valid %in% pickval)
    score <- predict/length(valid)
    return(score)
}

.sigCheckSignature <- function(expressionSet, signature, annotation,
                               bReturnFeatures=FALSE) {
    
    if(missing(annotation)) {
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

.sigCheckNA <- function(expressionSet) {
    nas <- apply(exprs(expressionSet),1,function(x)sum(is.na(x))>0)
    if(sum(nas)) {
        warning(sprintf("NOTE: %d features with NA values removed",
                        sum(nas)),call.=FALSE)
    }
    expressionSet <- expressionSet[!nas,]
    return(expressionSet)
}

