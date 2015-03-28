sigCheck <- function(expressionSet, classes, survival, signature, 
                     annotation, validationSamples, 
                     scoreMethod="PCA1", threshold=median,
                     classifierMethod=svmI, modeVal,
                     survivalLabel="Survival", timeLabel="Time",
                     plotTrainingKM=TRUE, plotValidationKM=TRUE){
    
    if(missing(validationSamples)) {
        validationSamples <- numeric(0)
    }
    if(missing(annotation)) {
        annotation <- ""
    }
    if(missing(survival)) {
        survival <- ""
    }
    if(missing(survivalLabel)) {
        survivalLabel <- survival
    }
    if(missing(timeLabel)) {
        timeLabel <- ""
    }
    if(missing(timeLabel)) {
        timeLabel <- "Time"
    }
    if(missing(modeVal)) {
        modeVal <- ""
    }
    
    result = as(expressionSet,"SigCheckObject")
    result@classes <- classes
    result@annotation <- annotation
    result@survival <- survival
    result@signature <- signature
    result@validationSamples <- validationSamples
    result@threshold <- threshold
    result@classifierMethod <- classifierMethod
    result@survivalLabel <- survivalLabel
    result@timeLabel <- timeLabel
    result@survivalMethod <- scoreMethod
    result@modeVal <- modeVal
    
    if(scoreMethod=="classifier") {
        if(missing(survival)) {
            result@survival <- ""
        } else {
            result@survival <- survival
        }
        
        result@checkType <- "Classifier"
        result <- sigCheckClassifier(result,
                                     plotTrainingKM=plotTrainingKM,
                                     plotValidationKM=plotValidationKM)
    } else if(scoreMethod!="classifier")  {
        if(survival=="") {
            stop("Must specify a survival label.")
        } else {
            result@survival  <- survival
            result@checkType <- "Survival"
        }
        result <- sigCheckSurvival(result,
                                   modeVal=modeVal,
                                   plotTrainingKM=plotTrainingKM,
                                   plotValidationKM=plotValidationKM)       
    } else {
        stop("Invalid scoreMethod: [",scoreMethod,"]")
    }
    
    return(result)
}


sigCheckAll <- function(check,
                        iterations=10, 
                        known="cancer",
                        plotResults=TRUE, ...){
    
    message('Check random signatures...')
    randomGeneOutput  <-
        sigCheckRandom(check,iterations=iterations)
    
    
    message('Check known signatures...')
    knownGenesOutput  <- 
        sigCheckKnown(check,known=known)
    
    toPermute <- "categories"
    if(check@checkType=="Survival" || check@survival!="") {
        nosurv <- FALSE
        toPermute = c(toPermute,"survival")
    } else {
        nosurv <- TRUE
        toPermute = c(toPermute,"features")
    }
    
    
    message(sprintf("Check permuted: [%s]...",toPermute[1]))   
    permute1Output <- 
        sigCheckPermuted(check, toPermute=toPermute[1],
                         iterations=iterations)
    message(sprintf("Check permuted: [%s]...",toPermute[2]))   
    permute2Output <- 
        sigCheckPermuted(check, toPermute=toPermute[2],
                         iterations=iterations) 
    
    if(toPermute[2]=="survival") {
        output <- list(checkRandom=randomGeneOutput,
                       checkKnown=knownGenesOutput,
                       checkPermutedCategories=permute1Output,
                       checkPermutedSurvival=permute2Output)
    } else {
        output <- list(checkRandom=randomGeneOutput,
                       checkKnown=knownGenesOutput,
                       checkPermutedCategories=permute1Output,
                       checkPermutedFeatures=permute2Output)        
    }
    
    if(plotResults) {
        savemfrow = par("mfrow")        
        par(mfrow=c(2,2))
        
        if(check@checkType=="Classifier") {
            sigCheckPlot(output,classifier=TRUE,...)
        }
        
        if(!nosurv) {
            sigCheckPlot(output,classifier=FALSE,...)
        }
        
        par(mfrow=savemfrow)
    }
    
    return(output)
}
