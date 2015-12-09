# Test Survival Profile  of Gene Signature
sigCheckSurvival <- function(SigCheck,
                             plotTrainingKM=TRUE, plotValidationKM=TRUE, 
                             timeLabel, survivalLabel,impute=FALSE, ...){
   
   category <- which(varLabels(SigCheck) %in% SigCheck@classes)
   if(length(unique(SigCheck[[category]])) != 2) {
      stop('ERROR: classes can specify only two categories for classification.')
   }
   if(class(SigCheck[[category]])!="factor") {
      SigCheck[[category]] <- as.factor(SigCheck[[category]])
   }
   classVals    <- .sigCheckPheno(SigCheck,SigCheck@classes,bAsFactor=TRUE)
   survivalVals <- .sigCheckPheno(SigCheck,SigCheck@survival,bAsFactor=FALSE)
   
   
   
   trainingSet <- 1:length(colnames(SigCheck))
   if(length(SigCheck@validationSamples)){
      trainingSet <- trainingSet[-SigCheck@validationSamples]
   } 
   
   SigCheck <- .sigCheckNA(SigCheck,impute=impute)
   
   #Subset Whole Matrix Based Off of GeneSigIndices
   saveSig <- SigCheck@signature
   SigCheck@signature <- .sigCheckSignature(SigCheck, SigCheck@signature, 
                                            SigCheck@annotation)
   
   
   survivalScores=.sigCheckSurvival(SigCheck, SigCheck@signature,
                                    SigCheck@survival, SigCheck@survivalMethod,
                                    trainingSet,SigCheck@validationSamples)
   
   ## confusion matrices: use classes to compute matrices for training/validation
   trainingPredictions <- 
      .sigCheckSurvivalClassify(survivalScores$trainingScores,
                                classVals[trainingSet],SigCheck@threshold)
   trainingConfusionMatrix <- 
      .sigCheckSurvivalConfuse(trainingPredictions,
                               classVals[trainingSet])
   
   if(length(SigCheck@validationSamples)) {
      validationPredictions <- 
         .sigCheckSurvivalClassify(survivalScores$validationScores,
                                   classVals[SigCheck@validationSamples],
                                   SigCheck@threshold)
      validationConfusionMatrix <-
         .sigCheckSurvivalConfuse(validationPredictions,
                                  classVals[SigCheck@validationSamples])
   }
   
   ## compute scores from matrix for training/validation
   trainingScore <- 
      .sigCheckSurvivalClassificationScore(trainingConfusionMatrix)
   if(length(SigCheck@validationSamples)) {
      validationScore <-  
         .sigCheckSurvivalClassificationScore(validationConfusionMatrix)
   }
   
   ## KM plots
   
   savemfrow=NULL
   if(length(SigCheck@validationSamples)) {
      if(plotTrainingKM & plotValidationKM) {
         savemfrow <- par("mfrow")
         par(mfrow=c(1,2))
      }
   } 
   if(length(SigCheck@validationSamples)) {
      mainstr <- "Survival: Training Set"
   } else {
      mainstr <- "Survival: Full Set"
   }
   categories <- 
      .sigCheckSurvivalCategories(survivalScores$trainingScores,
                                  SigCheck@threshold)    
   if(length(unique(categories))>1){
      trainingPval <- .sigCheckPlotKM(survivalVals[trainingSet],
                                      classVals[trainingSet],
                                      categories,
                                      pvalOnly=!plotTrainingKM, 
                                      xlab=SigCheck@timeLabel, 
                                      ylab=SigCheck@survivalLabel,
                                      main=mainstr)
      
   } else {
      if(plotTrainingKM) {
         warning("No KM plot for training samples: all samples in same category.")
      }
      trainingPval=1
   }
   if(length(SigCheck@validationSamples)){
      categories <- 
         .sigCheckSurvivalCategories(survivalScores$validationScores,
                                     SigCheck@threshold)
      if(length(unique(categories))>1){
         validationPval <- .sigCheckPlotKM(survivalVals[SigCheck@validationSamples],
                                           classVals[SigCheck@validationSamples],
                                           categories,
                                           pvalOnly=!plotValidationKM, 
                                           xlab=SigCheck@timeLabel, 
                                           ylab=SigCheck@survivalLabel,
                                           main="Survival: Validation Set") 
      } else {
         if(plotValidationKM) {
            warning("No KM plot for validation samples: all samples in same category.")
         }
         validationPval=1
      }
   } else validationPval=1
   
   if(!is.null(savemfrow)) {
      par(mfrow=savemfrow)
   }
   
   ## perfect survival: survivors vs. non-survivors
   
   SigCheck@checkType <- "Survival"
   
   if(length(SigCheck@validationSamples)) {
      SigCheck@survivalScores              <- survivalScores$validationScores
      SigCheck@survivalConfusionMatrix     <- validationConfusionMatrix
      SigCheck@survivalClassificationScore <- validationScore
      SigCheck@survivalPval                <- validationPval
      SigCheck@survivalTrainingScores              <- survivalScores$trainingScores
      SigCheck@survivalTrainingConfusionMatrix     <- trainingConfusionMatrix
      SigCheck@survivalTrainingClassificationScore <- trainingScore
      SigCheck@survivalTrainingPval                <- trainingPval
   } else {
      SigCheck@survivalScores              <- survivalScores$trainingScores
      SigCheck@survivalConfusionMatrix     <- trainingConfusionMatrix
      SigCheck@survivalClassificationScore <- trainingScore
      SigCheck@survivalPval                <- trainingPval
   }
   
   return(SigCheck)
}

.sigCheckSurvival <- function(expressionSet, signature,
                              survival,survivalMethod,
                              trainingSet,validationSamples) {
   
   survivalValues <- which(varLabels(expressionSet) %in% survival)
   subsetExpressionSet <- expressionSet[signature,trainingSet]
   if(!is.character(survivalMethod)){
      trainingScores <- 
         survivalMethod(subsetExpressionSet)
      if(length(validationSamples)) {
         subsetExpressionSet <- expressionSet[signature,validationSamples]
         validationScores <- 
            survivalMethod(subsetExpressionSet)
      } else {
         validationScores <- trainingScores
      }  
   } else {
      if(survivalMethod=="PCA1") {
         trainingScores <- 
            .sigCheckPCA1(subsetExpressionSet)
         if(length(validationSamples)) {
            subsetExpressionSet <- expressionSet[signature,validationSamples]
            validationScores <- 
               .sigCheckPCA1(subsetExpressionSet)
         } else {
            validationScores <- trainingScores
         }
      } else if (survivalMethod=="High") {
         trainingScores <- 
            .sigCheckHigh(subsetExpressionSet)
         if(length(validationSamples)) {
            subsetExpressionSet <- expressionSet[signature,validationSamples]
            validationScores <- 
               .sigCheckHigh(subsetExpressionSet)
         } else {
            validationScores <- trainingScores
         }    
      } else {
         stop("Invalid scoreMethod: ",survivalMethod)  
      }
   }
   return(list(trainingScores=trainingScores,
               validationScores=validationScores))
}

.sigCheckPCA1 <- function(expressionSet){
   pcs <- prcomp(t(exprs(expressionSet)))
   pc1 <- pcs$x[,1]
   return(pc1)
}


.sigCheckHigh <- function(expressionSet,fn=mean){
   means <- apply(exprs(expressionSet),2,fn)
   return(means)
}

.sigCheckSurvivalClassify <- function(scores,classVals,threshold,bFlip=TRUE){
   categories <- .sigCheckSurvivalCategories(scores,threshold)
   predict <- scores == "High"
   if(bFlip) {
      actual <- classVals == levels(classVals)[2]
      correct <- predict == actual
      if(sum(correct) > sum(!correct)) {
         return(predict)
      } else {
         return(!predict)
      }
   } else {
      return(predict)
   }
}

.sigCheckSurvivalConfuse <- function(predictions,classVals) {
   vals <- unique(classVals)
   wrong <- rep(vals[1],length(classVals))
   wrong[classVals==vals[1]] <- vals[2]
   predicted <- classVals
   predicted[!predictions] <- wrong[!predictions]
   tn <- sum(classVals==vals[1] & classVals==predicted)
   fn <- sum(classVals==vals[2] & classVals!=predicted)
   tp <- sum(classVals==vals[2] & classVals==predicted)
   fp <- sum(classVals==vals[1] & classVals!=predicted)
   
   confu <- matrix(c(tn,fn,fp,tp),2,2)
   rownames(confu) <- 0:1
   colnames(confu) <- 0:1
   return(confu)
}

.sigCheckSurvivalClassificationScore <- function(confusion) {
   total <- sum(confusion)
   correct <- confusion[1,1] + confusion[2,2]
   score <- correct/total
   return(score)
}

.sigCheckPlotKM <- function(y, status, x,
                            barcol=c("blue", "red", "grey"),
                            pvalOnly=FALSE,twogroups=TRUE, ...) {
   status <- as.numeric(status)-1
   
   if(twogroups) {
      mid <- x == "Mid"
   } else {
      mid <- rep(FALSE,length(x))
   }
   xx = x[!mid]
   yy = y[!mid]
   ss = status[!mid]
   
   
   pval <- pchisq(survdiff(Surv(yy, ss) ~xx)$chisq, 
                  length(unique(xx[!is.na(xx)]))-1, lower.tail=FALSE)
   
   if(!pvalOnly) {
      
      plot(survfit(Surv(y, status) ~ x), col=barcol,...)
      
      if (pval < 0.001) {
         pvalStr <- "<0.001"
      } else {
         pvalStr <- paste("=", round(pval, 3))
      }
      text(x=max(y)/12, y=0.1, paste("p", pvalStr))
      p1 <- table(x[!is.na(status)])
      p2 <- tapply(status, x, sum, na.rm=TRUE)
      p2[is.na(p2)] <- 0
      text.leg <- paste(levels(x) , ": ", p1, "(", p2, ")", sep="")
      legend("topright", lty=rep(1,3), lwd=2, col=barcol, 
             legend=text.leg, bty="n")
   }
   return(pval)
}


.sigCheckSurvivalCategories <- function(scores,threshold) {
   if(is.function(threshold)) {
      threshold <- threshold(scores)    
   } else {
      threshold = sort(threshold)
      threshold <- as.numeric(quantile(scores,threshold))
   }
   two <- rep(TRUE,length(scores))
   if(length(threshold)==1) {
      predict <- scores > threshold
      categories <- rep("Low",length(predict))
      categories[predict] <- "High"
   } else if(length(threshold)==2){
      low  <- scores <  threshold[1]
      high <- scores >= threshold[2]
      mid  <- !(low | high)
      categories <- rep("Low",length(scores))
      categories[high] <- "High"
      categories[mid] <- "Mid"
      two[mid] <- FALSE
   } else {
      stop("At most three categories supported")
      categories=rep(length(threshold)+1,length(scores))
      for(i in length(threshold):1) {
         categories[scores<threshold[i]]=i   
      }
   }
   categories <- factor(categories)   
   return(categories)
}

