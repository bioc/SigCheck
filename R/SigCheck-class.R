
setClass("SigCheckObject",contains = "ExpressionSet",
         slots=c(
             checkType="character",
             classes="character",
             annotation="character",
             survival="character",
             signature="vector",              
             validationSamples="vector",  
             survivalMethod="ANY",
             threshold="ANY", 
             survivalScores="numeric",  
             survivalConfusionMatrix="matrix",
             survivalClassificationScore="numeric",  
             survivalPval="numeric",  
             survivalTrainingScores="numeric",  
             survivalTrainingConfusionMatrix="matrix",
             survivalTrainingClassificationScore="numeric",  
             survivalTrainingPval="numeric",
             survivalLabel="character", 
             timeLabel="character",
             classifierMethod="learnerSchema",
             sigPerformance="numeric", 
             confusion="matrix",
             modeVal="character", 
             modePerformance="numeric", 
             classifier="classifierOutput"
         ))


setMethod("show","SigCheckObject",
          function (object){
              message(class(object)," ")
              message("Type: ",object@checkType," ")
              if(object@checkType=="Classifier") {
                  message("Classifier performance: ",
                          sprintf("%1.2f",object@sigPerformance),
                          ". Confusion matrix: ")
                  print(object@confusion)
              }
              if(length(object@survivalPval)) {
                  message("Survival p-value: ",
                          sprintf("%1.4e",object@survivalPval)," ")
              }
              print(as(object,"ExpressionSet"))
          }
)

setGeneric('checkType', function(object='SigCheckObject') 
    standardGeneric('checkType'))
setMethod('checkType', 'SigCheckObject', 
          function(object) object@checkType)
setGeneric('classes', function(object='SigCheckObject') 
    standardGeneric('classes'))
setMethod('classes', 'SigCheckObject', 
          function(object) object@classes)
setGeneric('annotation', function(object='SigCheckObject') 
    standardGeneric('annotation'))
setMethod('annotation', 'SigCheckObject', 
          function(object) object@annotation)
setGeneric('survival', function(object='SigCheckObject') 
    standardGeneric('survival'))
setMethod('survival', 'SigCheckObject', 
          function(object) object@survival)
setGeneric('signature', function(object='SigCheckObject') 
    standardGeneric('signature'))
setMethod('signature', 'SigCheckObject',
          function(object) object@signature)
setGeneric('validationSamples', function(object='SigCheckObject')
    standardGeneric('validationSamples'))
setMethod('validationSamples', 'SigCheckObject', 
          function(object) object@validationSamples)
setGeneric('survivalMethod', function(object='SigCheckObject') 
    standardGeneric('survivalMethod'))
setMethod('survivalMethod', 'SigCheckObject', 
          function(object) object@survivalMethod)
setGeneric('threshold', function(object='SigCheckObject') 
    standardGeneric('threshold'))
setMethod('threshold', 'SigCheckObject',
          function(object) object@threshold)
setGeneric('survivalType', function(object='SigCheckObject') 
    standardGeneric('survivalType'))
setMethod('survivalType', 'SigCheckObject',
          function(object) object@survivalType)
setGeneric('timeUnits', function(object='SigCheckObject') 
    standardGeneric('timeUnits'))
setMethod('timeUnits', 'SigCheckObject', 
          function(object) object@timeUnits)
setGeneric('survivalScores', function(object='SigCheckObject') 
    standardGeneric('survivalScores'))
setMethod('survivalScores', 'SigCheckObject', 
          function(object) object@survivalScores)
setGeneric('survivalConfusionMatrix', function(object='SigCheckObject')
    standardGeneric('survivalConfusionMatrix'))
setMethod('survivalConfusionMatrix', 'SigCheckObject',
          function(object) object@survivalConfusionMatrix)
setGeneric('survivalClassificationScore', function(object='SigCheckObject')
    standardGeneric('survivalClassificationScore'))
setMethod('survivalClassificationScore', 'SigCheckObject', 
          function(object) object@survivalClassificationScore)
setGeneric('survivalPval', function(object='SigCheckObject') 
    standardGeneric('survivalPval'))
setMethod('survivalPval', 'SigCheckObject', 
          function(object) object@survivalPval)
setGeneric('survivalTrainingScores', function(object='SigCheckObject') 
    standardGeneric('survivalTrainingScores'))
setMethod('survivalTrainingScores', 'SigCheckObject',
          function(object) object@survivalTrainingScores)
setGeneric('survivalTrainingConfusionMatrix', function(object='SigCheckObject')
    standardGeneric('survivalTrainingConfusionMatrix'))
setMethod('survivalTrainingConfusionMatrix', 'SigCheckObject', 
          function(object) object@survivalTrainingConfusionMatrix)
setGeneric('survivalTrainingClassificationScore', 
           function(object='SigCheckObject') 
               standardGeneric('survivalTrainingClassificationScore'))
setMethod('survivalTrainingClassificationScore', 'SigCheckObject',
          function(object) object@survivalTrainingClassificationScore)
setGeneric('survivalTrainingPval', function(object='SigCheckObject')
    standardGeneric('survivalTrainingPval'))
setMethod('survivalTrainingPval', 'SigCheckObject', 
          function(object) object@survivalTrainingPval)
setGeneric('classifierMethod', function(object='SigCheckObject') 
    standardGeneric('classifierMethod'))
setMethod('classifierMethod', 'SigCheckObject', 
          function(object) object@classifierMethod)
setGeneric('sigPerformance', function(object='SigCheckObject')
    standardGeneric('sigPerformance'))
setMethod('sigPerformance', 'SigCheckObject',
          function(object) object@sigPerformance)
setGeneric('confusion', function(object='SigCheckObject') 
    standardGeneric('confusion'))
setMethod('confusion', 'SigCheckObject',
          function(object) object@confusion)
setGeneric('modeVal', function(object='SigCheckObject')
    standardGeneric('modeVal'))
setMethod('modeVal', 'SigCheckObject', 
          function(object) object@modeVal)
setGeneric('modePerformance', function(object='SigCheckObject') 
    standardGeneric('modePerformance'))
setMethod('modePerformance', 'SigCheckObject',
          function(object) object@modePerformance)
setGeneric('classifier', function(object='SigCheckObject') 
    standardGeneric('classifier'))
setMethod('classifier', 'SigCheckObject', 
          function(object) object@classifier)




