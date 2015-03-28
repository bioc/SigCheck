# Plot check checkResults
sigCheckPlot = function(checkResults, classifier=FALSE, ...) {
    if( (length(checkResults)== 4) &&
            (all(unlist(lapply(checkResults,function(x)x$checkType)) == 
                     c("Random","Known",
                       "Permuted", "Permuted")))) {
        savemfrow = par("mfrow")
        par(mfrow=c(2,2))
        for(tocheck in checkResults) {
            sigCheckPlot(tocheck,classifier,...)
        }
        par(mfrow=savemfrow)
    } else {   
        if(classifier) {
            if(is.null(checkResults$sigPerformance)) {
                stop("No classification results present!")
            } else {
                sigCheckPlotClassifier(checkResults, ...)
            }
        } else {
            sigCheckPlotSurvival(checkResults, ...)
        }
    }
}

sigCheckPlotClassifier <- function(checkResults,...){
    
    if( (length(checkResults)== 4) &&
            (names(checkResults) == c("checkClassifier","checkRandom",
                                      "checkKnown","checkPermutedFeatures",
                                      "checkPermutedCategories")) ) {
        for(i in 2:5) {
            sigCheckPlot(checkResults[[i]])
        }
    } else {
        scores <- checkResults$performance
        scores <- sort(scores)
        uscores <- unique(scores)
        numscore <- NULL
        for(uscore in uscores) {
            numscore <- c(numscore,sum(scores %in% uscore))    
        }
        fields <- names(checkResults)
        
        if ("performanceRandom" %in% fields) {
            titlestring <- sprintf("Check: Random Signatures")
        }
        if ("performanceKnown" %in% fields) {
            titlestring <- sprintf("Check: Known Signatures [%s]",
                                   checkResults$known)
        }
        if ("performancePermuted" %in% fields) {
            titlestring <- sprintf("Check: Permuted Data: [%s]",
                                   checkResults$permute)
        }
        xmin <- min(uscores,checkResults$sigPerformance,
                    checkResults$modePerformance)
        xmax <- max(uscores,checkResults$sigPerformance,
                    checkResults$modePerformance)    
        pval <- 1 - (sum(checkResults$sigPerformance > scores)/length(scores))
        accuracy <- ceiling(log10(length(scores)))
        if(pval > 0) {
            rankstr <- sprintf("Percentile:%%.2f (Tests:%%d  p=%%1.%df)",
                               accuracy)
        } else {
            rankstr <- sprintf("Percentile:%%.2f (Tests:%%d  p<%%1.%df)",
                               accuracy)
            pval <- 1/(10^accuracy)
        }
        rankstr <- sprintf(rankstr,ecdf(scores)(checkResults$sigPerformance),
                           length(scores),pval)
        ymax <- max(numscore)
        if(length(unique(numscore))==1){
            ymax <- ymax*2
        }
        plot(uscores,numscore,xlab="Performance",ylab="Frequency",type="b",
             xlim=c(xmin,xmax), ylim=c(0,ymax),
             main=titlestring,sub=rankstr,...)
        abline(v=checkResults$sigPerformance,col="red",lwd=3)
        abline(v=checkResults$modePerformance,col="red",lty="dotted",lwd=2)
        legend("topleft",legend=c(sprintf("Signature [%0.3f]",
                                          checkResults$sigPerformance),
                                  sprintf("Mode [%0.3f]",
                                          checkResults$modePerformance)),
               col="red",lty=c("solid","dotted"),lwd=c(3,2))
    }
}

sigCheckPlotSurvival <- function(checkResults,...){
    
    if( (length(checkResults)== 5) &&
            (names(checkResults) == c("checkClassifier","checkRandom",
                                      "checkKnown","checkPermutedFeatures",
                                      "checkPermutedCategories")) ) {
        for(i in 2:5) {
            sigCheckPlotSurvival(checkResults[[i]])
        }
    } else {
        scores <- checkResults$survivalPvals
        if(is.null(scores)) {
            stop("No survival scores computed.")
        }
        scores <- sort(scores)
        uscores <- unique(scores)
        numscore <- NULL
        for(uscore in uscores) {
            numscore <- c(numscore,sum(scores %in% uscore))    
        }
        fields <- names(checkResults)
        
        if (checkResults$checkType == "Random") {
            titlestring <- sprintf("Survival: Random Signatures")
        }
        if (checkResults$checkType == "Known") {
            titlestring <- sprintf("Survival: Known Signatures [%s]",
                                   checkResults$known)
        }
        if (checkResults$checkType == "Permuted") {
            titlestring <- sprintf("Survival: Permuted Data: [%s]",
                                   checkResults$permute)
        }
        xmin <- min(uscores)
        xmax <- max(uscores)    
        pval <- 1 - (sum(checkResults$survivalPval < scores)/length(scores))
        accuracy <- ceiling(log10(length(scores)))
        
        if(pval > 0) {
            rankstr <- sprintf("Percentile:%%.2f (Tests:%%d  p=%%1.%df)",
                               accuracy)
            pvalstr <- sprintf("Signature p-val =%%1.%df",accuracy)
        } else {
            rankstr <- sprintf("Percentile:%%.2f (Tests:%%d  p<%%1.%df)",
                               accuracy)
            pvalstr <- sprintf("Signature   p-val <%%1.%df",accuracy)        
            pval <- 1/(10^accuracy)
        }
        
        if(checkResults$survivalPval > 0.001) {
            sigpvalstr <- "Signature   p-val =%1.3f (%1.2f)"
            pvalSurvival <- checkResults$survivalPval
        } else {
            sigpvalstr <- "Signature   p-val <%1.3f (%1.2f)"        
            pvalSurvival <- 1/(10^3)
        }
        
        rankstr <- sprintf(rankstr,1-ecdf(scores)(checkResults$survivalPval),
                           length(scores),pval)
        
        plot(density(log10(uscores)),
             #xlim=c(log10(xmin),log10(xmax)),
             main=titlestring,sub=rankstr,...)
        #if(checkResults$survivalPval>xmin) {
        abline(v=log10(checkResults$survivalPval),col="red",lwd=3)
        #}
        
        #if(0.05>xmin) {
        abline(v=log10(.05),col="red",lty="dotted",lwd=2)
        #} 
        legend("topleft",
               legend=c(sprintf(sigpvalstr,pvalSurvival, 
                                log10(checkResults$survivalPval)),
                        "Significant p-val <0.05    (-1.30)"),
               col="red",lty=c("solid","dotted"),lwd=c(3,2))
    }
}

.sigCheckPval = function(performance,scores) {
    pval <- 1 - (sum(performance < scores)/length(scores))
    accuracy <- ceiling(log10(length(scores)))
    pval = signif(pval,accuracy)
    return(pval)
}

