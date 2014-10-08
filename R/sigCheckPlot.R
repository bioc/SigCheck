# Plot check checkResults
# 
# Plot histogram of checkResults against signature checkResults

sigCheckPlot <- function(checkResults,...){
    
    if( (length(checkResults)== 5) &&
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
        rankstr <- sprintf("Percentile:%.2f (Tests:%d  p=%1.4f)",
                           ecdf(scores)(checkResults$sigPerformance),
                           length(scores),pval)
        plot(uscores,numscore,xlab="Performance",ylab="Frequency",type="b",
             xlim=c(xmin,xmax),
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
