#' A function to calculate and plot ROC curves
#'
#' @param Z.test the bipartite interaction matrix used for the test set
#' @param Z_est the estimated bipartite interaction matrix used
#' @param P the posterior probability matrix output by \code{\link{sample_parameter}}
#' @param plot TRUE/FALSE to plot the ROC curve.
#' @param bins the number of bins that the interval (0,1) is divided into (default is 400)
#' @param all TRUE/FALSE to calculate the ROC curve based on the whole dataset, or only the held-out portion
#' 
#' @description
#' 
#' Computes the ROC curves as x = false positive rate (FPR) and y = true positive rate (TPR)
#' 
#' FPR = False Positives / (False Positives + True Negatives)
#' TPR = True Positives / (True Positives + False Negatives) 
#'
#' @return 
#' Returns: 
#' 'auc': the maximum AUC value
#' 'threshold': the threshold where P > threshold has the maximum AUC value
#' 'roc': a matrix containing the threthold, FPR, and TPR 
#' 
#' @examples
#'\dontrun{}  
#' @export
#'
rocCurves <-
function(Z.test, Z.train, P, plot=TRUE,bins=400, all=FALSE){
    ## Computes the ROC curves as x = FPR and y = TPR
    ##  FPR = FP/(FP +TN)
    ##  TPR = TP/(TP+FN)
    ## Arguments:
    ## plot = TRUE to actually plot the ROC curve
    ## bin  = the number of intervals (0,1) is divided into
    ## all = TRUE, when TPR and FPR are calculated on the whole dataset, FALSE on the held out portion only
    ## Returns:
    ## auc = the maximum AUC value
    ## threshold = the threshold where (P>threshold) has maximum AUC
    ## roc = a matrix of (threshold, FPR, TPR)
    Z.test= 1*(Z.test>0)
    Z.train = 1*(Z.train>0)
    
    u = seq(0, 1, length.out = bins)    # binning
    Z1 = (Z.test - Z.train)==1
    if(all) Z1 = Z.test==1
    m = sum(1*(Z1))
    Z2 = Z.test==0
    n = sum(1*Z2)
    aux = sapply(u, function(r){
                     aux = 1*(P>=r)
                     TP = sum(aux[Z1])            # True positive
                     FN = m - TP                  # False negative
                     FP = sum(aux[Z2])           # False positive
                     TN = n - FP                # True negative
                     c(TP=TP, FN=FN, FP=FP, TN=TN)
                 })
    FPR = aux['FP',]/(aux['FP',] + aux['TN',])
    TPR = aux['TP',]/(aux['TP',] + aux['FN',])
    roc = data.frame(u=u, FPR = FPR, TPR = TPR)
    max.point =   which.min(abs(roc$TPR-1) + roc$FPR)
    threshold = roc$u[max.point]
    n = length(u)
    auc =  0.5*t(abs(FPR[2:n] - FPR[2:n -1]))%*% (TPR[2:n] + TPR[2:n -1])
    auc = round(100*auc,2)
    if(plot){
        plot(FPR, TPR, xlab='FPR (1-specifity)', ylab = 'TPR (sensitivity)',
             type ='l', col='red', lwd=2, main = paste('AUC ', auc),
             xlim = c(0,1), ylim = c(0,1), pch =6, lty=4)
        abline(a = 0, b=1,col='black',lty=2, lwd=2)
    }
    list(auc = auc, threshold = threshold, roc=roc)
}
