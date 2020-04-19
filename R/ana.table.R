ana.table <-
function(Z, ZCross, P, roc, plot=FALSE){
    Z = 1*(Z>0)
    ZCross  = 1*(ZCross>0)
    Zpost = 1*(P>roc$threshold)
    data.frame(auc = roc$auc/100, thresh= roc$threshold,
               tot.ones = sum(Z), held.out.ones= sum(abs(ZCross - Z)[Z==1]),
               pred.held.out.ones = 100*sum(Zpost[Z==1 & ZCross==0])/sum(abs(ZCross - Z)[Z==1]),
               pred.tot.ones = sum(Zpost[Z==1])/sum(Z)*100)
    
}
