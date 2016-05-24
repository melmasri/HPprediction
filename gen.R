#!/bin/R
##################
tol.err=1e-5
## Gibbs sampling based on sections 2.4 and 2.5
rExp<-function(l,a=1){
    ## Sampling from a truncated exponential distribution
    ##l = outer(c(1,2),c(2,3))
	unif = runif(length(l))
	-log(1- unif*(1-exp(-l*a)))/(l + tol.err)
}

rExp2<-function(l, g, Z, Z0){
    p = 1- exp(-l)
    unif = matrix(runif(length(l)), dim(p))
    ## Method one P(z=0|g) = 1\delta_{s<0} + g\delta_{s>=0}
    U = 1 + 0*p
    U[!Z0] = -log(1- unif[!Z0]*p[!Z0])/(l[!Z0] + tol.err)
    aa = (unif < g*p/(g*p + 1-p)) & Z0
    U[aa] =  -log(1 - unif[aa]*(g*p[aa] + 1-p[aa])/g)/l[aa]
    ## End of Method one
    ## ## U = -log(unif*(1-Z) + (1-unif)*Z )/(l+tol.err)
    ## aa = 1*(runif(sum(Z0))<= ((1-p)/(1-p + g*p))[Z0])
    ## U[Z0]<-    U[Z0]^(1-aa) 
    ## U
    ## Method two P(z=1|g) = g\delta_{s<0} + 1\delta_{s>=0}
    ## U = 1 + 0*p
    ## aa = (unif < p/(p + g*(1-p))) & !Z0
    ## U[aa] = -log(1- unif[aa]*(p[aa] + g*(1-p[aa])))/(l[aa])
    U
}

raffinity<-function(n, m1, m2 , hyper){
    a = hyper[1]
    b = hyper[2]
	rgamma(n,shape = a + m1, rate = b + m2)
}

raffinity.MH<-function(old, eta , m, Ud, sig=0.1, paramEta, hyper){
    a = hyper[1]
    b = hyper[2]
    if(!paramEta) eta = 1
    old= old^(1/eta)
    epsilon = rnorm(length(old), 0, sd = sig) 
    new = old +  -2*(old + epsilon <=0)*epsilon + epsilon    

    likeli = (eta*m+a-1)*log(new/old) -  (new - old)*b - (new^eta - old^eta)*Ud
    
    ratio = pmin(1, exp(likeli))
    u = 1*(runif(length(old))<=ratio)
    u*new^eta + (1-u)*old^eta
}

rHyper<-function(old,y, w, eta, U, pdist, mr, mc, type='host', sig = 0.1){
    new = old*exp(rnorm(length(old), 0, sd = sig ))
    new[2]<-1
    if(type=='host'){
        u = (U*pdist)%*%w
        r = length(y)*(new[1]*log(new[2]) - lgamma(new[1]) )  +
            sum(lgamma(mr+new[1]) - (mr + new[1]*log(new[2] + u))) +
                -(length(y)*(old[1]*log(old[2]) - lgamma(old[1]))  +
                      sum(lgamma(mr+old[1]) - (mr + old[1]*log(old[2] + u))))
    }
    if(type=='parasite'){
        u =  y%*%(U*pdist)
        r = length(w)*(new[1]*log(new[2]) - lgamma(new[1]) )  +
            sum(lgamma(mc+new[1]) - (mc + new[1]*log(new[2] + u))) +
            -(length(w)*(old[1]*log(old[2]) - lgamma(old[1]))  +
                  sum(lgamma(mc+old[1]) - (mc + old[1]*log(old[2] + u))))
    }
    ratio = min(1, exp(r));ratio
    u = 1*(runif(1)<=ratio)
    u*new + (1-u)*old
}

rEta<-function(eta.old, dist , pdist,Z, y, w, U,md, eta_sd =0.01, wEta, yEta){
    ## A function that generates the prior for the power of distance
    eta.prop = eta.old*exp(rnorm(1, 0, sd = eta_sd))
    pdist.old = (dist^eta.old)%*%Z  +  md^eta.old
    if(wEta) w.new= w^(eta.prop/eta.old) else w.new = w
    if(yEta) y.new =y^(eta.prop/eta.old) else y.new = y
    pdist.new = (dist^eta.prop)%*%Z + md^eta.prop
    likli = sum(log((pdist.new/pdist.old)^Z )) - sum(U*(outer(y.new,w.new)*pdist.new - pdist.old*outer(y,w)))
    
    ratio = min(1, exp(likli));ratio
    u = (runif(1)<=ratio)
    if(u) { eta.old  = eta.prop;pdist.old = pdist.new }
    list (eta=eta.old, dist=pdist.old)
}

rg<-function(Z,Y,l){
    ## Method 1 P(Z=0|g) = 1 \delta(s<0) + g \delta(s>=0)
    S = log(l) - 0.5772
    ZZ = 1*(S>=0)
    ## Method 1
    ##ZZ = 1*(Y<1)
    M = sum(ZZ*Z)         #Y >0 , Z=1,  N++ OR U<1(S>0) and Z=1
    N = sum((1-Z)*ZZ) # Y=0, Z =1 , N-+   OR U<1(S>0) and Z=0
    g = rbeta(1 , N + 1, M + 1)
    ## End of method 1
    ## Method 2
    ## ZZ = 1*(S<0)
    ## #ZZ = 1*(Y==1)
    ## M = sum(ZZ*Z)         #Y >0 , Z=1,  N++
    ## N = sum((1-Z)*ZZ) # Y=0, Z =1 , N-+
    ## g = rbeta(1 , M + 1, N+ 1)
    g
}

AdaptiveSigma<-function(param, ls, i, batch.size =50){
    batch = (floor(i/batch.size)-1)*batch.size + 2
    if(!is.null(dim(param))){
        pvar= apply(param[,batch:i],1,sd)
        ac = 1-rowMeans(1*abs(param[,batch:i]- param[,batch:i-1])<=pvar/20)
    }else{
        pvar= sd(param[batch:i])/20
        ac = 1-mean(1*abs(param[batch:i]- param[batch:i-1])<=pvar)
    }
    #print(cbind(i, range(ls), range(ac)))
    ls + sign(ac - 0.44)*(1*(abs(ac - 0.44)>0.03))
}

gibbs_one<-function(Z,y,w,dist, slice = 10, eta,hyper, uncertain =FALSE,wMH=FALSE, yMH=FALSE, wEta = TRUE, yEta =FALSE,updateHyper=TRUE, AdaptiveMC=TRUE){

	## A one step update in a Gibbs sampler.
	## ## initialize
    Z = 1*(Z>0)
    colnames(Z)<-NULL
    rownames(Z)<-NULL
    if(!missing(dist))
        colnames(dist)<-rownames(dist)<-NULL
	n_w = ncol(Z);n_y = nrow(Z)
    n = n_w*n_y
    burn_in = n_w*slice
    throw.out.slice = floor(0.2*slice)
    throw.out = throw.out.slice*ncol(Z)
    print(sprintf("Run for %i slices, and %i burn ins",slice, burn_in))
    print(sprintf("Settings: uncertain=%s, wMH=%s,yMH=%s,wEta=%s,yEta=%s .",uncertain, wMH, yMH,wEta,yEta))
    print(dim(Z))
    ## Hyper parameters
    a_w = 1;b_w = 1; a_y = 1; b_y = 1;

    y_sd = apply(Z, 1,sd)
    w_sd = apply(Z, 2,sd)
    eta_sd = 0.01
    
    beta= 0.05
    
    if(!missing(hyper)){
        print(hyper)
        a_y =  hyper[['host']][1]
        b_y =  hyper[['host']][2]
        y_sd = hyper[['host']][3]
        a_w =  hyper[['parasite']][1]
        b_w =  hyper[['parasite']][2]
        w_sd = hyper[['parasite']][3]
        eta_sd = hyper[['eta']][1]
    }
    
    
    ls = floor(log(tol.err + w_sd)/beta)
    lsy = floor(log(tol.err + y_sd)/beta)
    lseta = floor(log(tol.err + eta_sd)/beta)
    
    batch.size = 50
    
    ## Variable holders
    g0<-rep(0.5, burn_in)
    hh<-matrix(0, nrow=4, ncol=burn_in)
    if(!missing(y)) y0<-matrix( y, nrow= n_y, ncol =burn_in) else
    y0<-matrix(1,nrow= n_y, ncol =burn_in)
    w0<-matrix(1, nrow= n_w, ncol =burn_in)
	Z0 <- Z==0
	mc <- colSums(1*(Z>0))
	mr <- rowSums(1*(Z>0))
    if(!missing(dist))
        min_dist  = min(dist[upper.tri(dist)]) else min_dist=0
    if(!missing(eta)) peta = rep(eta, burn_in)  else
    peta = rep(1, burn_in)
    if(missing(dist)) pdist = 1 else
    pdist = (dist^peta[1])%*%Z
    tryCatch(
        
        for (i in 1:(burn_in-1)){
            ## for (i in 1:200){
            if(AdaptiveMC)
                if((i%%batch.size==0) & (i < 0.8*burn_in)){
                    if(wMH){
                        ls  = AdaptiveSigma(w0, ls, i)
                        w_sd = exp(beta*ls)
                    }
                    if(yMH){
                        lsy = AdaptiveSigma(y0, lsy, i)
                        y_sd = exp(beta*lsy)
                    }
                    if(!missing(eta) & !missing(dist)){
                        lseta = AdaptiveSigma(peta, lseta, i)
                        eta_sd = exp(beta*lseta)
                    }
                    
                }
            
            ## Updatting latent scores
            if(!uncertain){
                U0<- rExp(pdist*outer(y0[,i],w0[,i]))
                U0[Z0]<-1
            }else
                U0 <-rExp2(pdist*outer(y0[,i],w0[,i]), g0[i], Z, Z0)
            
            ## Updating parasite parameters
            if(wMH){
                w0[,i+1]<-raffinity.MH(w0[,i],peta[i],mc,y0[,i]%*%(U0*pdist),
                                       sig=w_sd, wEta, c(a_w, b_w))
            }else
                w0[,i+1]<-raffinity(n_w, mc, y0[,i]%*%(U0*pdist) ,c(a_w, b_w))
            
            ## Updating host parameters
            if(yMH)
                y0[,i+1]<-raffinity.MH(y0[,i],peta[i],mr, (U0*pdist)%*%w0[,i+1],
                                       sig=y_sd, yEta, c(a_y, b_y)) else
            y0[,i+1]<-raffinity(n_y, mr, (U0*pdist)%*%w0[,i+1] ,c(a_y, b_y))

            ## updating eta
            if(!missing(eta)){  
                new.eta = rEta(peta[i],dist,pdist,Z,y0[,i+1],
                    w0[,i+1], U0,min_dist, eta_sd=eta_sd, wEta, yEta)
                peta[i+1] = new.eta$eta
                if(wMH & wEta)
                    w0[,i+1]=w0[,i+1]^(peta[i+1]/peta[i])
                if(yMH & yEta)
                    y0[,i+1]=y0[,i+1]^(peta[i+1]/peta[i])
                pdist = new.eta$dist
                #if(i%%100==0) print(peta[i+1])
            }

            ## Uncertainty variable
            if(uncertain){
                g0[i+1] = rg(Z, U0,l=pdist*outer(y0[,i],w0[,i])) 
            }
            
            ## Updating Hyper parameters
            if(updateHyper){
                if(!yEta ){
                    new = rHyper(c(a_y, b_y),y0[,i+1],w0[,i+1],
                        peta[i+1],U0, pdist, mr, mc, type='host')
                    a_y = new[1]
                    b_y = new[2]
                hh[1,i+1]<-a_y
                    hh[2,i+1]<-b_y
                }
                
                if(!wEta){
                    new = rHyper(c(a_w, b_w),y0[,i+1],w0[,i+1],
                        peta[i+1], U0,pdist, mr, mc, type='parasite')
                    a_w = new[1]
                    b_w = new[2]
                    hh[3,i+1]<-a_w
                    hh[4,i+1]<-b_w
                }
                ##if(i%%100==0) print(c(a_y,b_y,a_w,b_w))
            }
        }
        
       ,warning = function(w)
           {print(c('warning at i:', i)); print(w);traceback()},
        error =function(e)
            {print(c('warning at i:', i)); print(e);traceback()} ,
        finally = print("Done!"))
    ## throwing out the burn_in stage.(approx 30%)
    if(throw.out==0) throw.out = c(1:burn_in) else throw.out = -c(1:throw.out)
    #w0 = w0[,throw.out]
    y0 = if(missing(y)) y0[,throw.out] else y
    w0 = if(missing(w)) w0[,throw.out] else w
    if(!missing(eta)){
        eta = peta[throw.out]
        if(wEta) w0 = w0^(1/eta)
        if(yEta) y0 = y0^(1/eta)
    }else eta = NULL
    g = if(uncertain) t(data.frame(g = g0[throw.out])) else NULL
    param_phy = list(w = w0, y = y0, burn_in = burn_in - max(-throw.out), throw.out = max(-throw.out),eta = eta, g=g, hh=hh, sd = list(w=w_sd, y = y_sd, eta= eta_sd))
    param_phy
}
##################################################
##################################################
## old functions

## rWmc<-function(w.old, y,eta=1, pdist, mc, U, w_sd=0.1, wEta, hyper){
##     a= hyper[1]
##     b = hyper[2]
##     if(!wEta) eta = 1
##     w.old = w.old^(1/eta)
##     epsilon = rnorm(length(w.old), 0, sd = w_sd) # 0.2 for regular 0.12 tested with AdaptiveMCMC
##     w.prop = w.old +  -2*(w.old + epsilon <=0)*epsilon + epsilon
##     likeli = (eta*mc+a-1)*log(w.prop/w.old) - (y%*%(U*pdist))*(w.prop^eta - w.old^eta) - (w.prop - w.old)*b
##                                         #    prior = dgamma(w.prop, shape=a, rate=b, log=TRUE) - dgamma(w.old, shape = a,rate = b, log=TRUE)
##     ratio = pmin(1, exp(likeli))
##     u = 1*(runif(length(w.old))<=ratio)
    
##     u*w.prop^eta + (1-u)*w.old^eta
## }

## rYmc <-function(y.old, w,eta=1, pdist, mr, U, y_sd, yEta, hyper){
##     ##a = a_y# shape parameter of a Gamma prior
##     ##b = b_y# Scale paramter of a Gamma prior
##     a = hyper[1];b=hyper[2]
##     if(!yEta) eta=1
##     y.old = y.old^(1/eta)
##     epsilon = rnorm(length(y.old), 0, sd = y_sd) # 0.2 for regular MH 0.6 tested with Adaptive MCMC
##     y.prop = y.old +  -2*(y.old + epsilon <=0)*epsilon + epsilon
    
##     likeli = (eta*mr+a-1)*log(y.prop/y.old) - (y.prop^eta - y.old^eta)*((U*pdist)%*%w) - (y.prop - y.old)*b
##  #   prior =  dgamma(y.prop, shape=a, rate=b, log=TRUE) - dgamma(y.old, shape = a,rate = b, log=TRUE)
##     ratio = pmin(1, exp(likeli));ratio
##     u = 1*(runif(length(y.old))<=ratio)
##     u*y.prop^eta + (1-u)*y.old^eta
## }
