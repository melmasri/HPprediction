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
    unif = runif(length(l))
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

raffinity.MH<-function(old, z, Ud, sig=0.1, hyper){
    a = hyper[1]
    b = hyper[2]
    epsilon = rnorm(length(old), 0, sd = sig)
    new = old +  -2*(old + epsilon <=0)*epsilon + epsilon    

    likeli = (z+a-1)*log(new/old) -  (new - old)*(b +Ud)
    
    ratio = pmin(1, exp(likeli))
    u = 1*(runif(length(old))<=ratio)
    u*new + (1-u)*old
}

rHyper.Diag<-function(old,Psi,m, sig = 0.1){
    new = old*exp(rnorm(length(old), 0, sd = sig))
    new[2]<-1
    r = length(old)*(new[1]*log(new[2]) - lgamma(new[1]) )  +
        sum(lgamma(m+new[1]) - (m + new[1])*log(new[2] + Psi)) +
            -(length(old)*(old[1]*log(old[2]) - lgamma(old[1]))  +
                  sum(lgamma(m+old[1]) - (m + old[1])*log(old[2] + Psi)))
    ratio = min(1, exp(r));ratio
    u = 1*(runif(1)<=ratio)
    u*new + (1-u)*old
}

rHyper.Horiz<-function(old,Psi, m, sig = 0.1){
    new = old*exp(rnorm(length(old), 0, sd = sig ))
    new[2]<-1
    r = sum(new[1]*log(new[2]) - lgamma(new[1]) + 
                lgamma(new[1]+m) - (new[1]+m)*log(new[2]+Psi)) -
                    sum(old[1]*log(old[2]) - lgamma(old[1]) + 
                            lgamma(old[1]+ m ) - (old[1]+m)*log(old[2] + Psi))
    ratio = min(1, exp(r));ratio
    u = 1*(runif(1)<=ratio)
    u*new + (1-u)*old
}

rEta.Diag<-function(eta.old,dist,pdist.old,zdiag, Z, y, w, U, eta_sd =0.01){
    ## A function that generates the prior for the power of distance
    eta.prop = eta.old*exp(rnorm(1, 0, sd = eta_sd))
    dist1=0*dist
    dist1[upper.tri(dist1)] =dist[upper.tri(dist)]^eta.prop
    dist1 = dist1+ t(dist1)
    pdist.new = dist1%*%Z +tol.err
    likeli = sum(log((pdist.new[zdiag]/pdist.old[zdiag])^Z[zdiag] ))-
        sum(U*y*w*(pdist.new[zdiag] - pdist.old[zdiag]))
    
    ratio = min(1, exp(likeli));ratio
    u = (runif(1)<=ratio)
    if(u) {eta.old  = eta.prop;pdist.old = pdist.new}
    list (eta=eta.old, dist=pdist.old)
}

rEta.Horiz<-function(eta.old,dist,pdist.old,i, Z, y, w, U, eta_sd =0.01){
    ## A function that generates the prior for the power of distance
    eta.prop = eta.old*exp(rnorm(1, 0, sd = eta_sd))
    pdist.new = (dist[i,]^eta.prop) %*% Z + tol.err
    likeli = sum(log((pdist.new/pdist.old)^Z[i,] ))-
        sum(U*y*w*(pdist.new - pdist.old))
    
    ratio = min(1, exp(likeli));ratio
    u = (runif(1)<=ratio)
    if(u) { eta.old  = eta.prop; pdist.old = c(pdist.new)}
    list (eta=eta.old, dist=pdist.old)
}

rg<-function(Z,Y,l){
    ## g0[i+1] = rg(Z, U0,l=pdist*outer(y0[,i],w0[,i])) 
    ## Method 1 P(Z=0|g) = 1 \delta(s<0) + g \delta(s>=0)
    S = log(l) #- 0.5772
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

ICM_est<-function(Z,dist, slice = 10, eta,hyper, uncertain =FALSE,updateHyper=FALSE, AdaptiveMC=FALSE, distOnly=FALSE, ICM.HORIZ=TRUE){

    ## Warnings
    if(missing(Z)) stop('Interaction matrix is missing!')

    ## Setting up parameters
	n_w = ncol(Z);n_y = nrow(Z)
    n = n_w*n_y
    burn_in = slice
    throw.out = floor(0.2*slice)
    print(sprintf("Run for %i slices, and %i burn ins",slice, throw.out))
    print(sprintf("Adaptive MC: %s",AdaptiveMC))
    print(sprintf("Updating Hyperparameters: %s",updateHyper))
    print(sprintf("Using uncertainity: %s",uncertain))
    print(sprintf("Using affinity parameters: %s",!distOnly))
    print(sprintf("Using Similarity matrix: %s",!missing(dist)))
    print(sprintf("Updating eta: %s",!missing(eta)))
    print(sprintf("ICM Horizontal(TRUE), Diagonal(FALSE): %s.",ICM.HORIZ))
    print(sprintf('Input matrix is Binary: %s, with dimension:', all(range(Z)==c(0,1))))
    print(dim(Z))
    
    ## Hyper parameters
    a_w = 1;b_w = 1; a_y = 1; b_y = 1;
    y=1;w=1;etaStart=1
    y_sd = apply(Z, 1,sd)
    w_sd = apply(Z, 2,sd)
    eta_sd = 0.01
    beta= 0.05
    
    if(!missing(hyper)){
        ## Setting up host and parasite hyper parameters
        if(!is.null(hyper[['hostHyper']])){
            a_y =  hyper[['hostHyper']][1]
            b_y =  hyper[['hostHyper']][2]
            y_sd = 0.1
        }
        if(!is.null(hyper['parasiteHyper'])){
            a_w =  hyper[['parasiteHyper']][1]
            b_w =  hyper[['parasiteHyper']][2]
            w_sd =0.1
        }
        ## Setting up starting values
        if(!is.null(hyper[['hostStart']])) y = hyper[['hostStart']] else y = 1
        if(!is.null(hyper[['parasiteStart']])) w = hyper[['parasiteStart']] else w = 1
        if(!is.null(hyper[['etaStart']])) etaStart = hyper[['etaStart']] else etaStart = 1
        ## Setting up hyper sampling tunning parameter
        if(!is.null(hyper[['hostSamplingSD']])) y_sd = hyper[['hostSamplingSD']]
        if(!is.null(hyper[['parasiteSamplingSD']])) w_sd = hyper[['parasiteSamplingSD']]
        if(!is.null(hyper[['etaSamplingSD']])) eta_sd = hyper[['etaSamplingSD']]
        print('host Hyper')
        print(hyper[['hostHyper']])
                print('parasite Hyper')
        print(hyper[['parasiteHyper']])
    }
    
    
    ls = floor(log(tol.err + w_sd)/beta)
    lsy = floor(log(tol.err + y_sd)/beta)
    lseta = floor(log(tol.err + eta_sd)/beta)
    batch.size = 50
    
    ## Variable holders
    hh<-matrix(0, nrow=4, ncol=burn_in+1)
    hh[1,]<-a_y; hh[2,]<-b_y
    hh[3,]<-a_w; hh[4,]<-b_w
    
    g0<-rep(0, burn_in+1)
    y0<-matrix(y, nrow= n_y, ncol =burn_in+1)
    w0<-matrix(w, nrow= n_w, ncol =burn_in+1)

    subItra = if(ICM.HORIZ) n_y else n_w
    g0in<-rep(0, subItra+1)
    y0in<-matrix(y, nrow= n_y, ncol =subItra+1)
    w0in<-matrix(w, nrow= n_w, ncol =subItra+1)
    petain = rep(etaStart, subItra+1)

	Z0 <- Z==0
    U0 <- Z*0
	mc <- colSums(1*(Z>0))
	mr <- rowSums(1*(Z>0))
    if(!missing(eta)) peta = rep(eta, burn_in)  else
    peta = rep(etaStart, burn_in)
    
    if(missing(dist)) pdist.right<-pdist<- matrix(1, nrow(Z), ncol(Z)) else
    pdist = (dist^peta[1])%*%Z
    tryCatch(
        for(s in 1:slice){
            if(s%%100==0)
                print(sprintf('Slice %d', s))
            g0in[1]= g0[s]
            y0in[,1] = y0[,s]
            w0in[,1] = w0[,s]
            petain[1] = peta[s]
            if(!missing(dist)){
                dist1=0*dist
                dist1[upper.tri(dist1)] =dist[upper.tri(dist)]^peta[s]
                dist1 = dist1 + t(dist1)
                pdist = (dist1)%*%Z
            }
            if(AdaptiveMC)              # Updaing tunning parameters
                if((((s+1)%%batch.size==0) &  (s < 0.5*burn_in))){
                    if(!distOnly){
                        ls  = AdaptiveSigma(w0, ls, s+1)
                        w_sd = exp(beta*ls)
                        
                        lsy = AdaptiveSigma(y0, lsy, s+1)
                        y_sd = exp(beta*lsy)
                    }
                    if(!missing(eta) & !missing(dist)){
                        lseta = AdaptiveSigma(peta, lseta, s+1)
                        eta_sd = exp(beta*lseta)
                    }
                }

            ## Updatting latent scores
            if(!uncertain){
                U0<-rExp(pdist*outer(y0[,s],w0[, s]))
                U0[Z0]<-1
            }else
                U0 <-rExp2(pdist*outer(y0[,s],w0[, s]),g0[s], Z, Z0)
            
            for (i in 1:subItra){
                ## needed only for digonal
                if(ICM.HORIZ){
                    if(!distOnly){                # Updating parasite parameters
                        ## Updating host parameters
                        y0in[i, 1]<-raffinity.MH(y0in[i,1],mr[i],
                                                 (U0[i,]*pdist[i,])%*%w0in[,i],
                                                 sig=y_sd, c(hh[1,s], hh[2,s]))
                        ## Updating the parasite parameters
                        w0in[, i+1]<-raffinity.MH(w0in[,i],Z[i,],
                                                  y0in[i,1]*(U0[i,]*pdist[i,]),
                                                  sig=w_sd, c(hh[3,s], hh[4,s]))
                    }
                    ## Updating similarity matix parameter
                    if(!missing(eta) & !missing(dist)){
                        new.eta = rEta.Horiz(petain[i],dist,pdist[i,],i,Z,y0in[i,1],
                            w0in[,i+1],U0[i,], eta_sd=eta_sd)
                        petain[i+1] = new.eta$eta
                        pdist[i,] = new.eta$dist
                    }
                    if(uncertain)
                        g0in[i+1] = rg(Z[i,], U0[i,],l=pdist[i,]*y0in[i]*w0in[,i+1])
                    
                }else{
                    Diag =  cbind(1:n_y, 1+(1:n_y-1 + (i-1)) %% n_w)
                    if(!distOnly){     
                        ## Updating host parameters
                        y0in[,i+1]<-raffinity.MH(y0in[,i],Z[Diag],
                                                 U0[Diag]*pdist[Diag]*w0in[Diag[,2],i],
                                                 sig=y_sd, c(hh[1,s], hh[2,s]))
                        ## Updating the parasite parameters
                        w0in[Diag[,2],i+1]<- raffinity.MH(w0in[Diag[,2],i],Z[Diag],
                                                          y0in[,i+1]*U0[Diag]*pdist[Diag],
                                                          sig=w_sd, c(hh[3,s], hh[4,s]))
                        w0in[-Diag[,2], i+1]<-0
                        
                    }
                    ## updating eta
                    if(!missing(eta) & !missing(dist)){
                        new.eta = rEta.Diag(petain[i],dist,pdist,Diag,Z,y0in[,i+1],
                            w0in[Diag[,2],i+1],U0[Diag], eta_sd=eta_sd)
                        petain[i+1] = new.eta$eta
                        pdist = new.eta$dist
                    }
                    ## Uncertainty variable
                    if(uncertain)
                        g0in[i+1] = rg(Z, U0,l=pdist*outer(y0in[,i+1],w0in[,i+1])) 
                }
                ## end of slice loop
            }

            if(!distOnly){
                y0[,s+1]<-if(ICM.HORIZ) y0in[,1] else rowSums(y0in[,-1])/subItra
                w0[,s+1]<- rowSums(w0in[,-1])/subItra

                y0in<-matrix(y, nrow= n_y, ncol =subItra+1)
                w0in<-matrix(w, nrow= n_w, ncol =subItra+1)
            }

            peta[s+1]<-mean(petain[-1])
            petain = rep(etaStart, subItra+1)

            g0[s+1]= mean(g0in[-1])
            g0in<-rep(0, subItra+1)

            
            ## Hyper parameters
            if(updateHyper & !distOnly){
                ## ## Updating hyper parameters for Parasite
                rHyper = if(ICM.HORIZ) rHyper.Horiz else rHyper.Diag
                new = rHyper(c(hh[3,s], hh[4,s]),y0[,s+1]%*%(U0*pdist), mc)
                hh[3, s+1] = new[1]
                hh[4, s+1] = new[2]
                ## Updating hyper parameters for Hosts
                new = rHyper(c(hh[1,s], hh[2,s]),(U0*pdist)%*%w0[,s+1],mr)
                hh[1, s+1] = new[1]
                hh[2, s+1] = new[2]
            }
        }
       ,warning = function(w)
           {print(c('warning at i:', i)); print(w);traceback()},
        error =function(e)
            {print(c('warning at i:', i)); print(e);traceback()} ,
        finally = print("Done!"))
    ## throwing out the burn_in stage.(approx 30%)
    if(throw.out==0) throw.out = c(1:burn_in) else throw.out = -c(1:throw.out)
    y0 =  y0[,throw.out] 
    w0 =  w0[,throw.out]
    if(!missing(eta)){
        eta = peta[throw.out]
    }else eta = NULL
    g = if(uncertain) t(data.frame(g = g0[throw.out])) else NULL
    param_phy = list(w = w0, y = y0, burn_in = burn_in - max(-throw.out), throw.out = max(-throw.out),eta = eta, g=g, hh=hh, sd = list(w=w_sd, y = y_sd, eta= eta_sd))
    param_phy
}
##################################################
##################################################
## old functions
