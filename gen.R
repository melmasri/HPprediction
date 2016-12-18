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


rHyper<-function(old,y, w, U, pdist, mr, mc, type='host', sig = 0.1){
    ## old= c(a_y, b_y);y=y0[,i+1];w=w0[,i+1];U= U0; pdist=pdist.right; type='host'
    ## sig=0.1
    new = old*exp(rnorm(length(old), 0, sd = sig ))
    new[2]<-1
    if(type=='host'){
        u = (U*pdist)%*%w
        r = length(y)*(new[1]*log(new[2]) - lgamma(new[1]) )  +
            sum(lgamma(mr+new[1]) - (mr + new[1])*log(new[2] + u)) +
                -(length(y)*(old[1]*log(old[2]) - lgamma(old[1]))  +
                      sum(lgamma(mr+old[1]) - (mr + old[1])*log(old[2] + u)))
    }
    if(type=='parasite'){
        u =  y%*%(U*pdist)
        r = length(w)*(new[1]*log(new[2]) - lgamma(new[1]) )  +
            sum(lgamma(mc+new[1]) - (mc + new[1])*log(new[2] + u)) +
            -(length(w)*(old[1]*log(old[2]) - lgamma(old[1]))  +
                  sum(lgamma(mc+old[1]) - (mc + old[1])*log(old[2] + u)))
    }
    ratio = min(1, exp(r));ratio
    u = 1*(runif(1)<=ratio)
    u*new + (1-u)*old
}

rEta<-function(eta.old,dist,pdist.old,zdiag, Z, y, w, U, eta_sd =0.01){
    ## A function that generates the prior for the power of distance
    eta.prop = eta.old*exp(rnorm(1, 0, sd = eta_sd))
    dist1=0*dist
    dist1[upper.tri(dist1)] =dist[upper.tri(dist)]^eta.prop
    dist1 = dist1+ t(dist1)
    pdist.new = dist1%*%Z
    likeli = sum(log((pdist.new[zdiag]/pdist.old[zdiag])^Z[zdiag] ))-
        sum(U*y*w*(pdist.new[zdiag] - pdist.old[zdiag]))
    
    ratio = min(1, exp(likeli));ratio
    u = (runif(1)<=ratio)
    if(u) {eta.old  = eta.prop;pdist.old = pdist.new}
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

gibbs_one<-function(Z,dist, slice = 10, eta,hyper, uncertain =FALSE,updateHyper=FALSE, AdaptiveMC=FALSE, distOnly=FALSE){
    ## Z=1*(com>0);slice=slice;dist=phy_dist; eta=1; hyper=hyper; updateHyper = updateHyper; AdaptiveMC=AdaptiveMC
    ## A one step update in a Gibbs sampler.
	## ## initialize
    #Z = 1*(Z>0)
    colnames(Z)<-NULL
    rownames(Z)<-NULL
    if(!missing(dist))
        colnames(dist)<-rownames(dist)<-NULL
	n_w = ncol(Z);n_y = nrow(Z)
    n = n_w*n_y
    burn_in = slice
    throw.out = floor(0.2*slice)
    print(sprintf("Run for %i slices, and %i burn ins",slice, burn_in))
    print(sprintf("Settings: uncertain=%s.",uncertain))
	print(sprintf("updateHyper= %s,AdaptiveMC= %s ", updateHyper, AdaptiveMC))
    print(dim(Z))
    ## Hyper parameters
    a_w = 1;b_w = 1; a_y = 1; b_y = 1;

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
        if(distOnly){
            y=1
            w=1
        }
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
    g0<-rep(0.5, burn_in+1)
    y0<-matrix(y,nrow= n_y, ncol =burn_in+1)
    w0<-matrix(w, nrow= n_w, ncol =burn_in+1)

    g0in<-rep(0.5, n_w+1)
    y0in<-matrix(0,nrow= n_y, ncol =n_w+1)
    w0in<-matrix(0, nrow= n_w, ncol =n_w+1)
    petain = rep(1, n_w+1)
	Z0 <- Z==0
    U0 <- Z*0
	mc <- colSums(1*(Z>0))
	mr <- rowSums(1*(Z>0))
    if(!missing(dist))
        min_dist  = min(dist[upper.tri(dist)]) else min_dist=0
    if(!missing(eta)) peta = rep(etaStart, burn_in)  else
    peta = rep(1, burn_in)
        
    if(missing(dist)) pdist.right<-pdist<- 1 else
    pdist = dist%*%Z
    tryCatch(
        for(s in 1:slice){
            print(sprintf('Slice %d', s))
            g0in[1]= g0[s]
            y0in[,1] = y0[,s]
            w0in[,1] = w0[,s]
            system.time(
            for (i in 1:n_w-1){
                Diag = cbind(1:n_y, 1+ (1:n_y-1 + i) %% n_w)
                if(AdaptiveMC)
                    if((((i+1)%%batch.size==0) & (s < 0.4*burn_in))){
                        if(!distOnly){
                            #ls  = AdaptiveSigma(w0, ls, i)
                            #w_sd = exp(beta*ls)
                            
                            lsy = AdaptiveSigma(y0in, lsy, i+1)
                            y_sd = exp(beta*lsy)
                        }
                        if(!missing(eta) & !missing(dist)){
                            lseta = AdaptiveSigma(petain, lseta, i+1)
                            eta_sd = exp(beta*lseta)
                        }
                    }

                ## Updatting latent scores
                if(!uncertain){
                    U0[Diag]<-rExp(pdist[Diag]*y0in[,i+1]*w0in[Diag[,2], i+1])
                    U0[Diag[Z0[Diag],]]<-1
                }else
                    U0[Diag] <-rExp2(
                        pdist[Diag]*y0in[,i+1]*w0in[Diag[,2],i+1],
                        g0in[i+1], Z[Diag], Z0[Diag])
                
                ## Updating parasite parameters
                if(!distOnly){
                    w0in[Diag[,2],i+2]<- raffinity.MH(w0in[Diag[,2],i+1],Z[Diag],
                                              y0in[,i+1]*U0[Diag]*pdist[Diag],
                                              sig=w_sd, c(a_w, b_w))
                    ## Updating host parameters
                    y0in[,i+2]<-raffinity.MH(y0in[,i+1],Z[Diag],
                                             U0[Diag]*pdist[Diag]*w0in[Diag[,2],i+2],
                                             sig=y_sd, c(a_y, b_y))
                }
                
                ## updating eta
                if(!missing(eta)){
                    new.eta = rEta(petain[i+1],dist,pdist,
                        Diag,Z,y0in[,i+2], w0in[Diag[,2],i+2],
                        U0[Diag], eta_sd=eta_sd)
                    petain[i+2] = new.eta$eta
                    pdist.right = new.eta$dist
                    ##pdist = (dist^petain[i+2])%*%Z
                    ##if(i%%100==0) print(peta[i+1])
                }
                
                ## Uncertainty variable
                ## if(uncertain){
                ##     g0in[i+1] = rg(Z, U0,l=pdist*outer(y0in[,i+1],w0in[,i+1])) 
                ## }
                
                ## end of slice loop    
            }
            )
            y0[,s+1]<-rowMeans(y0in[,-1])
            peta[s+1]<-mean(petain[-1])
            w0[,s+1]<- rowSums(w0in[,-1])/n_w
            
            g0in<-rep(0.5, n_w+1)
            y0in<-matrix(0,nrow= n_y, ncol =n_w+1)
            w0in<-matrix(0, nrow= n_w, ncol =n_w+1)
            petain = rep(1, n_w+1)
            
            ## ## Updating Hyper parameters
            ## if(updateHyper & !distOnly){ 
            ##     new = rHyper(c(a_y, b_y),y0[,i+1],w0[,i+1],
            ##         U0, pdist.right, mr, mc, type='host')
            ##     a_y = new[1]
            ##     b_y = new[2]
            ##     hh[1,i+1]<-a_y
            ##     hh[2,i+1]<-b_y
                
            ##     new = rHyper(c(a_w, b_w),y0[,i+1],w0[,i+1],
            ##         U0,pdist.right, mr, mc, type='parasite')
            ##     a_w = new[1]
            ##     b_w = new[2]
            ##     hh[3,i+1]<-a_w
            ##     hh[4,i+1]<-b_w
            ## }
            ##if(i%%100==0) print(c(a_y,b_y,a_w,b_w))
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
