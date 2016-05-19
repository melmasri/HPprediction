#!/bin/R
##################
tol.err=1e-5
## Gibbs sampling based on sections 2.4 and 2.5
rExp<-function(l,a=1){
	# Sampling from a truncated exponential distribution
	#l = outer(c(1,2),c(2,3))
	unif = runif(length(l))
	-log(1- unif*(1-exp(-l*a)))/(l + tol.err)
}

rExp2<-function(l, l1, Z, Z0){
    p = 1- exp(-l)
    unif = matrix(runif(length(l)), dim(p))
    ## Method one P(z=0|g) = 1\delta_{s<0} + g\delta_{s>=0}
    U = 1 + 0*p
    U[!Z0] = -log(1- unif[!Z0]*p[!Z0])/(l[!Z0] + tol.err)
    aa = (unif < l1*p/(l1*p + 1-p)) & Z0
    U[aa] =  -log(1 - unif[aa]*(l1*p[aa] + 1-p[aa])/l1)/l[aa]
    ## End of Method one
    ## ## U = -log(unif*(1-Z) + (1-unif)*Z )/(l+tol.err)
    ## aa = 1*(runif(sum(Z0))<= ((1-p)/(1-p + l1*p))[Z0])
    ## U[Z0]<-    U[Z0]^(1-aa) 
    ## U
    ## Method two P(z=1|g) = g\delta_{s<0} + 1\delta_{s>=0}
    ## U = 1 + 0*p
    ## aa = (unif < p/(p + l1*(1-p))) & !Z0
    ## U[aa] = -log(1- unif[aa]*(p[aa] + l1*(1-p[aa])))/(l[aa])
    U
}

rW<-function(m, U,y){
	#Gibbs sampling for a popularity parameters (W)
	#See section 2.5
	alpha = m
	alpha[m!=0] = alpha[m!=0]- sigma
	beta  = c(1/tau + y%*%U)	 # not sure it is 1/beta or beta
	rgamma(length(alpha),shape = alpha + tol.err, scale = 1/beta)
}

rWmc<-function(w.old, y,eta=1, pdist, mc, U, w_sd=0.1, wEta){
    a = a_w # 0.2 #shape parameter of a Gamma prior
    b = b_w  # 1 # Scale paramter of a Gamma prior
    if(wEta) eta = 1
    w.old = w.old^(1/eta)
    epsilon = rnorm(length(w.old), 0, sd = w_sd) # 0.2 for regular 0.12 tested with AdaptiveMCMC
    w.prop = w.old +  -2*(w.old + epsilon <=0)*epsilon + epsilon
    likeli = (eta*mc+a-1)*log(w.prop/w.old) - (y%*%(U*pdist))*(w.prop^eta - w.old^eta) - (w.prop - w.old)*b
    ratio = pmin(1, exp(likeli))
    u = 1*(runif(length(w.old))<=ratio)

    u*w.prop^eta + (1-u)*w.old^eta
}

rY<-function(w0,w_star, U0,m_r){
    ## Based on first formula of section 3.
    ## a_y =0.2;b_y = 1
    beta  = c(1/b_y + U0 %*% w0  + w_star)
    rgamma(nrow(U0), shape = a_y + m_r, scale = 1/beta)
}

rYmc<-function(y.old, w,eta=1, pdist, mr, U, y_sd, yEta){
    a = a_y# shape parameter of a Gamma prior
    b = b_y# Scale paramter of a Gamma prior
    if(yEta) eta=1
    y.old = y.old^(1/eta)
    epsilon = rnorm(length(y.old), 0, sd = y_sd) # 0.2 for regular MH 0.6 tested with Adaptive MCMC
    y.prop = y.old +  -2*(y.old + epsilon <=0)*epsilon + epsilon
    
    likeli = (eta*mr+a-1)*log(y.prop/y.old) - (y.prop^eta - y.old^eta)*((U*pdist)%*%w) - (y.prop - y.old)*b

    ratio = pmin(1, exp(likeli));ratio
    u = 1*(runif(length(y.old))<=ratio)
    u*y.prop^eta + (1-u)*y.old^eta
}

rGumble.single<-function(mu, beta=1){
    mu = log(mu) 
    u = matrix(runif(length(mu)), dim(mu))
    -beta*log(-log(u)) + mu
}
    
rGumble<-function(w, y,beta = 1){
    ## Here we set s_ij to a latent score of a Gumble distribution, s.t s_ij  ~ Gumble( mu =  log(w_j), Beta =1 )
	## Initially it was supposed to be a Truncated Gumble by x_ij = min(0, s_ij + log(y_i)).
	## Hence this function would return the x_ij and s_ij.
	## Truncated Gummble
	 l = outer(y,w);
	 unif = matrix(runif(length(l)), dim(l))
	 log(l)	 -  log(-log(unif +exp(-l)*(1-unif)  ))
 }

rEta<-function(eta.old, dist , pdist,Z, y, w, U,md, eta_sd =0.01, wEta, yEta){
    ## A function that generates the prior for the power of distance
    ## eta.old = peta[i];y=y0[,i+1];w= w0[,i+1];U= U0
    pdist.old = (dist^eta.old)%*%Z  +  md^eta.old
    ## for (t in 1:5){
    epsilon = rnorm(1, 0, sd = eta_sd) # Tested to be 0.005 with Adaptive MCMC
    eta.prop = eta.old + if( eta.old + epsilon <=0) -1*epsilon else epsilon
    if(wEta) w.new= w^(eta.prop/eta.old) else w.new = w
    if(yEta) y.new =y^(eta.prop/eta.old) else y.new = y
	## w.new =w # when eta is 1.
        pdist.new = (dist^eta.prop)%*%Z + md^eta.prop
        likli = sum(log((pdist.new/pdist.old)^Z )) - sum(U*(outer(y.new,w.new)*pdist.new - pdist.old*outer(y,w)))
        prior = dgamma(eta.prop, shape=a_e, scale=b_e, log=TRUE) - dgamma(eta.old, shape = a_e,scale = b_e, log=TRUE)
        ratio = min(1, exp(likli + prior));ratio
        u = (runif(1)<=ratio)
        if(u) { eta.old  = eta.prop;pdist.old = pdist.new }
    ##}
    list (eta=eta.old, dist=pdist.old)
}

rEtaOptim<-function(eta, eta.old, Z, dist,w,y,U,md){
    ##w= w0[,i+1];y = y0[,i+1]; U=U0; eta.old= peta[i];md = min_dist
    w.new= w^(eta/eta.old)
    y.new =y^(eta/eta.old)
    ## w.new = w  # to set eta =1
    pdist.new = (dist^eta)%*%Z + md^eta
    likli = sum(log(pdist.new^Z)) - y.new%*%(pdist.new*U)%*%w.new
    prior = dgamma(eta, 1,1, log=TRUE)
    -(likli + prior)
}

rL<-function(Z,m, Y, sumY, l){
    ## Method 1 P(Z=0|g) = 1 \delta(s<0) + g \delta(s>=0)
    ## S = rGumble.single(l)
    S = log(l) - 0.5772
    ZZ = 1*(S>=0)
    ## Regular method
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
    c(g, 1)
    ## End of method 2
    ## alpha11 = 1
    ## yy = Y*Z
    ## yy = sum(1/yy[yy!=0]) + 1/alpha11
    ## c(rexp(1, 1/yy ),1)
        
    
}

gibbs_one<-function(Z,y,w,dist, slice = 10, eta, uncertain =FALSE,wMH=FALSE, yMH=FALSE, wEta = TRUE, yEta =FALSE){
	## A one step update in a Gibbs sampler.
	## ## initialize
    ## Z = 1*(comCross>0); slice =5 ;dist= phy_dist;
    ##Z = com_pa; slice=5 ;dist= phy_dist; eta=1; wMH = TRUE;uncertain =FALSE; yMH=FALSE    
   
    latent='exp'
	n_w = ncol(Z);n_y = nrow(Z)
    n = n_w*n_y
    burn_in = n_w*slice
    throw.out.slice = floor(0.2*slice)
    throw.out = throw.out.slice*ncol(Z)
    print(sprintf("Run for %i slices, and %i burn ins",slice, burn_in))
    print(sprintf("Settings: uncertain=%s, wmH=%s,yMH=%s,wEta=%s,yEta=%s .",uncertain, wMH, yMH,wEta,yEta))
    print(dim(Z))
    
    ## binary = TRUE
    ## if(max(range(Z))>1) binary = FALSE
    ## if(!binary){
    ##     ## When count data is available
    ##     Y = Z
    ##     ## Z = matrix(runif(prod(dim(Z)))<=0.1 + 0, nrow(Z), ncol(Z))
    ##     Z = 1*(Z>0)
    ##     L0<-rep(1, burn_in)
    ##     L1<-rep(1, burn_in)
    ##     sumY = sum(Y)
    ## }
    L0<-rep(0.5, burn_in)
    L1<-rep(0.5, burn_in)
    
    if(!missing(y)) y0<-matrix( y, nrow= n_y, ncol =burn_in) else
    y0<-matrix(1,nrow= n_y, ncol =burn_in)
    if(!missing(w)) w0<-matrix(w, nrow= n_w, ncol =burn_in) else 
    w0<-matrix(1, nrow= n_w, ncol =burn_in)
	w_star = 0
	Z0 <- Z==0
	m_c <- colSums(1*(Z>0))
	m_r <- rowSums(1*(Z>0))
    if(!missing(dist))
        min_dist  = min(dist[upper.tri(dist)]) else min_dist=0
    if(!missing(eta))  petashoot=peta = rep(eta, burn_in)  else
    petashoot=peta = rep(1, burn_in)
    if(missing(dist)) pdist = 1 else
    pdist = (dist^peta[1])%*%Z
    ls = lsy = 0
    lseta = -30
    batch.size = 50
    tryCatch(
        
        for (i in 1:(burn_in-1)){
        ##for (i in 1:99){
            ## # Updating scores
            if(i%%batch.size==0){
                batch = (floor(i/batch.size)-1)*batch.size + 2
                pvar= apply(w0[,5:i],1,sd)
                ac = 1-rowMeans(1*abs(w0[,batch:i]- w0[,batch:i-1])<pvar/20)
                ls = ls + sign(ac - 0.44)*(1*(abs(ac - 0.44)>0.03))
                w_sd = exp(0.1*ls)
                
                pvar= apply(y0[,5:i],1,sd)
                ac = 1-rowMeans(1*abs(y0[,batch:i]- y0[,batch:i-1])<pvar/20)
                lsy = lsy + sign(ac - 0.44)*(1*(abs(ac - 0.44)>0.03))
                y_sd = exp(0.1*lsy)
                if(!missing(eta)){
                    pvar= sd(peta[5:i])
                    batch=5
                    ac = 1-mean(1*abs(peta[batch:i] - peta[batch:i -1])<1e-4)
                    lseta = lseta + sign(ac - 0.44)*(1*(abs(ac - 0.44)>0.03))
                    eta_sd = exp(0.2*lseta)
                    ## print(sprintf('AC %0.3f - L %0.0f eta_sd %f', ac, lseta, eta_sd))
                }
            }
            if(!uncertain){
                if(missing(dist)){
                    U0<- rExp(outer(y0[,i],w0[,i]))
                }else{
                    U0<- rExp(pdist*outer(y0[,i],w0[,i]))  
                }
                U0[Z0]<-1 # 0.03162  seconds
            }else
                U0 <-rExp2(pdist*outer(y0[,i],w0[,i]), L1[i], Z, Z0)
            
            if(wMH)
                w0[,i+1] <-rWmc(w0[,i], y0[,i],peta[i], pdist, m_c, U0, w_sd =w_sd, wEta) else 
            rW(m_c, U0,y0[,i])

            ##w_star <- rgamma(1, shape=a, scale = 1/(tau  + sum(y0[,i])))  #What is left of the mass G(\theta)*
            w_star<-0
            ## Updating col parameter
            ##if(missing(y))
            if(yMH) y0[,i+1] <-rYmc(y0[,i], w0[,i+1],peta[i], pdist, m_r, U0,y_sd,yEta) else 
            y0[,i+1]<-rY(w0=w0[,i+1], w_star = w_star, U0 =U0*pdist, m_r) # 0.01406 seconds

            if(!missing(eta)){  # 0.38804 seconds, 0.27028  for optim by itself, 0.11 seconds for what is left
                ## Testing for minimum 
                ## rEtaIn <-function(eta)
                ## rEtaOptim(eta, peta[i], Z, dist, w0[,i+1], y0[,i+1], U0,min_dist)
                ## ## new.eta= optim(peta[i], rEtaIn, method = "Brent",lower = tol.err, upper = 5)
                ## new.eta = nlm(rEtaIn, peta[i],gradtol=1e-2, stepmax = 0.5,steptol = 1e-4)
                ## petashoot[i+1]<-new.eta$par
                ## petashoot[i+1]<-new.eta$estimate
                new.eta = rEta(peta[i],dist,pdist,Z,y0[,i+1], w0[,i+1], U0,min_dist, eta_sd=eta_sd, wEta, yEta)
                peta[i+1] = new.eta$eta
                if(wMH & wEta)
                    w0[,i+1]=w0[,i+1]^(peta[i+1]/peta[i])
                if(yMH & yEta)
                    y0[,i+1]=y0[,i+1]^(peta[i+1]/peta[i])
                pdist = new.eta$dist
            }
            if(uncertain){
                L10 = rL(Z,n, U0, sumY, l=pdist*outer(y0[,i],w0[,i]))  # 0.00768 seconds
                L1[i+1]<-L10[1]
                L0[i+1]<-L10[2]
            }
        }
       
       ,warning = function(w){print(c('warning at i:', i)); print(w);traceback()},
	error =function(e) {print(c('warning at i:', i)); print(e);traceback()} , finally = print("Done!"))
    ## throwing out the burn_in stage.(approx 30%)
    if(throw.out==0) throw.out = c(1:burn_in) else throw.out = -c(1:throw.out)
    w_star = w_star[throw.out]
    w0 = w0[,throw.out]
    y0 = if(missing(y)) y0[,throw.out] else y
    etashoot = if(!missing(eta)) petashoot[throw.out] else NULL
    if(!missing(eta)){
        eta = peta[throw.out]
        if(wEta) w0 = w0^eta
        if(yEta) y0 = y0^eta
    }else eta = NULL
    L = if(uncertain) t(data.frame(l1 = L1[throw.out], l0 = L0[throw.out])) else NULL
    param_phy = list(w_star  =w_star, w = w0, y = y0, burn_in = burn_in - max(-throw.out), throw.out = max(-throw.out),eta = eta, L=L, etashoot = etashoot)
    param_phy
}
##################################################
##################################################
