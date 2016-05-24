param_phy = gibbs_one(1*(com>0),slice=5,dist= phy_dist, eta=1, uncertain=FALSE, yMH=TRUE, wMH = TRUE, hyper=hyper, wEta = FALSE, yEta=FALSE)



param_phy = gibbs_one(1*(com>0),slice=5,dist= phy_dist, eta=1, uncertain=FALSE, yMH=TRUE, wMH = TRUE, hyper=hyper, wEta = FALSE, yEta=FALSE)
BAD,


com_pa = 1*(com>0)
r = which.max(rowSums(com_pa));r
c = which.max(colSums(com_pa));c
burn = param_phy$burn_in - 20000:1
dim(param_phy$y)
burn = burn[burn>0]
range(burn)


plot(param_phy$w[c,])
abline(h=w[c])
w[c]


plot(param_phy$y[r,])
abline(h=y[r])
y[r]
plot(param_phy$hh[1,])

lines(param_phy$hh[3,], col='blue')


Vulpes_vulpes 
          219 
mean(param_phy$y[r,])
[1] 0.1066
mean(param_phy$w[c,])
[1] 31.89
> 

