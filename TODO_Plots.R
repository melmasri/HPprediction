# TODO Plots

rm(list=ls())
load('param.RData')

rho_col <- "blue"
gamma_col <- "red"
eta_col <- "darkgreen"


# Original Data
com_pa = 1*(com>0)

burn = param_phy$burn_in - 20000:1
burn = burn[burn>0]

# pdf(paste0(dataset, '_Z.pdf'))
plot_Z(1*(com>0) , xlab = 'parasites', ylab = 'hosts')
# dev.off()

# HIST: Average RHO
# pdf(paste0(dataset, '_rho_post.pdf'))
hist(rowMeans(param_phy$w[,burn]), breaks=50, col=rho_col,main ='', ylab = 'Frequency',
	xlab = expression(paste('Averge posterior of ', rho)), border=TRUE, xlim=c(23,32))
abline(v=quantile(rowMeans(param_phy$w[,burn]),probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
# box()
# dev.off()

# HIST: Average GAMMA
# pdf(paste0(dataset, '_gamma_post.pdf'))
hist(rowMeans(param_phy$y[,burn]), breaks=35, col=gamma_col,main ='', ylab = 'Frequency',
	xlab = expression(paste('Averge posterior of ', gamma)), border=TRUE)
abline(v=quantile(rowMeans(param_phy$y[,burn]),probs = 0.95), col="red", lty=2, lwd=2)
# box()
# dev.off()

# HIST: Average ETA
# pdf(paste0(dataset, '_eta_post.pdf'))
hist(param_phy$eta[burn], breaks=30, col=eta_col,main ='', ylab = 'Frequency',
	xlab = expression(paste('Posterior of ', eta)), border=TRUE, xlim=c(1.32, 1.48))
abline(v=quantile(param_phy$eta[burn],probs = c(0.05, 0.95)), col="red", lty=2, lwd=2)
# box()
# dev.off()


# TODO
# HIST: g for both datasets


#######################################
# BOXPLOT: Invidual RHOs
# pdf(paste0(dataset, '_boxplot_rho_100.pdf'))
aux = apply(param_phy$w[, burn],1, median )
aux = order(aux, decreasing = T);aux[1:100]
boxplot(data.frame(t(param_phy$w[aux[1:100], burn])), outline=F, 
	names = paste(1:100), ylab= 'Value', xlab = 'Ordered parameters',
	whiskcol=rho_col, staplecol=rho_col, col=rho_col, medcol="white")# BOXPLOT: Invidual RHOs
# dev.off()


# BOXPLOT: Invidual GAMMAs
# pdf(paste0(dataset, '_boxplot_gamma_100.pdf'))
aux = apply(param_phy$y[, burn],1, median )
aux = order(aux, decreasing = T);aux[1:100]
boxplot(data.frame(t(param_phy$y[aux[1:100], burn])), outline=F, 
	names = paste(1:100), ylab= 'Value', xlab = 'Ordered parameters',
	whiskcol=gamma_col, staplecol=gamma_col, col=gamma_col, medcol="white")
# dev.off()



##############################
# P HeatMap plot

# rm(list=ls())
load('gmp-01-05-00h46/param.RData')
# load('eid-01-05-00h46/param.RData')

com_pa = 1*(com>0)
burn = param_phy$burn_in - 20000:1
burn = burn[burn>0]
aux = getMean(param_phy)
com_pa = 1*(com>0)
P = 1-exp(-outer(aux$y^aux$eta, aux$w^aux$eta)*((phy_dist^aux$eta)%*%com_pa))

P_new <- P
P_new[com_pa==1] <- 1
reverse <- nrow(P_new) : 1
P_new <- P_new[reverse,]
image(t(P_new), col=colorRampPalette(c("white","black"))(400))
image(t(P_new), col=colorRampPalette(c("white","blue","brown"))(400))




######################################
# Tables
require(xtable)
rm(list=ls())
load("/home/max/GitHub/HP-prediction/April28/Uncertain-10foldCV-gmp-27-04-20h54/param.RData")
ls()

filename <- "test.tex"
roctemp= rocCurves(Z=1*(com10>0), Z_cross=com, P=P, all=TRUE, plot=FALSE, bins=400)
zz = 1*(roctemp$P>roctemp$threshold)
dim =dim(com)
tot.int = sum(com10>0)
zeros =sum(com10==0)
TP =sum(zz[com10>0])
FN =sum(1-zz[com10>0])
FP =sum(zz[com10==0])
TN = sum(1-zz[com10==0])
tb = matrix(c(TP, FN, FP, TN),2,2)
tb = cbind(tb, rowSums(tb))
tb = rbind(tb, colSums(tb))
colnames(tb)<-c('Positive', 'Negative', 'Total')
rownames(tb)<-c('Predicted positive', 'Predicted negative', 'Total')
print(xtable(tb,digits=0,  align='lccc'), file=filename)
