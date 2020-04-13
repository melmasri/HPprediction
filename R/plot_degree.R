#' A function to plot degree distributions  
#'
#' @param Z bipartite interaction matrix
#' @param Z_est the posterior bipartite interaction matrix (optional)
#' @param type whether to plot 'hosts', 'parasites', or 'both'
#' @param host.col  Colour to use for hosts (default is 'blue')
#' @param parasite.col Colour to use for parasites (default is 'red')
#' 
#' @description
#' 
#' This function generates a degree distribution plot for the binary interaction matrix 'Z'.
#' Optionally, you can compare the observed matrix to the posterior matrix by including it as 'Z_est'.
#' By default the function plots the degree for both rows (hosts) and columns (parasites), but you can modify this with the parameter 'type'.
#' 
#'
#' @examples
#' 
#' # Simluate a Z matrix and plot the degree distributions
#' 
#' Z <- matrix(rbinom(50*200, 1, 0.01), nrow=50, ncol=200)
#' Z <- Z[,colSums(Z)>0]
#' plot_degree(tree)
#'  
#' @export
#' 
plot_degree <-
function(Z, Z_est, type='both', host.col='blue', parasite.col='red', cex.lab=1.5, cex.axis=2, pt.cex=1.5, legend.cex=2){
    ## Plots the degree distribution per marginal on a bipartite biadjacency matrix
    ## input
    ## Z - interaction matrix
    ## Z_est = posterior interaction matrix (optional)
    ## type  = hosts, parasites, both
    ## extra input for colours 

    ## Optional estimated presence/absence matrix (Z_est) can be added to existing plot.
    para_degrees <- as.data.frame(table(colSums(Z)))
    para_degrees$Var1 <- as.numeric(para_degrees$Var1)
    ## para_degrees = para_degrees[-which(para_degrees$Var1<2),]
    host_degrees <- as.data.frame(table(rowSums(Z)))
    host_degrees$Var1 <- as.numeric(host_degrees$Var1)
    ## host_degrees = host_degrees[-which(host_degrees$Var1<2),]

    xlim = c(1, max(para_degrees$Var1,host_degrees$Var1)*1.5)
    ylim = c(1, max(para_degrees$Freq,host_degrees$Freq)*1.5)

    if (!missing(Z_est)){
        para_est <- as.data.frame(table(colSums(Z_est)))
        para_est$Var1 <- as.numeric(para_est$Var1)
        ## para_est = para_est[-which(para_est$Var1<2),]
        host_est <- as.data.frame(table(rowSums(Z_est)))
        host_est$Var1 <- as.numeric(host_est$Var1)
        ## host_est = host_est[-which(host_est$Var1<2),]
        xlim = c(1, max(para_degrees$Var1,host_degrees$Var1,
            para_est$Var1, host_est$Var1)*1.5)
        ylim = c(1, max(para_degrees$Freq,host_degrees$Freq,
            para_est$Freq, host_est$Freq)*1.5)
    }
    gpch = c('+', '*')
    if(type=='parasites'){
        plot((para_degrees), type="p", col=parasite.col, pch=gpch[2], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = cex.lab , cex.axis = cex.axis)
        legend(xlim[2]*0.2, ylim[2]*1.4, c("Parasites"), col = parasite.col,
               pch = gpch[2], bty="n",pt.cex=pt.cex, cex=legend.cex)
        if(!missing(Z_est)){
            points((para_est), type="p", col=parasite.col, pch=16)
            legend(xlim[2]*0.2, ylim[2]*0.8, c("Est"), col = parasite.col,
               pch = 16, bty='n', pt.cex=pt.cex, cex=legend.cex)
        }
    }
    if(type=='hosts'){
        plot((host_degrees), type="p", col=host.col, pch=gpch[1], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = cex.lab, cex.axis = cex.axis)
        legend(xlim[2]*0.2, ylim[2]*1.4, c("Hosts"), col = host.col,
               pch = gpch[1], bty="n", pt.cex=pt.cex, cex=legend.cex)
        if(!missing(Z_est)){
            points((host_est), type="p", col=host.col, pch=16)
            legend(xlim[2]*0.2, ylim[2]*0.8, c("Est"), col = host.col,
                   pch = 16, bty='n',pt.cex=pt.cex, cex=legend.cex)
        }
    }
    if(type=='both'){
        plot((para_degrees), type="p", col=parasite.col, pch=gpch[2], log="xy", xlim=xlim, ylim=ylim, ylab="Number of nodes", xlab="Degree", cex.lab = cex.lab ,cex.axis = cex.axis)
        points((host_degrees), type="p", col=host.col, pch=gpch[1])
    legend(xlim[2]*0.2, ylim[2]*1.4, c("Parasites", "Hosts"), col = c(parasite.col, host.col),
           pch = gpch[2:1], bty = 'n', pt.cex=pt.cex, cex=legend.cex)
        if (!missing(Z_est)) {
            points((para_est), type="p", col=parasite.col, pch=16)
            points((host_est), type="p", col=host.col, pch=16)
            legend(xlim[2]*0.2, ylim[2]*0.33, c("Est"), col = c("black"),
                   pch = 16, bty='n', pt.cex=pt.cex, cex=legend.cex)
        }
    }
}
