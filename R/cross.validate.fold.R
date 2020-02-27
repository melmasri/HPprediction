cross.validate.fold <-
function(Z, n= 10, min.per.col = 1, missing.pattern=c('random','prop.to.col.sums')){
    ## n-fold cross validation
    ## Returns a matrix of 3 columns, the first two are the (row,col) index of the pair,
    ## the third is the group
    missing.pattern = tolower(missing.pattern[1])
    if(max(range(Z))>1) Z[Z>0]<-1
    pairs = which(Z==1, arr.ind=T)
    colnames(pairs)<-c('row', 'col')
    
    if(length(which(colSums(Z)<min.per.col))>0){
        aux = which(pairs[,'col'] %in% which(colSums(Z)<min.per.col))
        if(length(aux))
            pairs = pairs[-aux,]
    }
    
    colm = pmax(colSums(Z) -min.per.col , 0)
    size = floor(sum(colm)/n)
    gr = rep(size, n)
    if(sum(colm) %% size!=0)
        gr[n] =  gr[n] + sum(colm) %% size
    
    group.colm = rep(1:n,times = gr)[sample.int(sum(colm), sum(colm))]
    pair.list = numeric(sum(colm))
    for(i in 1:sum(colm)){
        a = which(colm>0)
        if(missing.pattern=='random')
            b = a[sample.int(length(a),1)] else
        if (missing.pattern=='prop.to.col.sums')
            b = a[sample.int(length(a),1, prob=colm[a]/sum(colm[a]))] else
        stop('missing pattern has to be specified from selection!')
        colm[b] = colm[b]-1
        pair.list[i]<-b
    }
    pair.list= tapply(pair.list, group.colm,identity)
    
    gr.list= list()
    bank= c()
    for(i in 1:n){
        a= table(pair.list[[i]])
        gr.rows = unlist(sapply(1:length(a), function(r){
            b = which(pairs[,'col']== as.numeric(names(a[r])))
            b =setdiff(b, bank)
            b[sample.int(length(b), a[r])]
        }))
        bank = c(bank, gr.rows)
        gr.list[[i]]<-cbind(gr.rows, i)
    }

    aux = do.call('rbind', gr.list)
    pairs = cbind(pairs[aux[,1], ],gr= aux[,2])
    
    print(sprintf("Actual cross-validation rate is %0.3f" , table(pairs[,'gr'])/sum(1*(Z>0))))
    pairs[order(pairs[,'gr']),]
    
}
