library(Renvlp)
library(expm)
library(readr)
library(data.table)
library(ggplot2)
library(plyr)


rrenv.choose_du = function(X,Y, Beta=NULL){
  
  red.rank.envelope<-function(X,Y,u,d){
    E<-env(t(X),t(Y),u)
    Gamma<-E$Gamma
    Gamma0<-E$Gamma0
    N<-dim(Y)[2]
    ifelse(u!=dim(Y)[1],Y.g<-t(Gamma)%*%as.matrix(Y),Y.g<-Y)
    r.g<-dim(Y.g)[1]
    Sx<-X%*%t(X)/N
    Sxyg<-t(X%*%scale(t(Y.g), scale = FALSE))/N
    Syg<-Y.g%*%t(Y.g)/N
    Cyg.x<-sqrtm(solve(Syg))%*%Sxyg%*%sqrtm(solve(Sx))
    D.truncated<-diag(c(svd(Cyg.x)$d[1:d],rep(0,length(svd(Cyg.x)$d)-d)))
    Cyg.x.d<-svd(Cyg.x)$u%*%D.truncated%*%t(svd(Cyg.x)$v)
    beta.r.e<-Gamma%*%sqrtm(Syg)%*%Cyg.x.d%*%solve(sqrtm(Sx))
    return(beta.r.e)
  }
  p = dim(X)[1]
  r = dim(Y)[1]
  max.u = r
  max.d = min(p,r)
  c = max.u-max.d
  rows = ((max.u*(max.u+1))/2) - ((c*(c+1))/2) - 1
  err.table = matrix(0, ncol=3, nrow=rows)
  cnt = 0
  
  if(is.null(Beta)){
    for(u in 2:max.u){
      d=0
      while(d<max.d & d<u){
        d = d+1
        cnt = cnt+1
        est.Betas = red.rank.envelope(X,Y,u,d)
        est.Ys = est.Betas %*% X
        err = norm(est.Ys-Y,type = "F")
        err.table[cnt,1]=d; err.table[cnt,2]=u; err.table[cnt,3]=err
      }
    }
    err.table = as.data.frame(err.table)
    names(err.table) = c('d','u','Error_in_Y')
    m <- ddply(err.table, .(d), transform, rescale = scales::rescale(Error_in_Y))
    p <- ggplot(m, aes(d, u)) + geom_tile(aes(fill = rescale), colour = "white") + scale_fill_gradient(low = "white",  high = "steelblue")
    return(list(error_table = err.table, best_combination = err.table[which.min(err.table[,3]),], heatmap=p))
  } 
  
  
  else if(all(dim(Beta)==c(r,p))) {
    for(u in 2:max.u){
      d=0
      while(d<max.d & d<u){
        d = d+1
        cnt = cnt+1
        est.Betas = red.rank.envelope(X,Y,u,d)
        err = norm(est.Betas-Beta,type = "F")
        err.table[cnt,1]=d; err.table[cnt,2]=u; err.table[cnt,3]=err
      }
    }
    err.table = as.data.frame(err.table)
    names(err.table) = c('d','u','Error_in_Beta')
    m <- ddply(err.table, .(d), transform, rescale = scales::rescale(Error_in_Beta))
    p <- ggplot(m, aes(d, u)) + geom_tile(aes(fill = rescale), colour = "white") + scale_fill_gradient(low = "white",  high = "steelblue")
    return(list(error_table = err.table, best_combination = err.table[which.min(err.table[,3]),], heatmap=p))
  }
  else{
    print("Error: Invalid Beta Matrix!")}
}


X = as.matrix(read.table('/Users/shnlee/Desktop/615final/X.txt'))
Y = as.matrix(read.table('/Users/shnlee/Desktop/615final/Y.txt'))
Beta = as.matrix(read.table('/Users/shnlee/Desktop/615final/Beta.txt'))

a = rrenv.choose_du(X=X, Y=Y, Beta=Beta)
b = rrenv.choose_du(X=X, Y=Y, Beta=NULL)
c = rrenv.choose_du(X=X, Y=Y, Beta=Beta2)

a$error_table
a$best_combination
a$heatmap
