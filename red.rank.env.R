#Pedro Orozco del Pino
#Evaluate our simulations
library(Renvlp)
library(expm)
library(readr)
setwd("~/Box Sync/Dropbox/Michigan/Courses/Fall 2018/BIOSTAT 615/Final Prjoect")
#setwd("/Users/Loielaine/Desktop/Good_Good_Study/2018Fall/BIOSTAT615_Statistic_Computing/project/Rpackages-Envelopes_and_ReducedRankReg/")
#WARNING: To follow the articles notation individuals are keep like columns insteadas the usual notation where individuals are rows. This is why the envelope recibes t(X) and t(Y)
Y<-as.matrix(read.table("Y.txt"))#Y file from data_simulation.cpp
X<-as.matrix(read.table("X.txt"))#X file from data_simulation.cpp
Beta<-as.matrix(read.table("Betas.txt"))#Betas file from data_simulation.cpp
r<-dim(Y)[1]#number of variables as response
p<-dim(X)[1]#number of predictor
u<-3 #Dimension of the Envelope
d<-2 #Rank of the betas
red.rank.envelope<-function(X,Y,u,d){
  E<-env(t(X),t(Y),u)#Estimate the envelope given the dimension
  Gamma<-E$Gamma
  Gamma0<-E$Gamma0
  N<-dim(Y)[2]
  ifelse(u!=dim(Y)[1],Y.g<-t(Gamma)%*%as.matrix(Y),Y.g<-Y) #If dim(envelope) is equal to the initial space then keep same Y, this is traditional reduced-rank regression
  r.g<-dim(Y.g)[1]
  ########## ESTIMATE THE REDUCED RANK ENVELOPE ESTIMATOR
  Sx<-X%*%t(X)/N #Empirical variance for X
  Sxyg<-t(X%*%scale(t(Y.g), scale = FALSE))/N #Empirical covariance
  Syg<-Y.g%*%t(Y.g)/N #Empirical variance of Y
  Cyg.x<-sqrtm(solve(Syg))%*%Sxyg%*%sqrtm(solve(Sx)) #Canonical Covariance
  D.truncated<-diag(c(svd(Cyg.x)$d[1:d],rep(0,length(svd(Cyg.x)$d)-d))) #Tacking the highest d eigen values
  Cyg.x.d<-svd(Cyg.x)$u%*%D.truncated%*%t(svd(Cyg.x)$v) #Reduced Canonical Covariance 
  beta.r.e<-Gamma%*%sqrtm(Syg)%*%Cyg.x.d%*%solve(sqrtm(Sx)) #Estimated betas of the reduced rank regression
  return(beta.r.e)
}
OLS<-function(X,Y){
  r<-dim(Y)[1]
  p<-dim(X)[1]
  betas<-matrix(rep(0,r*p),ncol = p)
  for(i in 1:r)
    betas[i,]<-lm(Y[i,]~t(X))$coefficients[-1]
  return(betas)
}

#write files 
write_csv(data.frame(env(t(X),t(Y),u)$Gamma),"Gamma.csv",col_names = F)
write_csv(data.frame(red.rank.envelope(X,Y,u,d)),"Beta_R.csv",col_names = F)

print("Reduce-rank ")
norm(red.rank.envelope(X,Y,dim(Y)[1],2)-Beta,type = "F")#reduce Rank
print("Reduce-rank Envelope")
(norm(red.rank.envelope(X,Y,u,d)-Beta,type = "F"))#Reduced rank envelope
print("OLS")
(norm(OLS(X,Y)-Beta,type="F"))#OLS
print("Envelope")
(norm(env(t(X),t(Y),u)$beta-Beta,type="F"))#Envelope


