#R_reduced_rank_envelope.so
library(Renvlp)
library(expm)
library(readr)
setwd("/Users/Loielaine/Desktop/Good_Good_Study/2018Fall/BIOSTAT615_Statistic_Computing/project/Rpackages-Envelopes_and_ReducedRankReg/")
dyn.load("R_reduced_rank_envelope.so")
R_reduced_rank_envelope <- function(X,Y){
   #u = .C("R_reduced_rank_envelope",X,Y)$u 
   #d = .C("R_reduced_rank_envelope",X,Y)$d
   u = 30
   d = 10
   X = as.matrix(X)
   Y = as.matrix(Y)
   write_csv(data.frame(env(t(X),t(Y),u)$Gamma),"Gamma.csv",col_names = F)
   write_csv(data.frame(env(t(X),t(Y),u)$Gamma0),"Gamma0.csv",col_names = F)
   write_csv(data.frame(X),"X.csv",col_names = F)
   write_csv(data.frame(Y),"Y.csv",col_names = F)
   if(u==nrow(Y)){Gamma0 = matrix(0)}
   res = .C("R_reduced_rank_envelope")
   return(list(res$Beta_hat,res$FitValues,res$Residuals,res$ResidualVariance))
}
#run 
Y<-as.matrix(read.table("Y.txt"))#Y file from data_simulation.cpp
X<-as.matrix(read.table("X.txt"))#X file from data_simulation.cpp
res <- R_reduced_rank_envelope(X,Y)
