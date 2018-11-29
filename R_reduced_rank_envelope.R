#R_reduced_rank_envelope.so
dyn.load("R_reduced_rank_envelope.so")
R_reduced_rank_envelopw <- function(X,Y){
   u = .C("R_reduced_rank_envelope",X,Y)$u 
   d = .C("R_reduced_rank_envelope",X,Y)$d
   Gamma = env(t(X),t(Y),u)$Gamma
   Gamma0 = env(t(X),t(Y),u)$Gamma0
   X = as.numeric(X)
   Y = as.numeric(Y)
   res = .C("R_reduced_rank_envelope",X,Y,Gamma,Gamma0)
   return(list(res$Beta,res$FitValues,res$Residuals,res$ResidualVariance))
}


