#BIOSTAT 615 Final Project
#Simulation of data

#Setting initial values
p<-10# Number of rows in Beta
r<-6# Number of columns in Beta
d<-3# Rank of Beta
u<-5# Rank of the Envelope
Gamma<-matrix(runif(u*r),nrow=r)
eta<-matrix(runif(u*d),nrow=u)
B<-matrix(runif(d*p),nrow=d)
Gamma<-svd(Gamma)$u
beta<-Gamma%*%eta%*%B
beta<-Gamma%*%eta%*%B/norm(beta,type="f")
norm(beta,type = "f")


x<-seq(1,10)
y<-0.4+0.5*x+rnorm(10)
plot