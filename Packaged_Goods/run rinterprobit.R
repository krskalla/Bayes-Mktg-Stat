#
# run interdependent function on real data
#

source("c:/userdata/per/res/bayes book/R_package/bayes_book_r_functions.R")
dyn.load("c:/userdata/per/res/bayes book/c_code/bayesm.dll")
source("rinterprobit.R")


data=read.table("Interdependent.dat",header=TRUE)
data$price=data$price/1000
data$option=data$option/100
data$age=data$age/10
data$income=data$income/10000
data$ethnic=data$ethnic/100
data$education=data$education/100
long=data$long+100
lat=data$lat-30
X=cbind(data[,5:10],lat,long)
y=data$y

k=ncol(X)+1
n=nrow(X)
iota=array(0,n)+1
X=cbind(iota,X)
X=as.matrix(X)
y=as.matrix(y)

# construct geographic group W using zip codes
W=array(0,dim=c(n,n))
for (i in 1:n) { for (j in 1:n) {
 if(data$zip[i]==data$zip[j]) {W[i,j]=1} else {W[i,j]=0}
  } }
rsum=W%*%iota
W=W/as.vector(rsum)

Data=list(y=y,X=X,W=W)
Prior=list(betabar=rep(0,ncol(X)),A=diag(rep(.01,ncol(X))),s0=5,q0=10)
Mcmc=list(R=100,keep=1)

out=rinterprobit(Data,Prior,Mcmc)






