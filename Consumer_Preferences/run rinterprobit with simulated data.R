data=read.table("Simdata.dat",header=TRUE)
X=data[,3:4]
y=data[,2]

#-----------------------------------------------------
# extract data for testing 
# data simulated such that beta1=beta2=1.0, sigma2=4.0, rho=0.5
#
X=as.matrix(X[1:50,])
y=as.matrix(y[1:50])
#-----------------------------------------------------

k=ncol(X)
n=nrow(X)
iota=array(0,n)+1
# construct social group Ws
W=array(0,dim=c(n,n))
for (i in 2:n-1) {W[i,i+1]=1; W[i,i-1]=1}
W[1,2]=1; W[1,n]=1
W[n,n-1]=1; W[n,1]=1
rsum=W%*%iota
W=W/as.vector(rsum)

# initial values
Data=list(y=y,X=X,W=W)
Prior=list(betabar=rep(0,ncol(X)),A=diag(rep(.01,ncol(X))),s0=5,q0=10)
Mcmc=list(R=20000,keep=1)

out=rinterprobit(Data,Prior,Mcmc)




