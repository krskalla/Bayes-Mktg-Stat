rinterprobit=
function(Data,Prior,Mcmc){
#
# purpose: run interdependent preference model
#
# Arguments:
#  Data contains X (n x nxvar) and y (0,1) n vector
#       and W which is n x n matrix of interdependencies
#  Prior contains
#     betabar, A
#     beta ~ N(betabar,A-1)
#  Mcmc is a list
#      R number of McMc draws
#      keep every keepth draw
#
# initialize storage for draWs
X=Data$X
y=Data$y
W=Data$W

betabar=Prior$betabar
A=Prior$A
s0=Prior$s0
q0=Prior$q0

R=Mcmc$R
keep=Mcmc$keep

nxvar=ncol(X)
n=nrow(X)

oldbeta=matrix(double(nxvar),ncol=1)
oldtheta=matrix(double(n),ncol=1)
oldz=matrix(double(n),ncol=1) 
betadraw=matrix(double(floor(R/keep)*nxvar),ncol=nxvar)
thetadraw=matrix(double(floor(R/keep)*n),ncol=n)
sigma2draw=double(floor(R/keep))
rhodraw=double(floor(R/keep))
alphadraw=double(floor(R/keep))
llike=double(floor(R/keep))

oldrho=0.5
oldsigma2=4

lambda=eigen(W,only.values=TRUE)$values
lmin=min(lambda)
lmax=max(lambda)

# truncation points for draWs
        a=ifelse(y == 0,-100, 0)
        b=ifelse(y == 0, 0, 100)
# statistics for data augmentation
sigma=c(rep(1,n))
root=chol(chol2inv(chol(crossprod(X)+A)))
Abetabar=crossprod(A,betabar)


itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min)",fill=TRUE)
flush.console()
for (j in 1:R) {

# draw latent z
mu=X%*%oldbeta+oldtheta
oldz=rtrun(mu,sigma,a,b)

# draw beta
oldbeta=breg(oldz-oldtheta,X,betabar,A)

# draw theta
B=diag(n)-oldrho*W
omega=chol2inv(chol(diag(n)+crossprod(B)/oldsigma2))
nu=omega%*%(oldz-X%*%oldbeta)
oldtheta=nu+t(chol(omega))%*%rnorm(n)

# draw oldsigma2
aa=.5*(s0+n)
bb=1/(.5*(q0+crossprod(B%*%oldtheta)))
oldsigma2=1/rgamma(1,shape=aa,scale=bb)

# draw oldrho
newrho=rnorm(1,mean=oldrho,sd=.1)
if(newrho > 1./lmin & newrho < 1./lmax) {
Aold=crossprod(diag(n)-oldrho*W)
Anew=crossprod(diag(n)-newrho*W)
alphan=(det(Anew)^.5)*exp(-1/(2*oldsigma2)*(crossprod(oldtheta,Anew)%*%oldtheta))
alphad=(det(Aold)^.5)*exp(-1/(2*oldsigma2)*(crossprod(oldtheta,Aold)%*%oldtheta))
alpha=alphan/alphad
if(alpha=="NaN") alpha=-1
u = runif(n=1,min=0, max=1)
if(u < alpha) { oldrho = newrho} 
}

# update screen
if(j%%100==0)
          {
           ctime=proc.time()[3]
           timetoend=((ctime-itime)/j)*(R-j)
           cat(" ",j," (",round(timetoend/60,1),")",fill=TRUE)
           flush.console()
           }

mkeep=j/keep
if(mkeep*keep == (floor(mkeep)*keep))
          {betadraw[mkeep,]=oldbeta
           thetadraw[mkeep,]=oldtheta
           sigma2draw[mkeep]=oldsigma2
           rhodraw[mkeep]=oldrho
#           alphadraw[mkeep]=oldalpha
	mu=X%*%oldbeta+oldtheta
	llike[mkeep]=sum(log(pnorm(mu)*y + (1-pnorm(mu))*(1-y)+.00001))
	            }
}

ctime=proc.time()[3]
cat(" Total Time Elapsed: ",round((ctime-itime)/60,2),fill=TRUE)
list(betadraw=betadraw,sigma2draw=sigma2draw,thetadraw=thetadraw,rhodraw=rhodraw,llike=llike)
}

