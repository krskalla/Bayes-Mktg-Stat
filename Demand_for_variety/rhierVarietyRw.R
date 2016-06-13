rhierVarietyRw=
function(Data,Prior,Mcmc) {

varietydata=Data$varietydata
Z=Data$Z
Deltabar=Prior$Deltabar
A=Prior$A
nu=Prior$nu
V=Prior$V
s=Mcmc$s
sdelta=Mcmc$sdelta
inc.root=Mcmc$inc.root
delta0=Mcmc$delta0
lambda0=Mcmc$lambda0
Betas0=Mcmc$Betas0
r=Mcmc$r
R=Mcmc$R
keep=Mcmc$keep
#
# function to run hierarchical variety model
#
# nunits variety
#      ngoods q and price obs on each
#
#     parameters are beta(nvar x 1), delta, and lambda
#
#
# Priors:
#        beta_i ~ N(ZDelta[i,],V_beta)
#               Note:  ZDelta is the matrix Z * Delta; [i,] refers to ith row of this product!
#
#          vec(Delta) | V_beta ~ N(vec(Deltabar),A^-1 (x)  Vbeta)
#          V_beta ~ IW(nu,V)  or V_beta^-1 ~ W(nu,V^-1)
#              Delta, Deltabar are nz x nvar
#              A is nz x nz
#              Vbeta is nvar x nvar
#        
#          NOTE: if you don't have any z vars, set Z=iota (nunits x 1) 
#              
#
# Metropolis "tuning" parms
#       s is the scaling parameter for the RW inc covariance matrix.  (s^2)Var
#       R is number of draws
#
#	Vbetadraw= R+1 x nvar*nvar  matrix of Vbeta draws, each row is a draw
#	Deltadraw R+1  x nz*nvar matrix of draws of Delta, first row is initial value
#	betadraw is nunits x nvar x (R+1) array of draws of betas
#
nunits=length(varietydata)
nvar=ncol(V)
nz=ncol(Z)

#
#  create functions needed
#
# --------------------------------------------------------------------------
rw.metrop.once=
function(q,price,delta,lambda,r,oldbeta,oldlpost,s,inc.root,betabar,rootpi){ 
#
# function to exec RW metropolis for beta_is in Variety model given
#  lambda and delta
# RW increments are N(0,s^2*t(inc.root)%*%inc.root)
# prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
#	inc.root, rootpi are upper triangular
#	this means that we are using the UL decomp of Sigma^-1 for prior 
# oldbeta is the current
     stay=0
     betac=oldbeta + s*t(inc.root)%*%matrix(rnorm(length(betabar)),ncol=1)
     oldlpost=oldlpost+lndMvn(oldbeta,betabar,rootpi)
     clpost=llVarietyScale(c(betac,delta),lambda,X=q,P=price,r)+lndMvn(betac,betabar,rootpi)
     ldiff=clpost-oldlpost
     alpha=min(1,exp(ldiff))
     if(alpha < 1) {unif=runif(1)} else {unif=0}
     if (unif <= alpha)
             {betadraw=betac; oldlpost=clpost}
           else
             {betadraw=oldbeta; stay=1}
list(betadraw=betadraw,stay=stay,oldlpost=oldlpost)             
}
#
# ----------------------------------------------------------------------------
extract=function(A,B,nz,NZERO=TRUE){
# function to extract elements of B corresponding to non-zero(def) A elements
#   if NZERO=FALSE, extracts elements of B corresponding to zero A elements
# A,B must be of the same dimension
# in order for the resulting to be a matrix, there must be the same number
# of non-zero elements in every row of A
B=t(B)
if(NZERO)
  {return(t(matrix(B[t(A)!=0],nrow=nz)))}
else
  {return(t(matrix(B[t(A)==0],nrow=(ncol(A)-nz))))}
}
#
# -------------------------------------------------------------------------------

llVarietyScale = function(zeta,lambda,X,P,r) {
#
# likelihood for variety model with beta_1=0
#
# revision history:
#   P. Rossi 4/05
#   changed Jacobian 4/17/05
#   set beta_1 =0 to achieve ident 4/18/05
#   added lambda scaling argument on 4/19/05
# 
# arguments:
#   zeta is c(beta,delta)
#     beta is a parameter vector with element log(psi_j * alpha_j), j=2,...nbrand!
#   lambda scale cov matrix = lambda*Sigma
#   X is a matrix of observed quantities.  Each row is an observation dim(X)=c(nobs,nbrand)
#   P is a matrix of prices
#   delta is a vector of the reparameterized satiation parameter with element = alpha_j - 1 
#   r is number of GHK draws
#
# output: 
#   value of log-likelihood
#
#
ngoods=ncol(X)
beta=zeta[1:(ngoods-1)]; delta=zeta[ngoods:length(zeta)]
beta=c(0,beta)
iota=matrix(1,nrow=nrow(X))
V=iota%*%t(beta)+(iota%*%t(delta))*log((X+1)) - log(P)
J=(iota%*%t(delta))/(X+1)
C=cbind(matrix(1,nrow=ngoods-1),-diag(ngoods-1)) 

llike=0
nz=rowSums(X != 0)
#
# find those observations with only one non-zero demand
#
if(sum(nz==1) !=0)
{
Xone=X[nz==1,]
Vone=V[nz==1,]
dim(Xone)=c(sum(nz==1),ngoods)
dim(Vone)=c(sum(nz==1),ngoods)
VV=cbind(extract(Xone,Vone,1),extract(Xone,Vone,1,NZERO=FALSE))
ht=C%*%t(VV)
nzero=nrow(ht)
trunpt=as.vector(ht)
Sigma=lambda*(matrix(1,nrow=nzero,ncol=nzero)+diag(nzero))
Lj=t(chol(Sigma))  
above=rep(1,nzero)
llike=llike+sum(log(ghkvec(Lj,trunpt,above,r)))

remainobs=nz!=1
nremainobs=sum(remainobs)
X=X[remainobs,]
dim(X)=c(nremainobs,ngoods)
V=V[remainobs,]
dim(V)=c(nremainobs,ngoods)
J=J[remainobs,]
dim(J)=c(nremainobs,ngoods)
P=P[remainobs,]
dim(P)=c(nremainobs,ngoods)
nz=nz[remainobs]
}
#
# handle observations with two non-zero demands
#
if(sum(nz==2) !=0) 
{
Xtwo=X[nz==2,]
Vtwo=V[nz==2,]
Ptwo=P[nz==2,]
dim(Xtwo)=c(sum(nz==2),ngoods)
dim(Vtwo)=c(sum(nz==2),ngoods)
dim(Ptwo)=c(sum(nz==2),ngoods)
Delta=matrix(1,nrow=nrow(Xtwo))%*%t(delta)
VV=cbind(extract(Xtwo,Vtwo,2),extract(Xtwo,Vtwo,2,NZERO=FALSE))
Pord=cbind(extract(Xtwo,Ptwo,2),extract(Xtwo,Ptwo,2,NZERO=FALSE))
Dord=cbind(extract(Xtwo,Delta,2),extract(Xtwo,Delta,2,NZERO=FALSE))
Xord=cbind(extract(Xtwo,Xtwo,2),extract(Xtwo,Xtwo,2,NZERO=FALSE))
h=VV%*%t(C)
#
# compute density term for llike
#
ha=h[,1]
rooti=matrix((1/sqrt(lambda*2)),ncol=1,nrow=1)
llvec=sapply(ha,lndMvn,mu=0,rooti=rooti)
llike=llike+sum(llvec)
#
# compute Jacobian term
#
Jac=(Pord[,2]/Pord[,1])*(Dord[,1]/(1+Xord[,1])) + Dord[,2]/(1+Xord[,2])
llike=llike+sum(log(abs(Jac)))
#
# compute corner prob
#
a=1
b=ngoods-2
omega_aa=lambda*matrix(2,ncol=a,nrow=a)
omega_ab=lambda*matrix(1,nrow=a,ncol=b)
omega_ba=t(omega_ab)
omega_bb=lambda*(matrix(1,nrow=b,ncol=b)+diag(b))
omega_aa_inv=1/omega_aa
Sigma=omega_bb-omega_ba%*%omega_aa_inv%*%omega_ab 
ha=omega_aa_inv*ha
mu=omega_ab%x%ha
hb=h[,2:ncol(h)]
trunpt=hb-mu
trunpt=as.vector(t(trunpt))
L=t(chol(Sigma))
above=rep(1,b)
llike=llike+sum(log(ghkvec(L,trunpt,above,r)))

remainobs=nz >= 3
nremainobs=sum(remainobs)
X=X[remainobs,]
dim(X)=c(nremainobs,ngoods)
V=V[remainobs,]
dim(V)=c(nremainobs,ngoods)
J=J[remainobs,]
dim(J)=c(nremainobs,ngoods)
P=P[remainobs,]
dim(P)=c(nremainobs,ngoods)
}

if(nremainobs > 0)
{
#
# loop over observations and accumlate likelihood
for (i in 1:nrow(X)) {
  Vp=as.vector(V[i,X[i,]!=0])
  V0=as.vector(V[i,X[i,]==0])
  VV=matrix(as.numeric(c(Vp,V0)),ncol=1)
  h=C%*%VV

  a=length(Vp)-1     # number of non-zero demands -1
  b=length(V0)       # number of zeroes


  # term for non-zero demands
    indnz=which(X[i,]!=0)
    indfnz=indnz[1]
    Ji=matrix(0,ncol=a,nrow=a)
    diag(Ji)=as.vector(J[i,indnz[2:length(indnz)]]) 
    Ji=Ji +(delta[indfnz]/(1+X[i,indfnz]))*
       t(matrix((P[i,indnz[2:length(indnz)]]/P[i,indfnz]),nrow=a,ncol=a))
    jac=determinant(Ji,logarithm=TRUE)$modulus

    llike=llike+jac   #Jacobian

    ha=h[1:a]
    mu_a=matrix(0,nrow=a,ncol=1)
    omega_aa=lambda*(matrix(1,nrow=a,ncol=a)+diag(a))
    rooti=backsolve(chol(omega_aa),diag(ncol(omega_aa)))
    llike=llike+lndMvn(ha,mu_a,rooti)                 # density
   

  if(b>0) 
  {  # term for zeroes or corners
       # here we have to compute conditonal cov matrix 
       omega_ab=lambda*matrix(1,nrow=a,ncol=b)
       omega_ba=t(omega_ab)
       omega_bb=lambda*(matrix(1,nrow=b,ncol=b)+diag(b))
       omega_aa_inv=chol2inv(chol(omega_aa))
       Sigma=omega_bb-omega_ba%*%omega_aa_inv%*%omega_ab 
       mu=omega_ba%*%omega_aa_inv%*%as.vector(ha) 

       hb=h[(a+1):(a+b)]-mu 
       Lj=t(chol(Sigma))  
       above=rep(1,b)
       llike=llike+sum(log(ghkvec(Lj,hb,above,r)))  # corner prob mass
  }
}
}
return(llike)
}

#
#-------------------------------------------------------------------------------------
#

#
#
#  set up fixed parms for the draw of Vbeta,Delta
#
#  note: in the notation of the MVR  Y =    X      B  
#                                  n x m  n x k  k x m
#                           "n" = nunits
#                           "m" = nvar
#                           "k" = nz
#			general model: Beta = Z Delta + U
#
#
#  allocate space for the draws and set initial values of Vbeta and Delta
#
Vbetadraw=matrix(double((R/keep)*nvar*nvar),ncol=nvar*nvar)
Deltadraw=matrix(double((R/keep)*nz*nvar),ncol=nz*nvar)
deltadraw=matrix(double((R/keep)*(nvar+1)),ncol=(nvar+1))
betadraw=array(double((R/keep)*nunits*nvar),dim=c(nunits,nvar,R/keep))
oldlpost=double(nunits)
newlpost=double(nunits)
olddpost=double(nunits)
newdpost=double(nunits)
reject = array(0,dim=c(R/keep))
loglike = array(0,dim=c(R/keep))
#
#
#  set up intial values 
#
lambda=lambda0
delta=delta0
Betas=Betas0
Vbeta=1*diag(nvar)
Delta=matrix(c(rep(0,nz*nvar)),ncol=nvar)
#
#
#  intialize betabar
#
      betabar=Z%*%Delta
#
#  intialize the oldlpost vector
#
rootp=chol(Vbeta)
rootpi=backsolve(rootp,diag(ncol(rootp)))
for (unit in 1:nunits) 
{
  oldlpost[unit]= llVarietyScale(c(Betas[unit,],delta),lambda,X=varietydata[[unit]]$q, 
                    P=varietydata[[unit]]$price,r) 
}
       
#
#	start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration | Time to End (min)",fill=TRUE)
flush.console()
for(rep in 1:R)
{
#
#
  betabar=Z%*%matrix(Delta,ncol=nvar)
#
#       loop over all unit equations
#
       rootp=chol(Vbeta)
       rootpi=backsolve(rootp,diag(nvar))
	rej=0
  for (unit in 1:nunits) 
      {
#      draw beta_i | Vbeta,betabar
       metropout=rw.metrop.once(varietydata[[unit]]$q,varietydata[[unit]]$price,delta,lambda,r,
                                    Betas[unit,],oldlpost[unit],s,inc.root,betabar[unit,],rootpi)
       Betas[unit,]=metropout$betadraw
       oldlpost[unit]=metropout$oldlpost
	rej=rej+metropout$stay
      }
#
#          draw Vbeta, Delta | {beta_i}
#
   regout=rmultireg(Betas,Z,Deltabar,A,nu,V)
   Vbeta=regout$Sigma
   Delta=regout$B
#
#     draw satiation parms (delta) | {beta_i}
#
      deltanew=delta+sdelta*rnorm(length(delta))
   for (unit in 1:nunits)
       {
      olddpost[unit]=llVarietyScale(c(Betas[unit,],delta),lambda,X=varietydata[[unit]]$q,
                                 P=varietydata[[unit]]$price,r)
      newdpost[unit]=llVarietyScale(c(Betas[unit,],deltanew),lambda,X=varietydata[[unit]]$q,
                                 P=varietydata[[unit]]$price,r)
       }
      alpha=min(1,exp(sum(newdpost)-sum(olddpost)))
#     assume U(-5,+5) prior on delta so values outside this range have a zero likelihood
      if(max(deltanew) >= 0 | min(deltanew) <= -1) {alpha=0}
      if(alpha < 1) {unif=runif(1)} else {unif=0}
      if(unif <= alpha)
             {delta=deltanew; oldlpost=newdpost} else {oldlpost=olddpost}
  if(rep%%100 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/rep)*(R-rep)
    cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
    fsh()}


  if(rep%%keep == 0) 
    {mkeep=rep/keep
     Vbetadraw[mkeep,]=Vbeta
     Deltadraw[mkeep,]=Delta
     deltadraw[mkeep,]=delta
     reject[mkeep]=rej/nunits
     loglike[mkeep]=sum(olddpost)
     betadraw[,,mkeep]=Betas}

}
ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

list(Vbetadraw=Vbetadraw,Deltadraw=Deltadraw,betadraw=betadraw,deltadraw=deltadraw,reject=reject,loglike=loglike)
}
