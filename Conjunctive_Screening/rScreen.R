rScreen =
function(Data,Prior,Mcmc){
# revision history:
#   created 12/04 by Allenby and McCulloch
#
# purpose: estimate MNP model with conjunctive screening rules
#      and heterogeneity
#
# Arguments:
#   Data contains a list of (X, y, xatt, ind, nhh, nset, nsize,nxvar,ntheta)
#      X : design matrix
#          dim(X) = (nxvar,nsize,nset,nhh)
#          nxvar= number of attribute level variables (discrete attributes are dummied)
#          nsize= number of choice alternatives
#          nset= number of choice tasks
#          nhh= number of respondents
#          e.g. respondent #1, faces 14 choice tasks, each with 7 alternatives and 
#          11 attributes coded into 18 variables
#      y : matrix of indicators for responses
#          dim(y) = (nsize,nset,nhh)
#          y_i,j,k = 1 if alternative i is chosen in choice task j for respondent k
#      xatt : matrix of attribute level values 
#          dim(xatt) = (natvar,nsize,nset,nhh)
#          natvar= number of attribute variables (e.g. 11)
#          e.g. if attribute variable 1 takes on three levels (codes as 0,1,2) then xatt(1,j,k,l) 
#          takes on values (0,1,2)  for further explanation see table 5 in Gilbride and Allenby(04)
#      ind: an array of dim(natvar-1,3)
#          index into theta vector
#          ind[1,] is the discrete attribute variables
#          ind[i,2] is the beginning index for attribute var i
#          ind[i,3] is the end index for attribute var i
#          col1: discrete attribute number
#	ntheta = number of cut-offs that are possible
#   Prior is a list of (nu,V0,betabar,Abeta)
#         beta ~ N(mu,Vbeta)
#         mu ~ N(betabar,Abeta^-1_
#         Vbeta ~ IW(nu,V0)
#       note: dirchlet prior is set with a=6
#             normal prior on gamma for lnprice, N(mu,sig)
#                    mu ~ N(0,100); sig ~ IGamma(10,10)
#   Mcmc is a list of (R,keep)
#
# Output:
#  betabardraw: R/keep by nxvar
#  Vbetadraw:  R/keep by nxvar**2
#  betadraw: nhh x nxvar x R/keep
#  thetadraw: R/keep by ntheta
#  gammabardraw: R/keep by 1
#  gammasigdraw: R/keep by 1
# 
#
# Note:
#   program set up for natvar-1 discrete variables, 1 continuous variable (-ln price)
#   the -ln price variable is assumed to be the last variable.  
#   The price variable is coded as "-ln price" to allow for an upper threshold

X=Data$X
y=Data$y
xatt=Data$xatt
ind=Data$ind

nhh=Data$nhh
nset=Data$nset
nsize=Data$nsize
nxvar=Data$nxvar
natvar=Data$natvar
ntheta=Data$ntheta

V0=Prior$V0
nu=Prior$nu
betabarbar=Prior$betabarbar
Abeta=Prior$Abeta

R=Mcmc$R
keep=Mcmc$keep

#
# -----------------------------------------------------------------------------
#
# define functions needed
#
cchoiceset = function(xatt,gamma) {
dd = dim(xatt)
natvar = dd[1]; nsize = dd[2]; nset = dd[3]
.C('choiceset',as.integer(natvar),as.integer(nsize),as.integer(nset),as.double(xatt),as.double(gamma),cset = integer(nsize*nset))$cset
}


ctestgam = function(xatt,y,z,gamma) {
dd = dim(xatt)
natvar = dd[1]; nsize = dd[2]; nset = dd[3]
.C('testgam',as.integer(natvar),as.integer(nsize),as.integer(nset),as.double(xatt),as.integer(y),as.double(z),as.double(gamma),test = integer(1))$test
}


cgetbounds = function(z,c,y) {
n = length(z)
temp = .C('getbounds',as.integer(n),as.double(z),as.integer(c),as.integer(y),zc = double(1),zy=double(1))
list(zc = temp$zc, zy = temp$zy)
}

cdrawz = function(V,Y,C,Z) {
dd = dim(V)
array(.C('drawz',as.integer(dd[1]),as.integer(dd[2]),as.double(V),as.integer(Y),as.integer(C),z=as.double(Z))$z,dim=dd)
}


trunl <- function(lower,mu, sig){
#Draws from truncated normal distribution
# sig is the variance
stdev <- sqrt(sig)
prob <- pnorm(lower, mean=mu, sd=stdev, lower.tail=TRUE)
nprob <- runif(n=1, min=prob, max = 1)
draw <- qnorm(nprob, mean=mu, sd=stdev)
draw
}

#
#
#
trunu <- function(upper,mu, sig){
#Draws from truncated normal distribution
# sig is the variance
stdev <- sqrt(sig)
prob <- pnorm(upper, mean=mu, sd=stdev, lower.tail=TRUE)
nprob <- runif(n=1, min=0, max = prob)
draw <- qnorm(nprob, mean=mu, sd=stdev)
draw
}

#
#
#
trund <- function(lower,upper,mu, sig){
#Draws from truncated normal distribution
# sig is the variance
stdev <- sqrt(sig)
probu <- pnorm(upper, mean=mu, sd=stdev, lower.tail=TRUE)
probl <- pnorm(lower, mean=mu, sd=stdev, lower.tail=TRUE)
nprob <- runif(n=1, min=probl, max = probu)
draw <- qnorm(nprob, mean=mu, sd=stdev)
draw
}


#Initialize vector and matrix for mean and covariance for beta's
# e.g. the distribution of heterogeneity

bmean = array(0, dim=c(nxvar))
V = as.matrix(diag(1,nxvar))


z = array(0, dim=c(nsize,nset,nhh))
Vbetadraw = matrix(double(floor(R/keep)*nxvar*nxvar),ncol=nxvar*nxvar)
betabardraw = array(0, dim=c(R/keep,nxvar))
oldbetadraw  = matrix(double(nhh*nxvar),ncol=nxvar)
betadraw = array(0,dim=c(nhh,nxvar,floor(R/keep)))
thetadraw = array(0, dim=c(floor(R/keep),ntheta))
gammabardraw =array(0, dim=c(R/keep))
gammasigdraw =array(0, dim=c(R/keep))
bbar = matrix(0, ncol=1,nrow=nxvar)
Bdbar = matrix(0, nxvar,1)
Vh = matrix(0,nxvar,nxvar)
Vi = as.matrix(diag(1,nxvar))
sigma = 1
iota=matrix(rep(1,nhh),ncol=1)


#Initialize values for gamma's and theta's
#and other parameters for conjunctive screening rules

cset = array(1, dim=c(nsize,nset,nhh))
oldtheta= matrix(0, ntheta,1)
gamma= matrix(-0.5,nhh, natvar)
gamma[ ,natvar] = -7
gind= array(0, dim=c(nhh,ntheta))
oldgammabar = -7
oldgammasig = 1


for (i in 1:(natvar-1)){
	begin = ind[i,2]
	end = ind[i,3]
	zend = end-begin+1
	n = 0

	for (j in 1:zend){
		n = begin + j - 1
		oldtheta[n] = (1/(ind[i,3]-ind[i,2]+1))
		}
	}

itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min)",fill=TRUE)
flush.console()

for(o in 1:R){
	
#	Draw B-h|B-bar, V
#	Prior is distribition of heterogeneity

	for(i in 1:nhh){

#	z|Beta, sigma, y, cset

	XX=matrix(X[,,,i],nrow=nxvar,ncol=nsize*nset)
	VV=matrix(t(oldbetadraw[i,])%*%XX,nrow=nsize,ncol=nset)
	z[,,i]=cdrawz(VV,y[,,i],cset[,,i],z[,,i])

#	Beta|z, sigma
#	have to stack up the obs to do matrix manipulations

	Xd=t(XX)
	w=matrix(z[,,i],nrow=nset*nsize,ncol=1)

	XTw = t(Xd) %*% w
	Vh = chol2inv(chol(crossprod(Xd,Xd) + Vi))
	root=chol(Vh)
	Bdbar = Vh %*% ((XTw) + Vi%*%bbar)
	oldbetadraw[i,] = Bdbar + t(root)%*%rnorm(length(Bdbar))
	}	
	
# 	Draw B-bar|B-h,V
#	betabar ~ N(betabarbar,Abeta-1)
	vv=chol2inv(chol(nhh*Vi+Abeta))
	mu=vv%*%(nhh*Vi%*%(t(oldbetadraw)%*%iota/nhh)+Abeta%*%betabarbar)
	bbar=mu + t(chol(vv))%*%rnorm(nxvar)

#	Draw V|B-h,B-bar
#	Prior is IW(nu,V0)
		S=crossprod(oldbetadraw-t(array(bbar,dim=c(nxvar,nhh))))
		W=rwishart(nhh+nu,chol2inv(chol(V0+S)))
		V=W$IW
		Vi=W$W

#	Draw Gamma-h|Theta, z, y
#	Likelihood is indicator function
#	Prior is distribution of heterogeneity

	reject = 0

	for (i in 1:nhh){

	for (j in 1:(natvar-1)){

	size = ind[j,3] - ind[j,2] + 1
	begin = ind[j,2]
	end = ind[j,3]
	gammanew = array(0,dim=c(natvar))
	checkgam = array(0,dim=c(size))
	choice = array(0,dim=c(size))	

#	Determine which values of gamma are ok

		for (k in 1:size){
		gammanew = gamma[i,]
		gammanew[j] = (-.5)+(k-1)
		checkgam[k] = ctestgam(xatt[,,,i],y[,,i],z[,,i],gammanew)
		}

#	Choose among the gamma's that work using griddy gibbs
	
	thetat= oldtheta[begin:end]
	thetat = thetat*checkgam
	denom = sum(thetat)

	prob = array(0, dim=c(size))
	cprob = array(0, dim=c(size))

	prob = thetat / denom

	cprob[1] = prob[1]
	for(k in 2:size){
	cprob[k] = cprob[k-1] + prob[k]
	}

	a = runif(1, min=0, max=1)

	newgamma = 0
	for(k in 1:size){
	newgamma = newgamma+1
	if(a < cprob[k]) break
	}
	
	gamma[i,j] = (-0.5)+(newgamma-1)
	choice[newgamma] = 1
	gind[i,begin:end] = choice
	
	}

#	Now have to choose a gamma for the continuous variable

	gammanew = array(0,dim=c(natvar))
	gammanew = gamma[i,]
	gammanew[natvar] = gamma[i,natvar] + rnorm(1,mean=0, sd=0.5)

	cgam = ctestgam(xatt[,,,i],y[,,i],z[,,i],gammanew)

	znum = exp((-.5/oldgammasig)*(gammanew[natvar]-oldgammabar)^2)
	denom = exp((-.5/oldgammasig)*(gamma[i,natvar]-oldgammabar)^2)

 	zalpha = (znum/denom)     
   
	a = runif(1, min=0, max=1)

	if(zalpha > a & cgam > 0 ) gamma[i,natvar] = gammanew[natvar]
	if(zalpha < a | cgam < 1 ) reject = reject + 1

	}

#	DETERMINE CHOICE SETS|GAMMA, xatt
	
	for(i in 1:nhh){
	cset[,,i] = cchoiceset(xatt[,,,i],gamma[i,])
	}

#	DRAW THETA|GAMMA (GIND WHICH IS INDICATOR)
#	prior is dirichlet(6), use gamma dist to make draw
#	assumes natvar-1 discrete variables, 1 continuous varaible

	draw = array(0, dim=c(ntheta))

	for(i in 1:ntheta){
	cnt = sum(gind[,i])
	alpha = cnt + 6
	draw[i] = rgamma(1, shape=alpha, scale=1)
	}

	for(j in 1:(natvar-1)){	
	begin = ind[j,2]
	end = ind[j,3] 
	denom = sum(draw[begin:end])
	oldtheta[begin:end] = draw[begin:end] / denom
	}


#	DRAW GAMMABAR|GAMMASIG, GAMMA[natvar]
#	prior is n(0,100)

	gammabr = sum(gamma[,natvar])/nhh
	gammasg = oldgammasig/nhh
	gammab = (100/(100+gammasg))*gammabr
	stddev = sqrt((100*gammasg)/(100+gammasg))

	oldgammabar = rnorm(1, mean=gammab, sd=stddev)


#	DRAW GAMMASIG|GAMMABAR, GAMMA[natvar]
#	prior is inverse gamma(10,10)

	avg = array(oldgammabar, dim=c(nhh))
	sse = sum((gamma[,natvar]-avg)^2)
	nalpha = .5*(10+nhh)
	nbeta = 1/(.5*(10+sse))
	sigi = rgamma(1, shape=nalpha, scale=nbeta)
	oldgammasig = 1/sigi

# update screen
	if(o%%100==0)
          {
           ctime=proc.time()[3]
           timetoend=((ctime-itime)/o)*(R-o)
           cat(" ",o," (",round(timetoend/60,1),")",fill=TRUE)
           flush.console()
           }
	mkeep=o/keep
	if(mkeep*keep == (floor(mkeep)*keep))
          {betabardraw[mkeep,]=bbar
           Vbetadraw[mkeep,]=as.vector(V)
           betadraw[,,mkeep]=oldbetadraw
	   gammabardraw[mkeep]=oldgammabar
	   gammasigdraw[mkeep]=oldgammasig 
	   thetadraw[mkeep,]=oldtheta
          }

	
}
ctime=proc.time()[3]
cat(" Total Time Elapsed: ",round((ctime-itime)/60,2),fill=TRUE)

list(betabardraw=betabardraw, Vbetadraw=Vbetadraw, betadraw=betadraw,thetadraw=thetadraw, gammabardraw=gammabardraw,
	gammasigdraw=gammasigdraw)

}
