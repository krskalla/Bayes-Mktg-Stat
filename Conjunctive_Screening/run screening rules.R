# read data and make call to run for screening rules
# nhh = number of respondents (households)
# nset = number of choice tasks
# nsize = number of alternatives per choice task
# nxvar = number of attribute levels, dim(beta)
# natvar = number of attributes
# ntheta = total number of grid points needed for discrete attributes
# the file "camchoice.txt" contains the dependent variable (y=1 or 0)
# the file "camdesign.txt" contains the design matrix (dummy variable coding)
# the file "camatt3.dat" contains the design matrix (levels coding)
# the file "tindex.txt" contains the theta index for each attribute
#
#  note you must change path in the file names below to point to the correct
#  location
#

source("rScreen.R",local=TRUE)
dyn.load("screen.dll")


y=read.table("camchoice.txt",header=FALSE)
x=read.table("camdesign.txt",header=FALSE)
inp2=read.table("camatt3.txt",header=FALSE)
ind <- matrix(scan("tindex.txt",0), ncol=3, byrow=TRUE)

nhh=302
nset=14
nsize=7
nxvar=18
natvar=11
ntheta=36

y = array(t(y),dim=c(nsize,nset,nhh))
X = array(t(x),dim=c(nxvar,nsize,nset,nhh))
xatt=array(t(inp2),dim=c(natvar,nsize,nset,nhh))

Data=list(y=y,X=X,xatt=xatt,ind=ind,nhh=nhh,nset=nset,
          nsize=nsize,nxvar=nxvar,natvar=natvar,ntheta=ntheta)

nu=nxvar+5
Prior=list(nu=nu,V0=nu*diag(rep(1,nxvar)),
      betabarbar=as.vector(rep(0,nxvar)),
      Abeta=.01*diag(rep(1,nxvar)))
Mcmc=list(R=1000,keep=1)

out=rScreen(Data,Prior,Mcmc)


