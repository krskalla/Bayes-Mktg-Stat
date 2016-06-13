#include <R.h>
#include <Rmath.h>
#include <math.h>

double rtrun(double mu, double sigma,double trunpt, int above) 
{
/*	function to draw truncated normal
		above=1 means from above b=trunpt, a=-inf
		above=0 means from below a=trunpt, b= +inf   
*/
	double FA,FB,rnd,result ;
	if (above) {
		FA=0.0; FB=pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
			}
	else {
		FB=1.0; FA=pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
		}
	
	GetRNGstate();
	rnd=unif_rand();
	result = mu + sigma*qnorm((rnd*(FB-FA)+FA),0.0,1.0,1,0);
	PutRNGstate();
	return result;
}

void choiceset(int *natvar_i, int *nsize_i , int *nset_i, double *xatt, double *gamma, int *cset)
{
   int natvar = *natvar_i; int nsize = *nsize_i; int nset = *nset_i;
   int i,j,k;
   int cst;

   for(i=0;i<nset;i++) {
      for(j=0;j<nsize;j++) {
         cst = 1;
         k=0;
         while((k<natvar) && cst) {
            if((*(xatt + i*(natvar*nsize) + j*natvar + k)) > gamma[k]) k++;
            else cst = 0;
         }
         if(j == (nsize-1)) cst=1;
         *(cset + i*nsize + j) = cst;
      }
   }
}

void testgam(int *natvar_i, int *nsize_i , int *nset_i, double *xatt, int *y, double *z,double *gamma, int *test)
{
   int natvar = *natvar_i; int nsize = *nsize_i; int nset = *nset_i;
   int i,j,k;  
   double zmax;
   int cset;
   double zz;
   int yy;

   *test = 1;
   i=0;
   while((i<nset) && (*test)) {
      k=0;
      while((*(y + k + i*nsize)) != 1) k++;
      zmax = *(z + k + i*nsize);
      j=0;
      while((j<nsize) && (*test)) {
         cset = 1;
         k=0;
         while((k<natvar) && cset) {
            if((*(xatt + i*(natvar*nsize) + j*natvar + k)) > gamma[k]) k++;
            else cset=0;
         }
         yy = *(y+j+i*nsize);
         if(!cset && yy) *test = 0;
         zz = *(z+j+i*nsize);
         if(cset && (zz>zmax)) *test = 0;
         if(j==(nsize-1)) *test = 1;
         j++;
      }
      i++;
   }
}

void getbounds(int n,double *z, int *c, int *y, double *zc, double *zy)
{
   *zc = -1000;
	int k;
	for(k=0;k<n;k++) {
	   if(y[k]) *zy = z[k];
		else {
		   if(c[k] && (z[k] > (*zc))) *zc = z[k];
		}
	}
}

void drawz(int *nsize_i, int *nset_i, double *v, int *y, int *c, double *z)
{
	int nsize = *nsize_i; int nset = *nset_i;

	int j,k;
   double *vt;
	double *zt;
	int *ct;
	int *yt;
	double zc,zy,zold;

	GetRNGstate();

   for(j=0;j<nset;j++) {

      vt = v+j*nsize;
      zt = z+j*nsize;
	   ct = c+j*nsize;
	   yt = y+j*nsize;

      getbounds(nsize,zt,ct,yt,&zc,&zy);

	   for(k=0;k<nsize;k++) {
		   if(ct[k]) {
			   if(yt[k]) {
				   zt[k] = rtrun(vt[k],1.0,zc,0);
					zy = zt[k];
				}
				else {
					zold = zt[k];
				   zt[k] = rtrun(vt[k],1.0,zy,1);
					if(zt[k]>zc) zc = zt[k];
					else {
					   if(zold == zc) getbounds(nsize,zt,ct,yt,&zc,&zy);
					}
				}
			}
			else {
			   zt[k] = vt[k] + norm_rand(); 
			}
		}
	}
	PutRNGstate();
}





