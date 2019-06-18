#include "numchol.h"
void   dCopy(int,double*,double*);

static void SolFwdSnode(chfac  *sf,
                        int    snde,
                        int    f,
                        int    l,
                        double x[])
{
  int    i,t,sze,*ls,*subg=sf->subg,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead,
         *usub=sf->usub;
  double xi,*l1,*diag=sf->diag,*uval=sf->uval;

  f += subg[snde];
  l += subg[snde];

  for(i=f; i<l; ++i)
  {
    x[i] /= diag[i];
    xi    = x[i];

    ls    = usub+ujbeg[i];
    l1    = uval+uhead[i];
    sze   = l-i-1;

    for(t=0; t<sze; ++t)
      x[ls[t]] -= l1[t]*xi;
  }
} /* SolFwdSnode */

static void SolBward(int    nrow,
                     double diag[],
                     double uval[],
                     int    fir[],
                     double x[])
{
  int    i,t,sze;
  double x1,x2,rtemp,
         *x0,*l1,*l2;

  for(i=nrow; i;) {
    for(; i>1; --i) {
          -- i;
      l1   = uval+fir[i-1]+1;
      l2   = uval+fir[i  ]+0;
      sze  = nrow-i-1;
      x1   = 0.0;
      x2   = 0.0;
      x0   = x+1+i;

      for(t=0; t<sze; ++t)
      {
        rtemp = x0[t];

        x1   += l1[t]*rtemp;
        x2   += l2[t]*rtemp;
      }

      x[i]   -= x2/diag[i];
      x[i-1] -= (uval[fir[i-1]]*x[i]+x1)/diag[i-1];
    }

    for(; i;) {
          -- i;
      l1   = uval+fir[i];
      sze  = nrow-i-1;
      x1   = 0.0;
      x0   = x+1+i;

      for(t=0; t<sze; ++t)
        x1 += l1[t]*x0[t];

      x[i] -= x1/diag[i];
    }
  }
} /* SolBward */

void ChlSolveForwardPrivate(chfac  *sf,
			    double x[])
{
  int    k,s,t,sze,f,l,itemp,*ls,
         *subg=sf->subg,*ujsze=sf->ujsze,*usub=sf->usub,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead;
  double rtemp1,rtemp2,rtemp3,rtemp4,
         rtemp5,rtemp6,rtemp7,rtemp8,
         *l1,*l3,*l2,*l4,*l5,*l6,*l7,*l8,
         *uval=sf->uval;

  for(s=0; s<sf->nsnds; ++s) {
    f = subg[s];
    l = subg[s+1];

    SolFwdSnode(sf,s,0,l-f,x);

    itemp = l-f-1;
    ls    = usub+ujbeg[f]+itemp;
    sze   = ujsze[f]-itemp;
    k     = f;

    itemp = l-1;
    for(; k+7<l; k+=8) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);
      l5       = uval+uhead[k+4]+itemp-(k+4);
      l6       = uval+uhead[k+5]+itemp-(k+5);
      l7       = uval+uhead[k+6]+itemp-(k+6);
      l8       = uval+uhead[k+7]+itemp-(k+7);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];
      rtemp5   = x[k+4];
      rtemp6   = x[k+5];
      rtemp7   = x[k+6];
      rtemp8   = x[k+7];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t]
                    + rtemp5*l5[t]
                    + rtemp6*l6[t]
                    + rtemp7*l7[t]
                    + rtemp8*l8[t];
    }

    for(; k+3<l; k+=4) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t];
    }

    for(; k+1<l; k+=2) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t];
    }


    for(; k<l; ++k) {
      l1       = uval+uhead[k+0]+itemp-(k+0);

      rtemp1   = x[k+0];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t];
    }
  }
 
}

void ChlSolveBackwardPrivate(chfac  *sf,
			     double x[],
			     double b[])
{
  /* Note:  x: input, or left hand side    b: output, or solution */

  int    i,s,t,sze,f,l,*ls,
         *subg=sf->subg,*ujsze=sf->ujsze,*usub=sf->usub,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead;
  double x1,x2,*l1,*l2,rtemp1,
         *diag=sf->diag,*uval=sf->uval;


  if (sf->nsnds) {
    s = sf->nsnds - 1;
    f = subg[s];
    l = subg[s+1];
    
    dCopy(l-f,x+f,b+f);
    
    SolBward(l-f,diag+f,uval,uhead+f,b+f);
    
    s = sf->nsnds-1;
    
    for(; s>=1; --s) {
      f = subg[s-1];
      l = subg[s];
      i = l;
      
      for(; i>1+f; --i) {
	-- i;
        ls   = usub+ujbeg[i];
        l1   = uval+uhead[i-1]+1;
        l2   = uval+uhead[i  ]+0;
        sze  = ujsze[i];
        x1   = 0.0;
        x2   = 0.0;
	
        for(t=0; t<sze; ++t) {
          rtemp1 = b[ls[t]];
	  
          x1    += l1[t]*rtemp1;
          x2    += l2[t]*rtemp1;
        }
	
        b[i]   = x[i  ] -  x2  / diag[i];
        b[i-1] = x[i-1] - (x1 + uval[uhead[i-1]]*b[i]) / diag[i-1];
      }
      
      for(; i>f;) {
            -- i;
	    l1   = uval+uhead[i];
	    ls   = usub+ujbeg[i];
	    sze  = ujsze[i];
	    x1   = 0.0;
	    
	    for(t=0; t<sze; ++t)
	      x1+= l1[t]*b[ls[t]];
	    
	    b[i] = x[i] - x1/diag[i];
      }
    }
  }
   
}

/* Everything right for permuted  system */
void ChlSolveForward2(chfac  *sf,
		     double b[],
		     double x[]){
  int i,nrow=sf->nrow;
  double *sqrtdiag=sf->sqrtdiag;
  ChlSolveForwardPrivate(sf,b);
  for(i=0; i<nrow; ++i){
    x[i] = b[i]*sqrtdiag[i];   /* x[i] = b[i]*sqrt(sf->diag[i]); */
  }
}
 
void ChlSolveBackward2(chfac  *sf,
		     double b[],
		     double x[]){
  int i,nrow=sf->nrow;
  double *sqrtdiag=sf->sqrtdiag;
  for(i=0; i<nrow; ++i){
    x[i] = b[i]/(sqrtdiag[i]);   /*  x[i] = b[i]/sqrt(sf->diag[i]) ; */
  }
  ChlSolveBackwardPrivate(sf,x,b);
  memcpy(x,b,nrow*sizeof(double));
}

/* These routines together will solve an equation correctly */
void ChlSolveForward(chfac  *sf,
		     double b[],
		     double x[]){
  int i,nrow=sf->nrow,*perm=sf->perm;
  double *w=sf->rw,*sqrtdiag=sf->sqrtdiag;
  for(i=0; i<nrow; ++i)
    w[i] = b[perm[i]];
  ChlSolveForwardPrivate(sf,w);
  for(i=0; i<nrow; ++i){
    x[i] = w[i]*sqrtdiag[i];  /*   x[i] = w[i]*sqrt(sf->diag[i]); */
  }

}
 
void ChlSolveBackward(chfac  *sf,
		      double b[],
		      double x[]){
  int i,nrow=sf->nrow,*invp=sf->invp;
  double *w=sf->rw,*sqrtdiag=sf->sqrtdiag;
  for(i=0; i<nrow; ++i){
    x[i] = b[i]/sqrtdiag[i];
  }
  ChlSolveBackwardPrivate(sf,x,w);
  for(i=0; i<nrow; ++i)
    x[i] = w[invp[i]];

}
 
void ChlSolve(chfac  *sf,
              double b[],
              double x[]){
  int i,nrow=sf->nrow,*perm=sf->perm,*invp=sf->invp;
  double *rw=sf->rw;
  /*
  ChlSolveForward(sf,b,w,x);
  ChlSolveBackward(sf,w,x,b);
  */
  for(i=0; i<nrow; ++i)
    x[i] = b[perm[i]];

  ChlSolveForwardPrivate(sf,x);
  ChlSolveBackwardPrivate(sf,x,rw);

  for(i=0; i<nrow; ++i)
    x[i] = rw[invp[i]];  /* x[i] = b[i];  */
  return;
}

  
void ForwSubst(chfac  *sf,
               double b[],
               double x[])
{
  int    i,k,s,t,sze,f,l,itemp,*ls,
         *subg=sf->subg,*ujsze=sf->ujsze,*usub=sf->usub,
         *ujbeg=sf->ujbeg,*uhead=sf->uhead;
  double rtemp1,rtemp2,rtemp3,rtemp4,
         rtemp5,rtemp6,rtemp7,rtemp8,
         *l1,*l3,*l2,*l4,*l5,*l6,*l7,*l8,
         *diag=sf->diag,*uval=sf->uval;
   
  for(i=0; i<sf->nrow; ++i)
    x[i]  = b[sf->perm[i]]; 

  for(s=0; s<sf->nsnds; ++s) {
    f = subg[s];
    l = subg[s+1];

    SolFwdSnode(sf,s,0,l-f,x);

    itemp = l-f-1;
    ls    = usub+ujbeg[f]+itemp;
    sze   = ujsze[f]-itemp;
    k     = f;

    itemp = l-1;
    for(; k+7<l; k+=8) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);
      l5       = uval+uhead[k+4]+itemp-(k+4);
      l6       = uval+uhead[k+5]+itemp-(k+5);
      l7       = uval+uhead[k+6]+itemp-(k+6);
      l8       = uval+uhead[k+7]+itemp-(k+7);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];
      rtemp5   = x[k+4];
      rtemp6   = x[k+5];
      rtemp7   = x[k+6];
      rtemp8   = x[k+7];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t]
                    + rtemp5*l5[t]
                    + rtemp6*l6[t]
                    + rtemp7*l7[t]
                    + rtemp8*l8[t];
    }

    for(; k+3<l; k+=4) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);
      l3       = uval+uhead[k+2]+itemp-(k+2);
      l4       = uval+uhead[k+3]+itemp-(k+3);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];
      rtemp3   = x[k+2];
      rtemp4   = x[k+3];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t]
                    + rtemp3*l3[t]
                    + rtemp4*l4[t];
    }

    for(; k+1<l; k+=2) {
      l1       = uval+uhead[k+0]+itemp-(k+0);
      l2       = uval+uhead[k+1]+itemp-(k+1);

      rtemp1   = x[k+0];
      rtemp2   = x[k+1];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t]
                    + rtemp2*l2[t];
    }


    for(; k<l; ++k) {
      l1       = uval+uhead[k+0]+itemp-(k+0);

      rtemp1   = x[k+0];

      for(t=0; t<sze; ++t)
        x[ls[t]] -=   rtemp1*l1[t];
    }
  }
  
  for (i=0; i<sf->nrow; i++){
    x[i] = x[i] * sqrt( fabs(diag[i]) );
    }

} /* ForwSubst */



static void mulSnod(chfac  *sf,
                    int    snde,
                    int    f,
                    int    l,
                    double *b,
                    double *x)
{
  int    i,t,sze,*ls,*subg,*ujbeg,*uhead,*usub;
  double xi,*l1,*diag,*uval;

  subg =sf->subg;
  ujbeg=sf->ujbeg;
  uhead=sf->uhead;
  usub =sf->usub;
  diag =sf->diag;
  uval =sf->uval;
  
  f += subg[snde];
  l += subg[snde];

  for(i=f; i<l; ++i) {
    xi   =b[i];
    ls   =usub+ujbeg[i];
    l1   =uval+uhead[i];
    sze  =l-i-1;
    x[i]+=xi*diag[i];
    for(t=0; t<sze; ++t)
      x[ls[t]]+=l1[t]*xi;
  }
} /* mulSnod */

void GetUhat(chfac  *sf,
             double *b,
             double *x)
     /* If S = L L^T, then b = L x  */ 
{
  int    i,k,n,s,t,sze,f,l,itemp,*ls,
         *subg,*ujsze,*usub,*ujbeg,*uhead;
  double rtemp1,rtemp2,rtemp3,rtemp4,
         rtemp5,rtemp6,rtemp7,rtemp8,
         *l1,*l3,*l2,*l4,*l5,*l6,*l7,*l8,
         *diag,*uval;
  
  n    =sf->nrow; 
  subg =sf->subg;
  ujsze=sf->ujsze;
  usub =sf->usub;
  ujbeg=sf->ujbeg;
  uhead=sf->uhead;
  diag =sf->diag;
  uval =sf->uval;
  
  for (i=0; i<n; i++) {
    if (diag[i]>0)
      x[i]=b[i]/sqrt(diag[i]);
    else x[i]=b[i]/sqrt(-diag[i]);
    b[i]=0.0;
  }
  
  for (s=0; s<sf->nsnds; s++) {
    f=subg[s];
    l=subg[s+1];
    
    mulSnod(sf,s,0,l-f,x,b);
    
    itemp=l-f-1;  
    ls   =usub+ujbeg[f]+itemp;
    sze  =ujsze[f]-itemp;
    k    =f;
    
    itemp=l-1;
    for(; k+7<l; k+=8) {
      l1    =uval+uhead[k+0]+itemp-(k+0);
      l2    =uval+uhead[k+1]+itemp-(k+1);
      l3    =uval+uhead[k+2]+itemp-(k+2);
      l4    =uval+uhead[k+3]+itemp-(k+3);
      l5    =uval+uhead[k+4]+itemp-(k+4);
      l6    =uval+uhead[k+5]+itemp-(k+5);
      l7    =uval+uhead[k+6]+itemp-(k+6);
      l8    =uval+uhead[k+7]+itemp-(k+7);
        
      rtemp1=x[k+0];
      rtemp2=x[k+1];
      rtemp3=x[k+2];
      rtemp4=x[k+3];
      rtemp5=x[k+4];
      rtemp6=x[k+5];
      rtemp7=x[k+6];
      rtemp8=x[k+7];
       
      for(t=0; t<sze; ++t)
        b[ls[t]]+= rtemp1*l1[t]
                  +rtemp2*l2[t]
                  +rtemp3*l3[t]
                  +rtemp4*l4[t]
                  +rtemp5*l5[t]
                  +rtemp6*l6[t]
                  +rtemp7*l7[t]
                  +rtemp8*l8[t];
    }
    

    for(; k+3<l; k+=4) {
      l1    =uval+uhead[k+0]+itemp-(k+0);
      l2    =uval+uhead[k+1]+itemp-(k+1);
      l3    =uval+uhead[k+2]+itemp-(k+2);
      l4    =uval+uhead[k+3]+itemp-(k+3);
       
      rtemp1=x[k+0];
      rtemp2=x[k+1];
      rtemp3=x[k+2];
      rtemp4=x[k+3];

      for(t=0; t<sze; ++t)
        b[ls[t]]+= rtemp1*l1[t]
                  +rtemp2*l2[t]
                  +rtemp3*l3[t]
                  +rtemp4*l4[t];
    }

    for(; k+1<l; k+=2) {
      l1    =uval+uhead[k+0]+itemp-(k+0);
      l2    =uval+uhead[k+1]+itemp-(k+1);

      rtemp1=x[k+0];
      rtemp2=x[k+1];

      for(t=0; t<sze; ++t)
        b[ls[t]]+= rtemp1*l1[t]
                  +rtemp2*l2[t];
    }


    for(; k<l; ++k) {
      l1    =uval+uhead[k+0]+itemp-(k+0);
      
      rtemp1=x[k+0];

      for(t=0; t<sze; ++t)
        b[ls[t]]+=rtemp1*l1[t];
    }
  }
  
  for (i=0; i<n; i++)
    x[ sf->invp[i] ]=b[i];
} /* GetUhat */




int MatSolve4(chfac*sf, double *b, double *x,int n){

  memcpy(x,b,n*sizeof(double));
  ChlSolve(sf, b, x);
  
  return 0;
}

int Mat4GetDiagonal(chfac*sf, double *b,int n){

  int i,*invp=sf->invp;
  double *diag=sf->diag;

  for (i=0; i<n; i++, invp++){
    b[i]=diag[*invp];
  }
  return 0;
}

int Mat4SetDiagonal(chfac*sf, double *b,int n){

  int i,*invp=sf->invp;
  double *diag=sf->diag;
  for (i=0; i<n; i++){
    diag[invp[i]]=b[i];
  }
  return 0;
}

int Mat4AddDiagonal(chfac*sf, double *b,int n){

  int i,*invp=sf->invp;
  double *diag=sf->diag;
  for (i=0; i<n; i++){
    diag[invp[i]]+=b[i];
  }
  return 0;
}

int MatAddDiagonalElement(chfac*sf, int row, double dd){

  int *invp=sf->invp;
  double *diag=sf->diag;
  diag[invp[row]]+=dd;
  return 0;
}


int MatMult4(chfac *sf, double *x, double *y, int n){

  int i,j,*invp=sf->invp,*perm=sf->perm;
  int *usub=sf->usub,*ujbeg=sf->ujbeg,*uhead=sf->uhead, *ujsze=sf->ujsze;
  int *iptr,k1,k2;
  double dd,*sval,*diag=sf->diag,*uval=sf->uval;

  for (i=0; i<n; i++){
    y[i] = diag[invp[i]] * x[i];
  }
  for (i=0; i<n; ++i){
   
    iptr=usub + ujbeg[i];
    sval=uval + uhead[i];
    k1=perm[i];
    for (j=0; j<ujsze[i]; j++){
      dd=sval[j];
      if (fabs(dd)> 1e-15){
	k2=perm[iptr[j]];	
	y[k1] += dd * x[k2];
	y[k2] += dd * x[k1];
      }
    }
  }
  return 0;
}


static void setXYind2(int     nnz,
                     double* y,
                     double* x,
                     int*    s,
                     int*    invp)
{
  int i;
  
  for(i=0; i<nnz; ++i) {
    x[i]=y[ invp[ s[i] ] ];
    y[ invp[s[i]] ]=0.0;
  }
} /* setXYind */

static void setColi(chfac*  cl,
                    int     i,
                    double* ai)
{
  setXYind2(cl->ujsze[i],ai,
           cl->uval+cl->uhead[i],
           cl->usub+cl->ujbeg[i],
	   cl->perm);
} /* setColi */

static void setXYind2add(int     nnz,
			 double dd,
			 double* y,
			 double* x,
			 int*    s,
			 int*    invp)
{
  int i;
  
  for(i=0; i<nnz; ++i) {
    x[i]+=dd*y[ invp[ s[i] ] ];
    y[ invp[s[i]] ]=0.0;
  }
} /* setXYind */

static void setColi2(chfac*  cl,
		     int     i,double dd,
                    double* ai)
{
  setXYind2add(cl->ujsze[i],dd,ai,
	       cl->uval+cl->uhead[i],
	       cl->usub+cl->ujbeg[i],
	       cl->perm);
} /* setColi */


static void getXYind2(int     nnz,
                     double* y,
                     double* x,
                     int*    s,
                     int*    invp)
{
  int i;
  
  for(i=0; i<nnz; ++i) {
    y[ invp[s[i]] ]=x[i];
  }
} /* setXYind */


int Mat4View(chfac *sf){

  int i,j,n=sf->nrow;
  double *v=sf->rw;
  for (i=0; i<n; i++){
    for (j=0;j<n;++j) v[j]=0.0;
    getXYind2(sf->ujsze[i],v,
	      sf->uval+sf->uhead[i],
	      sf->usub+sf->ujbeg[i],
	      sf->perm);
    v[i]=sf->diag[sf->invp[i]];
    printf("Row %d, ",i);
    for (j=0;j<n;j++){
      if (v[j]!=0) printf(" %d: %4.4e ",j,v[j]);
    }
    printf("\n");
  }

  return 0;

}

int MatZeroEntries4(chfac *sf){

  int i,n=sf->n;
  double *rw=sf->rw;
  memset((void*)(sf->diag),0,n*sizeof(double));
  memset((void*)(rw),0,n*sizeof(double));
  
  for (i=0; i<n; i++){
    setColi(sf,i,rw);
  }

  return 0;
}



int MatSetColumn4(chfac *cl, double *val, int col){
  
  int pcol=cl->invp[col];

  cl->diag[pcol]=val[col];
  val[col]=0.0;
  setColi(cl,pcol,val);

  return 0;
} /* SetColumn */

int MatAddColumn4(chfac *cl, double dd, double *val, int col){
  
  int pcol=cl->invp[col];

  cl->diag[pcol]+=dd*val[col];
  val[col]=0.0;
  setColi2(cl,pcol,dd,val);

  return 0;
} /* SetColumn */


int MatSetValue4(chfac *cl, int row,int col,double val, int setmode){
  
  int i;
  double* x=cl->uval+cl->uhead[col];
  int*    s=cl->usub+cl->ujbeg[col];
  int nnz=cl->ujsze[col];  
  int insertmode=1,addmode=2;

  if (row<0 || col<0 || row>=cl->n || col>=cl->n){
    printf("CHol set Value error: Row: %d, COl: %d \n",row,col);
    return 1;
  }

  if (setmode==insertmode&&row==col){
    cl->diag[cl->invp[col]]=val;
  } else if (setmode==addmode&&row==col) {
    cl->diag[cl->invp[col]]+=val;
  } else if (setmode==insertmode){
    for(i=0; i<nnz; ++i) {
      if (s[i]==row){
	x[i]=val;
      }
    }
  } else if (setmode==addmode){
    for(i=0; i<nnz; ++i) {
      if (s[i]==row){
	x[i]+=val;
      }
    }
  } else {
    return 1;
  }

  return 0;
} /* SetValue */


int Mat4DiagonalShift(chfac*sf, double shift){
  int i,n=sf->nrow;
  double *diag=sf->diag;
  for (i=0; i<n; i++){
    diag[i] += shift;
  }
  return 0;
}

int Mat4LogDet(chfac*sf, double *dd){
  int i,n=sf->nrow;
  double *diag=sf->diag,ddd=0;
  for (i=0; i<n; i++){
    if (diag[i]<=0) return 1;
    ddd+=log(diag[i]);
  }
  *dd=ddd;
  return 0;
}

