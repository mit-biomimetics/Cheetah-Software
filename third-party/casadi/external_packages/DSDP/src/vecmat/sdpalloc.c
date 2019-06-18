#include "numchol.h"

int    iSum(int,int*);
void ShutDown();
int ExitProc(int,char *);

int iAlloc(int  len,
            char *info,
	    int **r)
{
  *r=NULL;
  
  if (len) 
  {
    *r=(int*)calloc(len,sizeof(int));
    if (!(*r)){ ExitProc(OutOfSpc,info); return 1;}
  }
  return 0;
} /* iAlloc */

void iFree(int **x)
{
  int *r=*x;
  
  if (r) 
  {
    free(r);
    *x=NULL;
  }
} /* iFree */

void cFree(char **x)
{
 char *r = *x;
 /*
 printf("\n r %d &r %d",r,&r);
 getchar();
 */
   if(r)
   { free(r);
     *x = NULL;
   }

}

int dAlloc(int  len,
               char *info,
	       double **rr)
{
  double *r=NULL;

  if (len) {
    r=(double*)calloc(len,sizeof(double));
    if (!r){ ExitProc(OutOfSpc,info); return 1;}
  }
  *rr=r;
  return 0;
} /* dAlloc */

void dFree(double **x)
{
  double *r=*x;
  
  if (r) 
  {
    free(r);
    *x=NULL;
  }
} /* dFree */



int LvalAlloc(chfac *sf,
              char  *info)
{
  int ierr=0,nnz;
  
  nnz=iSum(sf->nrow,sf->ujsze);
  if ( nnz<=sf->unnz )
    return 1;
  
  sf->unnz=0;
  if (sf->uval) dFree(&sf->uval);
  ierr=dAlloc(nnz,info,&sf->uval);
  
  sf->unnz=nnz;
  if (ierr) return 1;
  return 0;
} /* LvalAlloc */

int CfcAlloc(int  maxrow,
                char *info,chfac**rr)
{
  chfac *r=NULL;
  int ierr=0;

  if (maxrow) {
    r=(chfac*)calloc(1,sizeof(chfac));
    if (!r) ExitProc(OutOfSpc,info);
    
    r->mrow =maxrow;
    r->nrow =maxrow;
      
    r->snnz =0;
    ierr=iAlloc(maxrow,info,&r->shead); if(ierr) return 1;
    ierr=iAlloc(maxrow,info,&r->ssize); if(ierr) return 1;
    r->ssub =NULL;
    ierr=dAlloc(maxrow,info,&r->diag); if(ierr) return 1;
    ierr=dAlloc(maxrow,info,&r->sqrtdiag); if(ierr) return 1;
    r->unnz =0;
    r->ujnz =0;
    ierr=iAlloc(maxrow,info,&r->ujbeg); if(ierr) return 1;
    ierr=iAlloc(maxrow,info,&r->uhead); if(ierr) return 1;
    ierr=iAlloc(maxrow,info,&r->ujsze); if(ierr) return 1;
    r->usub =NULL;
    r->uval =NULL;
    ierr=iAlloc(maxrow,info,&r->perm); if(ierr) return 1;
    ierr=iAlloc(maxrow,info,&r->invp); if(ierr) return 1;
    r->nsnds=0;
    ierr=iAlloc(maxrow+1,info,&r->subg);  if(ierr) return 1;     
    r->n=maxrow;
    r->alldense=0;
    r->tolpiv=1.0e-13; /* Standard */
    r->tolpiv=1.0e-35;
    r->cachesize   =256;
#ifdef DSDPCACHESIZE
    if (DSDPCACHESIZE>0){
      r->cachesize   = (int)DSDPCACHESIZE;
    }
#endif
    r->cacheunit   =1000;

  }
  *rr=r;
  return 0;
} /* SchlAlloc */

void CfcFree(chfac **sf)
{
  chfac *r=*sf;
  
  if (*sf) {
    iFree(&r->shead);
    iFree(&r->ssize);
    iFree(&r->ssub);
    dFree(&r->diag);
    dFree(&r->sqrtdiag);
    iFree(&r->uhead);
    iFree(&r->ujsze);
    dFree(&r->uval);
    iFree(&r->perm);
    iFree(&r->subg);
    iFree(&r->dhead);
    iFree(&r->dbeg);
    iFree(&r->dsub);
    iFree(&r->iw);
    dFree(&r->rw);
    if (r->alldense){
      r->invp=0;
      r->ujbeg=0;
      r->usub=0;
    }else{
      iFree(&r->invp);
      iFree(&r->ujbeg);
      iFree(&r->usub);
    }
    free(r);
  }
  *sf=NULL;
} /* CfcFree */

int dPtAlloc(int  n,
                  char *info,double ***rr)
{
  int    ierr,i;
  double **r;
  
  r=NULL;
  *rr=NULL;
  if (!n) return 0;
  
  r=(double **)calloc(n,sizeof(double*));
  if (!r){ 
    ExitProc(OutOfSpc,info);
    return 1;
  }
  ierr=dAlloc(n*(n-1)/2,info,&r[0]);
  if(ierr) return 1;
  for (i=1; i<n; i++)
    r[i]=r[i-1]+n-i;
  
  *rr=r;
  return 0;
} /* dPtAlloc */

void dPtFree(double ***x)
{
  double **r=*x;
  
  if (r) {
    if (r[0])
      dFree(&r[0]);
    free(r);
    *x=NULL;
  }
} /* dPtFree */

