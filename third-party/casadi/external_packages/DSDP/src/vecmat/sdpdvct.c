#include "numchol.h"
#include <stdio.h>

int ExitProc(int,char *);

void iZero(int n,
           int *x,
           int *s)
{
  int i;

  if (s) {
    for(i=0; i<n; ++i)
      x[s[i]]=0;
  }
  else
    memset(x,0,n*sizeof(int));
} /* iZero */

void iSet(int n,
          int a,
          int *x,
          int *id)
{
  int i;

  if (!id)
    for(i=0; i<n; ++i)
      x[i]=a;
  else
    for(i=0; i<n; ++i)
      x[id[i]]=a;
} /* iSet */

void iSwap(int i1,
           int i2,
           int *v)
{
  int temp;

  if (i1<0||i2<0)
    ExitProc(SysError,"index error");
  temp  = v[i1];
  v[i1] = v[i2];
  v[i2] = temp;
} /* iSwap */

void iCopy(int n,
           int *s,
           int *d)
{
  memcpy(d,s,n*sizeof(int));
} /* iCopy */

int iSum(int n,
         int *x)
{
  int i,sum=0;
  
  for (i=0; i<n; i++)
    sum+=x[i];
    
  return sum;
} /* iSum */

void dCopy(int    n,
           double *s,
           double *d)
{
  if (n) memcpy(d,s,n*sizeof(double));
} /* dCopy */

void dCat(int    n,
          int    *ix,
          double *s,
          double *d)
{
  int i;
  
  for (i=0; i<n; i++) {
    d[i]=s[ix[i]];
    s[ix[i]]=0.0;
  }
} /* dCat */

double dDot(double *x,
            double *y,
            int    n)
{
  int    i;
  double r;
  
  r=0.0;
  for (i=0; i<n; i++)
    r+=x[i]*y[i];
  
  return r;
}/* dDot */

void plusXs(int n,
            int *x,
            int *s)
{
  int i;
  
  if (!s) {
    for (i=0; i<n; i++)
      x[i]++;
  }
  else {
    for (i=0; i<n; i++)
      x[s[i]]++;
  }
} /* plusXs */
