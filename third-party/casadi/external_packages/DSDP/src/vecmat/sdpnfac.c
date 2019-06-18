#include "numchol.h"

static void UpdSnode(int    m,
                     int    n,
                     int    s,
                     double diaga[],
                     double *a,
                     int    fira[],
                     double diagb[],
                     double *b,
                     int    firb[])
{
  int    i,k,t,sze;
  double rtemp1, rtemp2, rtemp3, rtemp4,
         rtemp5, rtemp6, rtemp7, rtemp8,
         rtemp9, rtemp10,rtemp11,rtemp12,
         rtemp13,rtemp14,rtemp15,rtemp16,
         *a1,*a2, *a3, *a4, *a5, *a6, *a7, *a8,
         *a9,*a10,*a11,*a12,*a13,*a14,*a15,*a16,
         *b00,*b1,*b0;
    
  for(i=0; i+1<s; i+=2) {
    b00 = b+firb[i];
    b0  = b00+1;
    b1  = b+firb[i+1];
    sze = m-i-2;
    k   = 0;

    for(; k+3<n; k+=4) {
      a1          = a+(fira[k+0 ]+i);
      a2          = a+(fira[k+1 ]+i);
      a3          = a+(fira[k+2 ]+i);
      a4          = a+(fira[k+3 ]+i);

      rtemp1      = a1[0]/diaga[k+0];
      rtemp2      = a2[0]/diaga[k+1];
      rtemp3      = a3[0]/diaga[k+2];
      rtemp4      = a4[0]/diaga[k+3];

      diagb[i]   -=   rtemp1  * a1[0]
                    + rtemp2  * a2[0]
                    + rtemp3  * a3[0]
                    + rtemp4  * a4[0];

                 ++ a1;
                 ++ a2;
                 ++ a3;
                 ++ a4;

      rtemp5      = a1[0]/diaga[k+0];
      rtemp6      = a2[0]/diaga[k+1];
      rtemp7      = a3[0]/diaga[k+2];
      rtemp8      = a4[0]/diaga[k+3];

      b00[0]     -=   rtemp1 * a1[0]
                    + rtemp2 * a2[0]
                    + rtemp3 * a3[0]
                    + rtemp4 * a4[0];

      diagb[i+1] -=   rtemp5  * a1[0]
                    + rtemp6  * a2[0]
                    + rtemp7  * a3[0]
                    + rtemp8  * a4[0];

               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;

      for(t=0; t<sze; ++t) {
        b0[t] -=   rtemp1  * a1[t]
                 + rtemp2  * a2[t]
                 + rtemp3  * a3[t]
                 + rtemp4  * a4[t];

        b1[t] -=   rtemp5  * a1[t]
                 + rtemp6  * a2[t]
                 + rtemp7  * a3[t]
                 + rtemp8  * a4[t];
      }
    }

    for(; k+1<n; k+=2) {
      a1          = a+(fira[k+0 ]+i);
      a2          = a+(fira[k+1 ]+i);

      rtemp1      =    a1[0] / diaga[k+0];
      rtemp2      =    a2[0] / diaga[k+1];

      diagb[i]   -=   a1[0] * rtemp1
                    + a2[0] * rtemp2;

                 ++ a1;
                 ++ a2;

      rtemp5      =   a1[0] / diaga[k+0];
      rtemp6      =   a2[0] / diaga[k+1];

      b00[0]     -=   a1[0] * rtemp1
                    + a2[0] * rtemp2;

      diagb[i+1] -=   a1[0] * rtemp5
                    + a2[0] * rtemp6;


               ++ a1;
               ++ a2;

      for(t=0; t<sze; ++t) {
        b0[t] -=   a1[t] * rtemp1
                 + a2[t] * rtemp2;

        b1[t] -=   a1[t] * rtemp5
                 + a2[t] * rtemp6;
      }
    }

    for(; k<n; ++k) {
      a1          = a+(fira[k+0 ]+i);

      rtemp1      = a1[0] / diaga[k+0];

      diagb[i]   -= a1[0] * rtemp1;

                 ++ a1;

      rtemp5      = a1[0] / diaga[k+0];

      diagb[i+1] -= a1[0] * rtemp5;

      b00[0]     -= a1[0] * rtemp1;

                 ++ a1;

      for(t=0; t<sze; ++t)
      {
        b0[t] -= rtemp1 * a1[t];
        b1[t] -= rtemp5 * a1[t];
      }
    }
  }

  for(; i<s; ++i) {

    b0  = b+firb[i];
    sze = m-i-1;
    k   = 0;

    for(; k+15<n; k+=16) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);
      a3       = a+(fira[k+2 ]+i);
      a4       = a+(fira[k+3 ]+i);
      a5       = a+(fira[k+4 ]+i);
      a6       = a+(fira[k+5 ]+i);
      a7       = a+(fira[k+6 ]+i);
      a8       = a+(fira[k+7 ]+i);
      a9       = a+(fira[k+8 ]+i);
      a10      = a+(fira[k+9 ]+i);
      a11      = a+(fira[k+10]+i);
      a12      = a+(fira[k+11]+i);
      a13      = a+(fira[k+12]+i);
      a14      = a+(fira[k+13]+i);
      a15      = a+(fira[k+14]+i);
      a16      = a+(fira[k+15]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];
      rtemp3    = *a3/diaga[k+2];
      rtemp4    = *a4/diaga[k+3];
      rtemp5    = *a5/diaga[k+4];
      rtemp6    = *a6/diaga[k+5];
      rtemp7    = *a7/diaga[k+6];
      rtemp8    = *a8/diaga[k+7];
      rtemp9    = *a9/diaga[k+8];
      rtemp10   = *a10/diaga[k+9];
      rtemp11   = *a11/diaga[k+10];
      rtemp12   = *a12/diaga[k+11];
      rtemp13   = *a13/diaga[k+12];
      rtemp14   = *a14/diaga[k+13];
      rtemp15   = *a15/diaga[k+14];
      rtemp16   = *a16/diaga[k+15];

      diagb[i] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2)
                  + rtemp3  * (*a3)
                  + rtemp4  * (*a4)
                  + rtemp5  * (*a5)
                  + rtemp6  * (*a6)
                  + rtemp7  * (*a7)
                  + rtemp8  * (*a8)
                  + rtemp9  * (*a9)
                  + rtemp10 * (*a10)
                  + rtemp11 * (*a11)
                  + rtemp12 * (*a12)
                  + rtemp13 * (*a13)
                  + rtemp14 * (*a14)
                  + rtemp15 * (*a15)
                  + rtemp16 * (*a16);


               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;
               ++ a5;
               ++ a6;
               ++ a7;
               ++ a8;
               ++ a9;
               ++ a10;
               ++ a11;
               ++ a12;
               ++ a13;
               ++ a14;
               ++ a15;
               ++ a16;

      for(t=0; t<sze; ++t)
        b0[t] -=   rtemp1  * a1[t]
                 + rtemp2  * a2[t]
                 + rtemp3  * a3[t]
                 + rtemp4  * a4[t]
                 + rtemp5  * a5[t]
                 + rtemp6  * a6[t]
                 + rtemp7  * a7[t]
                 + rtemp8  * a8[t]
                 + rtemp9  * a9[t]
                 + rtemp10 * a10[t]
                 + rtemp11 * a11[t]
                 + rtemp12 * a12[t]
                 + rtemp13 * a13[t]
                 + rtemp14 * a14[t]
                 + rtemp15 * a15[t]
                 + rtemp16 * a16[t];

    }

    for(; k+11<n; k+=12) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);
      a3       = a+(fira[k+2 ]+i);
      a4       = a+(fira[k+3 ]+i);
      a5       = a+(fira[k+4 ]+i);
      a6       = a+(fira[k+5 ]+i);
      a7       = a+(fira[k+6 ]+i);
      a8       = a+(fira[k+7 ]+i);
      a9       = a+(fira[k+8 ]+i);
      a10      = a+(fira[k+9 ]+i);
      a11      = a+(fira[k+10]+i);
      a12      = a+(fira[k+11]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];
      rtemp3    = *a3/diaga[k+2];
      rtemp4    = *a4/diaga[k+3];
      rtemp5    = *a5/diaga[k+4];
      rtemp6    = *a6/diaga[k+5];
      rtemp7    = *a7/diaga[k+6];
      rtemp8    = *a8/diaga[k+7];
      rtemp9    = *a9/diaga[k+8];
      rtemp10   = *a10/diaga[k+9];
      rtemp11   = *a11/diaga[k+10];
      rtemp12   = *a12/diaga[k+11];

      diagb[i] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2)
                  + rtemp3  * (*a3)
                  + rtemp4  * (*a4)
                  + rtemp5  * (*a5)
                  + rtemp6  * (*a6)
                  + rtemp7  * (*a7)
                  + rtemp8  * (*a8)
                  + rtemp9  * (*a9)
                  + rtemp10 * (*a10)
                  + rtemp11 * (*a11)
                  + rtemp12 * (*a12);

               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;
               ++ a5;
               ++ a6;
               ++ a7;
               ++ a8;
               ++ a9;
               ++ a10;
               ++ a11;
               ++ a12;

      for(t=0; t<sze; ++t)
        b0[t] -=   rtemp1  * a1[t]
                 + rtemp2  * a2[t]
                 + rtemp3  * a3[t]
                 + rtemp4  * a4[t]
                 + rtemp5  * a5[t]
                 + rtemp6  * a6[t]
                 + rtemp7  * a7[t]
                 + rtemp8  * a8[t]
                 + rtemp9  * a9[t]
                 + rtemp10 * a10[t]
                 + rtemp11 * a11[t]
                 + rtemp12 * a12[t];

    }
    
    for(; k+7<n; k+=8) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);
      a3       = a+(fira[k+2 ]+i);
      a4       = a+(fira[k+3 ]+i);
      a5       = a+(fira[k+4 ]+i);
      a6       = a+(fira[k+5 ]+i);
      a7       = a+(fira[k+6 ]+i);
      a8       = a+(fira[k+7 ]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];
      rtemp3    = *a3/diaga[k+2];
      rtemp4    = *a4/diaga[k+3];
      rtemp5    = *a5/diaga[k+4];
      rtemp6    = *a6/diaga[k+5];
      rtemp7    = *a7/diaga[k+6];
      rtemp8    = *a8/diaga[k+7];

      diagb[i] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2)
                  + rtemp3  * (*a3)
                  + rtemp4  * (*a4)
                  + rtemp5  * (*a5)
                  + rtemp6  * (*a6)
                  + rtemp7  * (*a7)
                  + rtemp8  * (*a8);

               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;
               ++ a5;
               ++ a6;
               ++ a7;
               ++ a8;

      for(t=0; t<sze; ++t)
        b0[t] -=   rtemp1  * a1[t]
                 + rtemp2  * a2[t]
                 + rtemp3  * a3[t]
                 + rtemp4  * a4[t]
                 + rtemp5  * a5[t]
                 + rtemp6  * a6[t]
                 + rtemp7  * a7[t]
                 + rtemp8  * a8[t];

    }

    for(; k+3<n; k+=4) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);
      a3       = a+(fira[k+2 ]+i);
      a4       = a+(fira[k+3 ]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];
      rtemp3    = *a3/diaga[k+2];
      rtemp4    = *a4/diaga[k+3];

      diagb[i] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2)
                  + rtemp3  * (*a3)
                  + rtemp4  * (*a4);

               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;

      for(t=0; t<sze; ++t)
        b0[t] -=   rtemp1  * a1[t]
                 + rtemp2  * a2[t]
                 + rtemp3  * a3[t]
                 + rtemp4  * a4[t];

    }

    for(; k+1<n; k+=2) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];

      diagb[i] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2);

               ++ a1;
               ++ a2;

      for(t=0; t<sze; ++t)
        b0[t] -=   rtemp1  * a1[t]
                 + rtemp2  * a2[t];
    }

    for(; k<n; ++k) {
      a1        = a+(fira[k+0 ]+i);

      rtemp1    = *a1/diaga[k+0];

      diagb[i] -= rtemp1  * (*a1);

               ++ a1;

      for(t=0; t<sze; ++t)
        b0[t] -=   rtemp1  * a1[t];
    }
  }
} /*  UpdSnode */

static void iUpdSnode(chfac  *cf,
                        int    snde,
                        int    f,
                        int    l,
                        int    uf,
                        int    ul,
                        int    iw[])
{
  int    k,
         *ujsze=cf->ujsze,*uhead=cf->uhead,*subg=cf->subg;
  double *diag=cf->diag,*uval=cf->uval;

  if (f==l || uf==ul)
    return;

  f  += subg[snde];
  l  += subg[snde];
  uf += subg[snde];
  ul += subg[snde];

  for(k=f; k<l; ++k)
    iw[k-f] = uhead[k]+uf-k-1;

  UpdSnode(1+ujsze[uf],l-f,ul-uf,diag+f,uval,iw,diag+uf,uval,uhead+uf);
} /* iUpdSnode */

static int DiagUpdate(double *dii,
                      int    mode)
{
  int    id=TRUE;
  
  if (mode) {
    if (*dii<1.0e-13)
      return FALSE;
  }
  else {
    if (fabs(*dii)<1.0e-35) {
      printf(" diagonal nearly zero: %5.1e.\n",(*dii));
      return FALSE;
    }
  }
  return id;
} /* DiagUpdate */

static int FacSnode(chfac  *cf,
                    int    snde,
                    int    f,
                    int    l,
                    int    iw[],
                    int    mode)
{
  int    itemp,k;

  if (f==l)
    return (CfcOk);

  itemp = cf->subg[snde]+f;

  if (!DiagUpdate(cf->diag+itemp,
                  mode))
    return (CfcIndef);
  
  if (fabs(cf->diag[itemp])<=cf->tolpiv) {
      printf("Singular d[%d]=%e\n",
      cf->subg[snde]+f,cf->diag[cf->subg[snde]+f]);
      return (CfcIndef);
  }
   
  for(k=f+1; k<l; ++k) {
    iUpdSnode(cf,snde,f,k,k,k+1,iw);
     
    itemp = cf->subg[snde]+k;

    if (!DiagUpdate(&cf->diag[itemp],
                    mode))
      return (CfcIndef);
    
    if (fabs(cf->diag[itemp])<=cf->tolpiv) {
      printf(" singular d[%d]=%e\n",
             cf->subg[snde]+k,cf->diag[cf->subg[snde]+k]);

      return (CfcIndef);
   }
  }
  
  return CfcOk;
} /* FacSnode */

static void UpdSnodes(int    m,
                      int    n,
                      int    s,
                      double diaga[],
                      double *a,
                      int    fira[],
                      double diagb[],
                      double *b,
                      int    firb[],
                      int    subb[])
{
  int    i,j,k,t,u,sze,delay,
          *ls;
  double rtemp1,rtemp2,rtemp3,rtemp4,
         rtemp5,rtemp6,rtemp7,rtemp8,
         rtemp9,rtemp10,rtemp11,rtemp12,
         rtemp13,rtemp14,rtemp15,rtemp16,
         *a1,*a2,*a3,*a4,*a5,*a6,*a7,*a8,
         *a9,*a10,*a11,*a12,*a13,*a14,*a15,*a16,
         *b0;

  if (m<s)
    exit(0);

  if (m==0 || n==0)
    return;

  for(i=0; i<s; ++i) {
    j     = subb[i];
    u     = j-subb[0];

    b0    = b+firb[u];

    delay = j+1;
    sze   = m-i-1;
    ls    = subb+i+1;

    k     = 0;

    for(; k+15<n; k+=16) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);
      a3       = a+(fira[k+2 ]+i);
      a4       = a+(fira[k+3 ]+i);
      a5       = a+(fira[k+4 ]+i);
      a6       = a+(fira[k+5 ]+i);
      a7       = a+(fira[k+6 ]+i);
      a8       = a+(fira[k+7 ]+i);
      a9       = a+(fira[k+8 ]+i);
      a10      = a+(fira[k+9 ]+i);
      a11      = a+(fira[k+10]+i);
      a12      = a+(fira[k+11]+i);
      a13      = a+(fira[k+12]+i);
      a14      = a+(fira[k+13]+i);
      a15      = a+(fira[k+14]+i);
      a16      = a+(fira[k+15]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];
      rtemp3    = *a3/diaga[k+2];
      rtemp4    = *a4/diaga[k+3];
      rtemp5    = *a5/diaga[k+4];
      rtemp6    = *a6/diaga[k+5];
      rtemp7    = *a7/diaga[k+6];
      rtemp8    = *a8/diaga[k+7];
      rtemp9    = *a9/diaga[k+8];
      rtemp10   = *a10/diaga[k+9];
      rtemp11   = *a11/diaga[k+10];
      rtemp12   = *a12/diaga[k+11];
      rtemp13   = *a13/diaga[k+12];
      rtemp14   = *a14/diaga[k+13];
      rtemp15   = *a15/diaga[k+14];
      rtemp16   = *a16/diaga[k+15];

      diagb[u] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2)
                  + rtemp3  * (*a3)
                  + rtemp4  * (*a4)
                  + rtemp5  * (*a5)
                  + rtemp6  * (*a6)
                  + rtemp7  * (*a7)
                  + rtemp8  * (*a8)
                  + rtemp9  * (*a9)
                  + rtemp10 * (*a10)
                  + rtemp11 * (*a11)
                  + rtemp12 * (*a12)
                  + rtemp13 * (*a13)
                  + rtemp14 * (*a14)
                  + rtemp15 * (*a15)
                  + rtemp16 * (*a16);


               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;
               ++ a5;
               ++ a6;
               ++ a7;
               ++ a8;
               ++ a9;
               ++ a10;
               ++ a11;
               ++ a12;
               ++ a13;
               ++ a14;
               ++ a15;
               ++ a16;

      for(t=0; t<sze; ++t)
        b0[ls[t]-delay] -=   rtemp1  * a1[t]
                           + rtemp2  * a2[t]
                           + rtemp3  * a3[t]
                           + rtemp4  * a4[t]
                           + rtemp5  * a5[t]
                           + rtemp6  * a6[t]
                           + rtemp7  * a7[t]
                           + rtemp8  * a8[t]
                           + rtemp9  * a9[t]
                           + rtemp10 * a10[t]
                           + rtemp11 * a11[t]
                           + rtemp12 * a12[t]
                           + rtemp13 * a13[t]
                           + rtemp14 * a14[t]
                           + rtemp15 * a15[t]
                           + rtemp16 * a16[t];
    }
    
    for(; k+11<n; k+=12) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);
      a3       = a+(fira[k+2 ]+i);
      a4       = a+(fira[k+3 ]+i);
      a5       = a+(fira[k+4 ]+i);
      a6       = a+(fira[k+5 ]+i);
      a7       = a+(fira[k+6 ]+i);
      a8       = a+(fira[k+7 ]+i);
      a9       = a+(fira[k+8 ]+i);
      a10      = a+(fira[k+9 ]+i);
      a11      = a+(fira[k+10]+i);
      a12      = a+(fira[k+11]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];
      rtemp3    = *a3/diaga[k+2];
      rtemp4    = *a4/diaga[k+3];
      rtemp5    = *a5/diaga[k+4];
      rtemp6    = *a6/diaga[k+5];
      rtemp7    = *a7/diaga[k+6];
      rtemp8    = *a8/diaga[k+7];
      rtemp9    = *a9/diaga[k+8];
      rtemp10   = *a10/diaga[k+9];
      rtemp11   = *a11/diaga[k+10];
      rtemp12   = *a12/diaga[k+11];

      diagb[u] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2)
                  + rtemp3  * (*a3)
                  + rtemp4  * (*a4)
                  + rtemp5  * (*a5)
                  + rtemp6  * (*a6)
                  + rtemp7  * (*a7)
                  + rtemp8  * (*a8)
                  + rtemp9  * (*a9)
                  + rtemp10 * (*a10)
                  + rtemp11 * (*a11)
                  + rtemp12 * (*a12);

               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;
               ++ a5;
               ++ a6;
               ++ a7;
               ++ a8;
               ++ a9;
               ++ a10;
               ++ a11;
               ++ a12;

      for(t=0; t<sze; ++t)
        b0[ls[t]-delay] -=   rtemp1  * a1[t]
                           + rtemp2  * a2[t]
                           + rtemp3  * a3[t]
                           + rtemp4  * a4[t]
                           + rtemp5  * a5[t]
                           + rtemp6  * a6[t]
                           + rtemp7  * a7[t]
                           + rtemp8  * a8[t]
                           + rtemp9  * a9[t]
                           + rtemp10 * a10[t]
                           + rtemp11 * a11[t]
                           + rtemp12 * a12[t];
    }
    
    for(; k+7<n; k+=8) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);
      a3       = a+(fira[k+2 ]+i);
      a4       = a+(fira[k+3 ]+i);
      a5       = a+(fira[k+4 ]+i);
      a6       = a+(fira[k+5 ]+i);
      a7       = a+(fira[k+6 ]+i);
      a8       = a+(fira[k+7 ]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];
      rtemp3    = *a3/diaga[k+2];
      rtemp4    = *a4/diaga[k+3];
      rtemp5    = *a5/diaga[k+4];
      rtemp6    = *a6/diaga[k+5];
      rtemp7    = *a7/diaga[k+6];
      rtemp8    = *a8/diaga[k+7];

      diagb[u] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2)
                  + rtemp3  * (*a3)
                  + rtemp4  * (*a4)
                  + rtemp5  * (*a5)
                  + rtemp6  * (*a6)
                  + rtemp7  * (*a7)
                  + rtemp8  * (*a8);

               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;
               ++ a5;
               ++ a6;
               ++ a7;
               ++ a8;

      for(t=0; t<sze; ++t)
        b0[ls[t]-delay] -=   rtemp1  * a1[t]
                           + rtemp2  * a2[t]
                           + rtemp3  * a3[t]
                           + rtemp4  * a4[t]
                           + rtemp5  * a5[t]
                           + rtemp6  * a6[t]
                           + rtemp7  * a7[t]
                           + rtemp8  * a8[t];

    }

    for(; k+3<n; k+=4) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);
      a3       = a+(fira[k+2 ]+i);
      a4       = a+(fira[k+3 ]+i);

      rtemp1   = *a1/diaga[k+0];
      rtemp2   = *a2/diaga[k+1];
      rtemp3   = *a3/diaga[k+2];
      rtemp4   = *a4/diaga[k+3];

      diagb[u]-=   rtemp1  * (*a1)
                  +rtemp2  * (*a2)
                  +rtemp3  * (*a3)
                  +rtemp4  * (*a4);

               ++ a1;
               ++ a2;
               ++ a3;
               ++ a4;

      for(t=0; t<sze; ++t)
        b0[ls[t]-delay] -=   rtemp1  * a1[t]
                           + rtemp2  * a2[t]
                           + rtemp3  * a3[t]
                           + rtemp4  * a4[t];

    }

    for(; k+1<n; k+=2) {
      a1       = a+(fira[k+0 ]+i);
      a2       = a+(fira[k+1 ]+i);

      rtemp1    = *a1/diaga[k+0];
      rtemp2    = *a2/diaga[k+1];

      diagb[u] -=   rtemp1  * (*a1)
                  + rtemp2  * (*a2);

               ++ a1;
               ++ a2;

      for(t=0; t<sze; ++t)
        b0[ls[t]-delay] -=   rtemp1  * a1[t]
                           + rtemp2  * a2[t];
    }

    for(; k<n; ++k) {
      a1        = a+(fira[k+0 ]+i);

      rtemp1    = *a1/diaga[k+0];

      diagb[u] -= rtemp1  * (*a1);

               ++ a1;

      for(t=0; t<sze; ++t)
        b0[ls[t]-delay] -=   rtemp1  * a1[t];
    }
  }
} /*  UpdSnodes */

static void ExtUpdSnode(chfac  *cf,
                        int    snde,
                        int    usnde,
                        int    f,
                        int    l,
                        int    start,
                        int    iw[])
{
  int    k,sze,
          *ls,
          *subg=cf->subg,
          *ujsze=cf->ujsze,*usub=cf->usub,*ujbeg=cf->ujbeg,*uhead=cf->uhead;
  double *diag=cf->diag,*uval=cf->uval;

  f += subg[snde];
  l += subg[snde];

  if (usnde==cf->nsnds-1) {
    if (usub[ujbeg[f]+start]<subg[usnde]) {
      printf("\n Index error");
      exit(0);
    }

    if (cf->sdens)
      exit(0);

    ls  = usub+ujbeg[f]+start;
    sze = ujsze[f]-start;

    for(k=f; k<l; ++k)
      iw[k-f] =uhead[k]+start-(k-f);

    UpdSnodes(sze,l-f,sze,diag+f,uval,iw,diag+ls[0],uval,uhead+ls[0],ls);
  }
  else
    exit(0);
} /* ExtUpdSnode */

static void PushFward(chfac  *cf,
                      int    snde,
                      int    f,
                      int    l,
                      int    iw[])
{
  int    j,s,t,u,k,stops,offset,sze,itemp,
          *ls0,*ls1,
          *map=iw,*subg=cf->subg,
          *ujsze=cf->ujsze,*uhead=cf->uhead,*ujbeg=cf->ujbeg,*usub=cf->usub;
  double rtemp1,*l0,*l1,
          *diag=cf->diag,*uval=cf->uval;

  if (f>subg[snde+1]-subg[snde]) {
    printf("\n PushFward");
    exit(0);
  }

  if (f==l)
    return;

  f      += subg[snde];
  l      += subg[snde];

  offset  = subg[snde+1]-f-1;
  sze     = ujsze[f] - offset;
  ls1     = usub+ujbeg[f]+offset;

  if (f+1==l) {
    l1     = uval+uhead[f]+offset;

    stops = sze;
    for(t=0; t<sze; ++t) {
      j = ls1[0];

      if (j>=subg[cf->nsnds-1])
        break;

      rtemp1   = l1[0]/diag[f];
      diag[j] -= rtemp1*l1[0];
              ++ l1;

      l0       = uval+uhead[j];
      ls0      = usub+ujbeg[j];

              ++ ls1;
              -- stops;

      if (stops && ls1[stops-1]==ls0[stops-1]) {
        for(s=0; s<stops; ++s)
          l0[s] -= rtemp1 * l1[s];
      }
      
      else {
        for(s=0, u=0; s<stops; ++u) {
          if (ls0[u]==ls1[s]) {
            l0[u] -= rtemp1 * l1[s];
                  ++ s;
          }
        }
      }
    }

    if (t<sze && !cf->sdens)
      ExtUpdSnode(cf,snde,cf->nsnds-1,f-subg[snde],l-subg[snde],t,iw);
  }
  
  else {
    stops = sze;
    for(t=0; t<sze; ++t, ++offset) {
      j      = ls1[0];

      if (j>=subg[cf->nsnds-1]) {
        if (!cf->sdens)
          ExtUpdSnode(cf,snde,cf->nsnds-1,f-subg[snde],
                      l-subg[snde],offset,iw);
        break;
      }

      ls0   = usub+ujbeg[j];
      l0    = uval+uhead[j];

           ++ ls1;
           -- stops;

      k     = f;
      itemp = offset+f;

      if (stops && ls1[stops-1]==ls0[stops-1]) {
        for(k=f; k<l; ++k)
          map[k-f] = uhead[k]+itemp-k;

        UpdSnode(1+stops,l-f,1,diag+f,uval,map,diag+j,uval,uhead+j);
      }
      
      else {
        map[l] = 0;
        for(s=0, u=0; s<stops; ++u) {
          if (ls1[s]==ls0[u]) {
            map[1+l+s]  = 1+u;
                       ++ s;
          }
        }

        for(k=f; k<l; ++k)
          map[k-f] = uhead[k]+itemp-k;

        UpdSnodes(1+stops,l-f,1,diag+f,uval,map,diag+j,uval,uhead+j,map+l);
      }
    }
  }
} /* PushFwardard */

static int FacDenNode(chfac  *cf,
                      int    iw[],
                      double rw[],
                      int    mode)
{
  int    c,d,j,s,t,sncl,k,k0,m,cacsze,sze,offset,
         *subg=cf->subg,*ujsze=cf->ujsze,
         *ujbeg=cf->ujbeg,*uhead=cf->uhead,
         *usub=cf->usub,*dhead=cf->dhead,
         *dsub=cf->dsub,*dbeg=cf->dbeg,*ls;
  double *diag=cf->diag,*uval=cf->uval;
  int    sresp;

  cacsze = cf->cachesize*cf->cacheunit;

  if (cf->sdens) {
    for(d=0; d<cf->ndens; ++d) {
      c = 0;
      for(k=dhead[d]; k<dhead[d+1]; ++k) {
        offset = dbeg[k];
        s      = dsub[k];
        if (usub[ujbeg[subg[s]]+offset]<subg[cf->nsnds-1]) {
          printf("\nindex error1");
          exit(0);
        }

        for(j=subg[s]; j<subg[s+1]; ++c, ++j) {
          rw[c] = diag[j];
          iw[c] = uhead[j]+offset-(j-subg[s]);

          if (usub[ujbeg[j]+offset-(j-subg[s])]<subg[cf->nsnds-1]) {
            printf("\nindex error");
            exit(0);
          }
        }
      }

      if (c) {
        k  = dhead[d];
        s  = dsub[k];
        m  = ujsze[subg[s]]-dbeg[k];
        ls = usub+ujbeg[subg[s]]+dbeg[k];
        if (m) {
          for(k=0; k<c;) {
            t = cacsze/(m*sizeof(double));
            t = max(t,1);
            t = min(c-k,t);

            UpdSnodes(m,t,m,rw+k,uval,iw+k,
                      diag+ls[0],uval,uhead+ls[0],ls);
            k+=t;
          }
        }
      }
    }
  }

  s = cf->nsnds-1;

  sncl = cf->subg[s+1]-cf->subg[s];
  for(k=0; k<sncl;) {
    k0  = k;
    for(sze=0; sze<cacsze && k<sncl; ++k)
      sze += ujsze[subg[s]+k] * sizeof(double);

    if (k==k0)
      ++k;
    else if (k>=k0+2 && sze>cacsze)
      --k;

    if (k>sncl)
      exit(0);

    sresp = FacSnode(cf,s,k0,k,iw,mode);
    if (sresp!=CfcOk)
      return (sresp);

    iUpdSnode(cf,s,k0,k,k,sncl,iw);

  }
  return (CfcOk);
} /* FacDenNode */

int ChlFact(chfac  *cf,
            int    *iw,
            double *rw,
            int    mode)
{
  int s,sncl,k,k0,cacsze,sze,
      *subg=cf->subg,*ujsze=cf->ujsze;
  int cid;

  cacsze=cf->cachesize*cf->cacheunit;

  for(s=0; s+1<cf->nsnds; ++s) {
    sncl = cf->subg[s+1]-cf->subg[s];

    for(k=0; k<sncl;) {
      k0  = k;
      for(sze=0; sze<=cacsze && k<sncl; ++k)
        sze += ujsze[subg[s]+k]*sizeof(double);
       
      if (k==k0)
        ++k;
      else if (k>=k0+2 && sze>cacsze)
        --k;
      
      if (k>sncl)
        exit(0);

      cid=FacSnode(cf,s,k0,k,iw,mode);
      if (cid!=CfcOk)
        return (cid);

      iUpdSnode(cf,s,k0,k,k,sncl,iw);

      PushFward(cf,s,k0,k,iw);
    }
  }

  cid=FacDenNode(cf,iw,rw,mode);

  for (k=0;k<cf->nrow;k++){
    cf->sqrtdiag[k]=sqrt(fabs(cf->diag[k]));
  }

  return cid;
} /* ChlFact */
