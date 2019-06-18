#include "dsdp5.h"
/*!
\file printsdpa.c 
\brief Print data or solution in SDPA format
*/
static void DprintfD(FILE*fp, double d1){
  int i1=(int)d1,i2=(int)(d1*100),i3=(int)(d1*10000),i4=(int)(d1*10000000);
  if ( d1==(double)i1){ fprintf(fp,"%2.0f ",d1);
  } else if (d1==(double)i2/100.0){ fprintf(fp,"%4.2f ",d1);
  } else if (d1==(double)i3/10000.0){ fprintf(fp,"%6.4f ",d1);
  } else if (d1==(double)i4/1000000.0){ fprintf(fp,"%8.6f ",d1);
  } else { fprintf(fp,"%22.22e ",d1);
  }
  return;
}

static void Dprintf(FILE*fp, int i1, int i2, int i3, int i4, double d1){
  if (fabs(d1)<1e-30){return;}
  fprintf(fp,"%d %d %d %d ",i1,i2,i3+1,i4+1);
  if (i1==0)
    DprintfD(fp,-d1);
  else 
    DprintfD(fp,d1);
  fprintf(fp,"\n");
  return;
}


#undef __FUNCT__  
#define __FUNCT__ "DPrintUpperRightMatrix"
static void DPrintUpperRightMatrix(int constraint, int block, double amat[], int n, FILE*fp){
  int col,row;
  for (row=0;row<n;row++){
    for (col=0;col<=row;col++){
      if (fabs(amat[col])>1.0e-20){
	Dprintf(fp,constraint,block,col,row,amat[col]);
      }
    }
    amat+=(row+1);
  }
  return;
}

#undef __FUNCT__  
#define __FUNCT__ "DPrintUpperFullMatrix"
static void DPrintUpperFullMatrix(int constraint, int block, double amat[], int n, FILE*fp){
  int col,row;
  for (row=0;row<n;row++){
    for (col=0;col<=row;col++){
      if (fabs(amat[col])>1.0e-20){
	Dprintf(fp,constraint,block,col,row,amat[col]);
      }
    }
    amat+=n;
  }
  return;
}

#undef __FUNCT__  
#define __FUNCT__ "DPrintMatrix"
static void DPrintMatrix(char UPLQ, int constraint, int block, double amat[], int n, FILE*fp){
  if (UPLQ=='P'){
    DPrintUpperRightMatrix(constraint,block,amat,n,fp);
  } else if (UPLQ=='U'){
    DPrintUpperFullMatrix(constraint,block,amat,n,fp);
  }
  return;
}

#undef __FUNCT__  
#define __FUNCT__ "DPrintLPArray"
static int DPrintLPArray(int cc, int block, double *vv, int n, FILE *fp ){
  int i;
  for (i=0;i<n;i++){
    if ( fabs(vv[i]) > 0){
      Dprintf(fp,cc,block,i,i,vv[i]);
    }
  }
  return 0;
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPPrintSolution"
/*!
\fn int DSDPPrintSolution(FILE *fp,DSDP dsdp,SDPCone sdpcone, LPCone lpcone);
\brief Print solution in SDPA format
\param fp file pointer
\param dsdp the solver
\param sdpcone semidefinite cone
\param lpcone LP cone
\ingroup Examples
 */
int DSDPPrintSolution(FILE *fp,DSDP dsdp,SDPCone sdpcone, LPCone lpcone){
  int     i,kk,info,n,nn,lpn=0,nblocks,nvars;
  double  *ss,*xx,*y,*lparray;
  char UPLQ;

  info=DSDPGetNumberOfVariables(dsdp,&nvars);DSDPCHKERR(info);
  DSDPCALLOC2(&y,double,(nvars+2),&info);DSDPCHKERR(info);
  info=SDPConeGetNumberOfBlocks(sdpcone,&nblocks);DSDPCHKERR(info);

  if (lpcone){info=LPConeGetXArray(lpcone,&xx,&lpn);DSDPCHKERR(info);nblocks--;}
  DSDPCALLOC2(&lparray,double,(lpn+1),&info);DSDPCHKERR(info);
  /* Deleted at Borcher's request.
  fprintf(fp,"%d \n%d \n",nvars,nblocks);
  for (i=0;i<nblocks; i++){
    info=SDPConeGetBlockSize(sdpcone,i,&n);DSDPCHKERR(info);
    fprintf(fp,"%d ",n);
  }
  if (lpcone){ fprintf(fp,"%d ",-lpn); }
  fprintf(fp," \n");
  */

  info=DSDPGetY(dsdp,y+1,nvars);DSDPCHKERR(info);
  y[0]=1.0;y[nvars+1]=0;
  info=DSDPGetR(dsdp,y+nvars+1);DSDPCHKERR(info);
  for (i=0; i<nvars; i++){ DprintfD(fp,-y[i+1]);}
  fprintf(fp," \n");

  /* Print Dual Matrix Solution */
  for (kk=0;kk<nblocks;kk++){
    info=SDPConeGetBlockSize(sdpcone,kk,&n);DSDPCHKERR(info);
    info=SDPConeGetXArray(sdpcone,kk,&ss,&nn);DSDPCHKERR(info);
    info=SDPConeComputeS(sdpcone,kk,y[0],y+1,nvars,y[nvars+1],n,ss,nn);DSDPCHKERR(info);
    info=SDPConeGetStorageFormat(sdpcone,kk,&UPLQ);DSDPCHKERR(info);
    DPrintMatrix(UPLQ,1,kk+1,ss,n,fp);
    info=SDPConeRestoreXArray(sdpcone,kk,&ss,&nn);DSDPCHKERR(info);
  }
  if (lpcone){
    info=LPConeCopyS(lpcone,lparray,lpn);DSDPCHKERR(info);
    info=DPrintLPArray(1,nblocks+1,lparray,lpn,fp);DSDPCHKERR(info);
  }

  info=DSDPComputeX(dsdp);DSDPCHKERR(info);
  /* Print Primal Solution */ 
  for (kk=0; kk<nblocks; kk++){
    info=SDPConeGetBlockSize(sdpcone,kk,&n);DSDPCHKERR(info);
    info=SDPConeGetStorageFormat(sdpcone,kk,&UPLQ);DSDPCHKERR(info);
    info=SDPConeGetXArray(sdpcone,kk,&xx,&nn);DSDPCHKERR(info);
    DPrintMatrix(UPLQ,2,kk+1,xx,n,fp);
    info=SDPConeRestoreXArray(sdpcone,kk,&xx,&nn);DSDPCHKERR(info);
  }
  if (lpcone){
    info=LPConeGetXArray(lpcone,&xx,&nn);DSDPCHKERR(info);
    info=DPrintLPArray(2,nblocks+1,xx,nn,fp);DSDPCHKERR(info);
  }

  DSDPFREE(&y,&info);
  return 0;
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPPrintData"
/*!
\fn int DSDPPrintData(DSDP dsdp,SDPCone sdpcone, LPCone lpcone);
\brief Print data in SDPA format to a file named "output.sdpa"
\param dsdp the solver
\param sdpcone semidefinite cone
\param lpcone LP cone
\ingroup Examples
 */
int DSDPPrintData(DSDP dsdp, SDPCone sdpcone, LPCone lpcone){

  int info,nblocks,i,nvars,n,nn,kk,lpblock=0,lpn=0;
  double *ss,*y=0,*vv=0;
  char filename[100]="";
  char UPLQ;
  FILE *fp;

  info=DSDPGetNumberOfVariables(dsdp,&nvars);DSDPCHKERR(info);
  DSDPCALLOC2(&y,double,(nvars+3),&info);DSDPCHKERR(info);
  info=SDPConeGetNumberOfBlocks(sdpcone,&nblocks);DSDPCHKERR(info);
  strncat(filename,"output.sdpa",50);
  /*  fp=fopen(filename,"w"); */
  fp=fopen("input.sdpa","w");
  if (lpcone){ 
    info=LPConeGetDimension(lpcone,&lpn);DSDPCHKERR(info);
    DSDPCALLOC2(&vv,double,lpn,&info);DSDPCHKERR(info);
    lpblock=1;
    info=SDPConeGetBlockSize(sdpcone,nblocks-1,&n);DSDPCHKERR(info);
    if (n==0){nblocks--;}
  }
  fprintf(fp,"%d \n%d\n",nvars,nblocks+lpblock);
  for (kk=0;kk<nblocks;kk++){ 
    info=SDPConeGetBlockSize(sdpcone,kk,&n);DSDPCHKERR(info);
    fprintf(fp,"%d ",n);
  }
  if (lpcone){
    fprintf(fp,"%d ",-lpn);
  }
  fprintf(fp,"\n");
  info=DSDPCopyB(dsdp,y,nvars);
  for (i=0;i<nvars;i++){ 
    DprintfD(fp,y[i]);
  }
  fprintf(fp,"\n");

  for (i=0;i<=nvars;i++){
    for (kk=0;kk<nvars+2;kk++) y[kk]=0.0;
    if (i==0){y[i]=1.0;} else {y[i]=-1.0; }
    for (kk=0;kk<nblocks;kk++){
      info=SDPConeGetBlockSize(sdpcone,kk,&n);DSDPCHKERR(info);
      info=SDPConeGetXArray(sdpcone,kk,&ss,&nn);DSDPCHKERR(info);
      info=SDPConeComputeS(sdpcone,kk,y[0],y+1,nvars,y[nvars+1],n,ss,nn);DSDPCHKERR(info);
      info=SDPConeGetStorageFormat(sdpcone,kk,&UPLQ);DSDPCHKERR(info);
      DPrintMatrix(UPLQ,i,kk+1,ss,n,fp);
    }
  }
  if (lpcone && lpn>0){
    info=LPConeGetDimension(lpcone,&lpn);DSDPCHKERR(info);
    for (i=0;i<=nvars;i++){
      info=LPConeGetData(lpcone,i,vv,lpn);DSDPCHKERR(info);
      info=DPrintLPArray(i,nblocks+1,vv,lpn,fp);DSDPCHKERR(info);
    }
  }
  DSDPFREE(&y,&info);
  if (vv){
    DSDPFREE(&vv,&info);
  }
  fclose(fp);

  return 0;
}
