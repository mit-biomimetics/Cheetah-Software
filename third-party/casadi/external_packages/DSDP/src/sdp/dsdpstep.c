#include "dsdpdualmat.h"
#include "dsdpdsmat.h"
#include "dsdpxmat.h"
#include "dsdpsys.h"
#include "dsdplanczos.h"
#include "dsdplapack.h"

/*! 
\file dsdpstep.c
\brief Lanczos procedure determines the maximum step length
 */

typedef struct _P_Mat3* Mat3;

static int MatMult3(Mat3 A, SDPConeVec x, SDPConeVec y);
static int ComputeStepROBUST(Mat3 A, SDPConeVec *Q, int m, SDPConeVec W, SDPConeVec R,double*, SDPConeVec QAQTv, double *dwork, double *maxstep, double *mineig);
static int ComputeStepFAST(Mat3 A, SDPConeVec *Q, int m, SDPConeVec W, double *dwork, int*iwork,double *maxstep, double *mineig);

extern int DSDPGetEigsSTEP(double[],int,double[],int,long int[],int, 
		       double[],int,double[],int,int[],int);

int DSDPGetTriDiagonalEigs(int n,double *D,double *E, double*WORK2N,int*);

struct _P_Mat3{
  int type;
  DSDPDualMat ss;
  DSDPDSMat ds;
  SDPConeVec V;
  DSDPVMat x;
};


int DSDPGetTriDiagonalEigs(int N,double D[],double E[], double WORK[],int IIWORK[]){
  ffinteger LDZ=DSDPMax(1,N),INFO,NN=N;
  ffinteger M,IL=N-1,IU=N,*ISUPPZ=0;
  ffinteger *IWORK=(ffinteger*)IIWORK,LIWORK,LWORK;
  double WW[2],VL=-1e10,VU=1e10,*Z=0,ABSTOL=0;
  char JOBZ='N', RANGE='I';
  if (N<50){
    dstev(&JOBZ,&NN,D,E,Z,&LDZ,WORK,&INFO);
  } else {
    
    if (N<=1) IL=1;
    if (N<=1) IU=1;
    
    LWORK=20*N+1;
    LIWORK=10*N+1;

    dstevr(&JOBZ,&RANGE,&NN,D,E,&VL,&VU,&IL,&IU,&ABSTOL,&M,WW,Z,&LDZ,ISUPPZ,WORK,&LWORK,IWORK,&LIWORK,&INFO);
    
    if (N==1){
      D[0]=WW[0];
    } else if (N>=2){
      D[N-2]=WW[0];
      D[N-1]=WW[1];
    }

  }
  return INFO;
}

/* static int id1=0, id2=0; */
#undef __FUNCT__  
#define __FUNCT__ "MatMult3"
static int MatMult3(Mat3 A, SDPConeVec X, SDPConeVec Y){

  int info=0;
  double minus_one=-1.0;

  DSDPFunctionBegin;
  /*  DSDPEventLogBegin(id2); */
  if (A->type==2){
    info=DSDPVMatMult(A->x,X,Y);DSDPCHKERR(info);
  } else {
    info=DSDPDualMatCholeskySolveBackward(A->ss,X,Y); DSDPCHKERR(info);
    info=DSDPDSMatMult(A->ds,Y,A->V); DSDPCHKERR(info);
    info=DSDPDualMatCholeskySolveForward(A->ss,A->V,Y); DSDPCHKERR(info);
    info=SDPConeVecScale(minus_one,Y); DSDPCHKERR(info);
  }
  /*  DSDPEventLogEnd(id2);*/
  DSDPFunctionReturn(0); 
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPLanczosInitialize"
/*!
\fn int DSDPLanczosInitialize( DSDPLanczosStepLength *LZ );
\brief Initialize Lanczos structure.
\param LZ structure.
*/
int DSDPLanczosInitialize( DSDPLanczosStepLength *LZ ){
  DSDPFunctionBegin;
  /* Build Lanczos structures */
  LZ->n=0;
  LZ->lanczosm=0;
  LZ->maxlanczosm=20;
  LZ->type=0;
  LZ->dwork4n=0;
  LZ->Q=0;
  LZ->darray=0;
  /*
  if (id1==0 && id2==0){
  DSDPEventLogRegister("STEP EIGS",&id1); printf("ID1: %d\n",id1);
  DSDPEventLogRegister("STEP MULT",&id2); printf("ID2: %d\n",id2);
  }
  */
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPSetMaximumLanczosIterations"
/*!
\fn int DSDPSetMaximumLanczosIterations( DSDPLanczosStepLength *LZ, int maxlanczos );
\brief Set parameter.
\param LZ structure.
\param maxlanczos a parameter.
*/
int DSDPSetMaximumLanczosIterations( DSDPLanczosStepLength *LZ, int maxlanczos ){
  DSDPFunctionBegin;
  if (maxlanczos>0) LZ->lanczosm=maxlanczos;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPFastLanczosSetup"
/*!
\fn int DSDPFastLanczosSetup( DSDPLanczosStepLength *LZ, SDPConeVec V );
\brief Use Lanczos procedure. Assume off tridiagonal entries are zero.
\param LZ structure.
\param V work vector.
*/
int DSDPFastLanczosSetup( DSDPLanczosStepLength *LZ, SDPConeVec V ){
  int i,n,info;
  DSDPFunctionBegin;
  /* Build Lanczos structures */
  info=SDPConeVecGetSize(V,&n);DSDPCHKERR(info);
  LZ->lanczosm=DSDPMin(LZ->maxlanczosm,n);
  LZ->type=1;
  LZ->n=n;
  if (LZ->lanczosm<50){
    DSDPCALLOC2(&LZ->dwork4n,double,4*(LZ->lanczosm)+2,&info); DSDPCHKERR(info);
    DSDPCALLOC2(&LZ->iwork10n,int,1,&info); DSDPCHKERR(info);
  } else {
    DSDPCALLOC2(&LZ->dwork4n,double,23*(LZ->lanczosm)+2,&info); DSDPCHKERR(info);
    DSDPCALLOC2(&LZ->iwork10n,int,10*(LZ->lanczosm),&info); DSDPCHKERR(info);
  }
  DSDPCALLOC2(&LZ->Q,SDPConeVec,2,&info);DSDPCHKERR(info);
  for (i=0;i<2;i++){
    info = SDPConeVecDuplicate(V,&LZ->Q[i]);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPRobustLanczosSetup"
/*!
\fn int DSDPRobustLanczosSetup( DSDPLanczosStepLength *LZ, SDPConeVec V );
\brief Use slowerer but more robust method.
\param LZ structure.
\param V work vector.
*/
int DSDPRobustLanczosSetup( DSDPLanczosStepLength *LZ, SDPConeVec V ){
  int i,n,leig,info;
  DSDPFunctionBegin;
  /* Build Lanczos structures */
  info=SDPConeVecGetSize(V,&n);DSDPCHKERR(info);
  leig=DSDPMin(LZ->maxlanczosm,n);
  LZ->n=n;
  LZ->lanczosm=leig;
  LZ->type=2;

  DSDPCALLOC2(&LZ->dwork4n,double,3*(LZ->lanczosm)+1,&info); DSDPCHKERR(info);
  DSDPCALLOC2(&LZ->darray,double,(leig*leig),&info); DSDPCHKERR(info);
  DSDPCALLOC2(&LZ->Q,SDPConeVec,(leig+1),&info);DSDPCHKERR(info);

  for (i=0;i<=leig;i++){
    info = SDPConeVecDuplicate(V,&LZ->Q[i]);DSDPCHKERR(info);
  }
  info = SDPConeVecCreate(leig,&LZ->Tv);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPLanczosDestroy"
/*!
\fn int DSDPLanczosDestroy( DSDPLanczosStepLength *LZ);
\brief Free data structure.
\param LZ structure.
*/
int DSDPLanczosDestroy( DSDPLanczosStepLength *LZ){
  int i,info;
  DSDPFunctionBegin;
  if (LZ->type==2){
    for (i=0;i<=LZ->lanczosm;i++){
      info = SDPConeVecDestroy(&LZ->Q[i]);DSDPCHKERR(info);
    }
    info=SDPConeVecDestroy(&LZ->Tv);DSDPCHKERR(info);
    DSDPFREE(&LZ->darray,&info);DSDPCHKERR(info);
  } else if (LZ->type==1){
    info = SDPConeVecDestroy(&LZ->Q[1]);DSDPCHKERR(info);
    info = SDPConeVecDestroy(&LZ->Q[0]);DSDPCHKERR(info);
    DSDPFREE(&LZ->iwork10n,&info);DSDPCHKERR(info);
  }
  DSDPFREE(&LZ->Q,&info);DSDPCHKERR(info);
  DSDPFREE(&LZ->dwork4n,&info);DSDPCHKERR(info);
  info=DSDPLanczosInitialize(LZ);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPLanczosMinXEig"
int DSDPLanczosMinXEig( DSDPLanczosStepLength *LZ, DSDPVMat X, SDPConeVec W1, SDPConeVec W2, double *mineig ){
  int info,m;
  double smaxstep;
  struct _P_Mat3 PP;
  Mat3 A=&PP; 

  DSDPFunctionBegin;
  A->type=2;
  A->x=X;
  A->V=W2;
  m=LZ->lanczosm;

  if (LZ->type==1){
    info = ComputeStepFAST(A,LZ->Q,m,W1,LZ->dwork4n,LZ->iwork10n,&smaxstep,mineig);DSDPCHKERR(info);
  } else if (LZ->type==2){
    info = ComputeStepROBUST(A,LZ->Q,m,LZ->Q[m],W1,LZ->darray/*LZ->TT*/,LZ->Tv,LZ->dwork4n,&smaxstep,mineig);DSDPCHKERR(info);
  } else {
    DSDPSETERR1(1,"Lanczos Step Length Has not been SetUp. Type: %d\n",LZ->type);
  }
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPLanczosStepSize"
/*!
\fn int DSDPLanczosStepSize( DSDPLanczosStepLength *LZ, SDPConeVec W1, SDPConeVec W2, DSDPDualMat S, DSDPDSMat DS, double *maxstep );
\brief Compute distance to boundary
\param LZ structure.
\param W1 work vector
\param W2 work vector
\param S Current point in factored form.
\param DS Step direction.
\param maxstep output is distance to boundary.
*/
int DSDPLanczosStepSize( DSDPLanczosStepLength *LZ, SDPConeVec W1, SDPConeVec W2, DSDPDualMat S, DSDPDSMat DS, double *maxstep ){
  int info,m;
  double smaxstep,mineig;
  struct _P_Mat3 PP;
  Mat3 A=&PP; 

  DSDPFunctionBegin;
  A->ss=S;
  A->ds=DS; A->V=W2;
  A->type=1;
  m=LZ->lanczosm;

  if (LZ->type==1){
    info = ComputeStepFAST(A,LZ->Q,m,W1,LZ->dwork4n,LZ->iwork10n,&smaxstep,&mineig);DSDPCHKERR(info);
    *maxstep=smaxstep;
  } else if (LZ->type==2){
    info = ComputeStepROBUST(A,LZ->Q,m,LZ->Q[m],W1,LZ->darray/*LZ->TT*/,LZ->Tv,LZ->dwork4n,&smaxstep,&mineig);DSDPCHKERR(info);
    *maxstep=smaxstep;
  } else {
    DSDPSETERR1(1,"Lanczos Step Length Has not been SetUp. Type: %d\n",LZ->type);
  }
  DSDPFunctionReturn(0); 
}



#undef __FUNCT__  
#define __FUNCT__ "ComputeStepROBUST"
static int ComputeStepROBUST(Mat3 A, SDPConeVec *Q, int m, SDPConeVec W, SDPConeVec R, double*darray, SDPConeVec QAQTv, double *dwork, double *maxstep , double *mineig){

  int i,j,n,info;
  double tt,wnorm, phi;
  double lambda1=0,lambda2=0,delta=0;
  double res1,res2,beta;
  double one=1.0;
  double *eigvec;
  int N, LDA, LWORK;

  DSDPFunctionBegin;

  memset((void*)darray,0,m*m*sizeof(double));
  if (A->type==1){
    for (i=0; i<m; i++){ darray[i*m+i]=-1.0;}
  } else {
    for (i=0; i<m; i++){ darray[i*m+i]=1.0;}
  }
  
  info = SDPConeVecSet(one,Q[0]);DSDPCHKERR(info);
  info = SDPConeVecNormalize(Q[0]);DSDPCHKERR(info);

  for (i=0; i<m; i++){
    info = MatMult3(A,Q[i],W);DSDPCHKERR(info);
    info = SDPConeVecNorm2(W,&phi);DSDPCHKERR(info);
    if (phi!=phi){ *maxstep = 0.0;  return 0;} 
    if (i>0){
      tt=-darray[i*m+i-1];
      info = SDPConeVecAXPY(tt,Q[i-1],W);DSDPCHKERR(info);
    }
    info = SDPConeVecDot(W,Q[i],&tt);DSDPCHKERR(info);
    darray[i*m+i]=tt;
    tt*=-1.0;
    info = SDPConeVecAXPY(tt,Q[i],W);DSDPCHKERR(info);
    info = SDPConeVecNorm2(W,&wnorm);DSDPCHKERR(info);
    if (wnorm <= 0.8 * phi){
      for (j=0;j<=i;j++){
	info = SDPConeVecDot(W,Q[j],&tt);DSDPCHKERR(info);
	if (tt==tt){tt*=-1.0;} else {tt=0;}
	info = SDPConeVecAXPY(tt,Q[j],W);DSDPCHKERR(info);
	darray[j*m+i]-=tt;
	if (i!=j){ darray[i*m+j]-=tt; }
      }
    }

    info = SDPConeVecNorm2(W,&wnorm);DSDPCHKERR(info);
    if (i<m-1){
      darray[i*m+i+1]=wnorm;
      darray[i*m+m+i]=wnorm;
    }
    if (fabs(wnorm)<=1.0e-14) break;
    info=SDPConeVecCopy(W,Q[i+1]);DSDPCHKERR(info);
    info=SDPConeVecNormalize(Q[i+1]); DSDPCHKERR(info);

  }
  /*
  DSDPLogInfo("Step Length: Lanczos Iterates: %2.0f, ",i);
  DSDPLogInfo("VNorm: %3.1e, ",wnorm);
  */
  /*
  printf("  ---  TRI DIAGONAL MATRIX ---- \n");
  */


  LWORK=DSDPMax(3*m,1); LDA=DSDPMax(1,m); N=m;
  info=SDPConeVecGetArray(QAQTv,&eigvec);DSDPCHKERR(info);
  info=DSDPGetEigsSTEP(darray,m,0,0,0,0,eigvec,m,dwork,LWORK,0,0);DSDPCHKERR(info);
  info=SDPConeVecRestoreArray(QAQTv,&eigvec);DSDPCHKERR(info);

  if (N==0){
    lambda1=-0.0;
    delta=1.0e-20;
    *mineig=0;
  } else if (N==1){
    info=SDPConeVecGetArray(QAQTv,&eigvec);DSDPCHKERR(info);
    lambda1=-eigvec[0];
    info=SDPConeVecRestoreArray(QAQTv,&eigvec);DSDPCHKERR(info);
    delta=1.0e-20;
    *mineig=lambda1;
  } else if (N>1){
    info=SDPConeVecGetArray(QAQTv,&eigvec);DSDPCHKERR(info);
    *mineig=eigvec[0];
    lambda1=-eigvec[N-1];
    lambda2=-eigvec[N-2];
    info=SDPConeVecRestoreArray(QAQTv,&eigvec);DSDPCHKERR(info);

    info = SDPConeVecZero(W);DSDPCHKERR(info);
    for (i=0;i<m;i++){
      tt=darray[(N-1)*m+i];
      info=SDPConeVecAXPY(tt,Q[i],W);DSDPCHKERR(info);
    }
    info = MatMult3(A,W,R);DSDPCHKERR(info);
    info = SDPConeVecAXPY(lambda1,W,R);DSDPCHKERR(info);
    info = SDPConeVecNorm2(R,&res1);DSDPCHKERR(info);
        
    info = SDPConeVecZero(W);DSDPCHKERR(info);
    for (i=0;i<m;i++){
      tt=darray[(N-2)*m+i];
      info=SDPConeVecAXPY(tt,Q[i],W);DSDPCHKERR(info);
    }
    info = MatMult3(A,W,R);DSDPCHKERR(info);
    info = SDPConeVecAXPY(lambda2,W,R);DSDPCHKERR(info);
    info = SDPConeVecNorm2(R,&res2);DSDPCHKERR(info);
     
    tt = -lambda1 + lambda2 - res2;
    if (tt>0) beta=tt;
    else beta=1.0e-20;
    delta = DSDPMin(res1,sqrt(res1)/beta);
        
  }


  if (delta-lambda1>0)
    *maxstep = 1.0/(delta-lambda1);
  else
    *maxstep = 1.0e+30;

  info=SDPConeVecGetSize(W,&n);DSDPCHKERR(info);
  DSDPLogInfo(0,19,"Robust Lanczos StepLength: Iterates %d, Max: %d, BlockSize: %d, Lambda1: %4.2e, Res1: %4.2e, Lambda2: %4.2e, Res2: %4.2e, Delta: %4.2e, MaxStep: %4.2e\n",i,m,n,lambda1,res1*res1,lambda2,res2*res2,delta,*maxstep);

  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "ComputeStepFAST"
static int ComputeStepFAST(Mat3 A, SDPConeVec *Q, int m, SDPConeVec W, double *dwork, int *iwork,double *maxstep ,double *mineig){

  int i,j,n,info;
  double tt,wnorm, phi;
  double lambda1=0,lambda2=0,delta=0;
  double res1,res2,beta;
  double one=1.0;
  int N=m;
  double *diag,*subdiag,*ddwork;

  DSDPFunctionBegin;
  diag=dwork;
  subdiag=dwork+m;
  ddwork=dwork+2*m;

  if (A->type==1){
    for (i=0; i<m; i++){ diag[i]=-1; subdiag[i]=0;}
  } else {
    for (i=0; i<m; i++){ diag[i]=1.0; subdiag[i]=0;}
  }
  info = SDPConeVecSet(one,Q[0]);DSDPCHKERR(info);
  info = SDPConeVecNormalize(Q[0]);DSDPCHKERR(info);

  for (i=0; i<m; i++){
    info = MatMult3(A,Q[0],W);DSDPCHKERR(info);
    info = SDPConeVecNorm2(W,&phi);DSDPCHKERR(info);
    if (phi!=phi){ *maxstep = 0.0;  return 0;} 
    if (i>0){
      tt=-subdiag[i-1];
      info = SDPConeVecAXPY(tt,Q[1],W);DSDPCHKERR(info);
    }
    info = SDPConeVecDot(W,Q[0],&tt);DSDPCHKERR(info);
    diag[i]=tt;
    tt*=-1.0;
    info = SDPConeVecAXPY(tt,Q[0],W);DSDPCHKERR(info);
    info = SDPConeVecNorm2(W,&wnorm);DSDPCHKERR(info);
    if (wnorm <= 1.0 * phi){
      for (j=0;j<=i;j++){
	if (j==i-1){
	  info = SDPConeVecDot(W,Q[1],&tt);DSDPCHKERR(info);
	  if (tt==tt){tt*=-1.0;} else {tt=0;}
	  info = SDPConeVecAXPY(tt,Q[1],W);DSDPCHKERR(info);
	  subdiag[i-1]-=tt;
	} else if (j==i){
	  info = SDPConeVecDot(W,Q[0],&tt);DSDPCHKERR(info);
	  if (tt==tt){tt*=-1.0;} else {tt=0;}
	  info = SDPConeVecAXPY(tt,Q[0],W);DSDPCHKERR(info);
	  diag[i]-=tt; 
	}

      } 
    }
    
    info = SDPConeVecNorm2(W,&wnorm);DSDPCHKERR(info);
    /*    printf("PHI: %4.4e, VNORM: %4.2e Diag: %4.2e\n",phi,wnorm,diag[i]); */
    if (i<m-1){
      subdiag[i]=wnorm;
    }
    if (fabs(wnorm)<=1.0e-10){i++;break;}
    info=SDPConeVecCopy(Q[0],Q[1]);DSDPCHKERR(info);
    info=SDPConeVecCopy(W,Q[0]);DSDPCHKERR(info);
    info=SDPConeVecNormalize(Q[0]); DSDPCHKERR(info);
    
  }
  
  /*  DSDPEventLogBegin(id1); */
  info=DSDPGetTriDiagonalEigs(m,diag,subdiag,ddwork,iwork);  DSDPCHKERR(info);
  /*  DSDPEventLogEnd(id1); */
  if (N==0){
    lambda1=-0.0;
    delta=1.0e-20;
    *mineig=0;
  } else if (N==1){
    lambda1=-diag[0];
    delta=1.0e-20;
    *mineig=diag[0];
  } else if (N>1){
    lambda1=-diag[N-1];
    lambda2=-diag[N-2];

    res1=1.0e-8;
    res2=1.0e-8;
     
    tt = -lambda1 + lambda2 - res2;
    if (tt>0) beta=tt;
    else beta=1.0e-20;
    delta = DSDPMin(res1,sqrt(res1)/beta);
    
    *mineig=diag[0];
  }

  
  if (delta-lambda1>0)
    *maxstep = 1.0/(delta-lambda1);
  else
    *maxstep = 1.0e+30;

  info=SDPConeVecGetSize(W,&n);DSDPCHKERR(info);
  DSDPLogInfo(0,19,"Step Length: Fast Lanczos Iterates: %2d, Max: %d, Block Size: %d, VNorm: %3.1e, Lambda1: %4.4e, Lambda2: %4.4e, Delta: %4.2e, Maxstep: %4.2e\n",
	      i,m,n,wnorm,lambda1,lambda2,delta,*maxstep);

  DSDPFunctionReturn(0); 
}
