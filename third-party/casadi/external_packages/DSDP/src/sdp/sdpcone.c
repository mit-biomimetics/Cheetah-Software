#include "dsdpsdp.h"
#include "dsdpsys.h"
/*!
\file sdpcone.c
\brief Implement operations on the SDPCone object.
*/

#undef __FUNCT__
#define __FUNCT__ "SDPConeComputeSS"
/*!
\fn int SDPConeComputeSS(SDPCone sdpcone, int blockj, DSDPVec Y, DSDPVMat SS);
\brief Sum the data matrices.
\param sdpcone semidefinite cone object
\param blockj block number
\param Y scalar multiples of A matrices.
\param SS destination array.
*/
int SDPConeComputeSS(SDPCone sdpcone, int blockj, DSDPVec Y, DSDPVMat SS){
  int info;
  DSDPFunctionBegin;
  info=DSDPVMatZeroEntries(SS); DSDPCHKBLOCKERR(blockj,info);
  info=DSDPBlockASum(&sdpcone->blk[blockj].ADATA,1,Y,SS); DSDPCHKBLOCKERR(blockj,info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SDPConeComputeS"
/*!
\fn int SDPConeComputeS(SDPCone sdpcone, int blockj, double cc,double y[], int nvars, double r, int n, double s[], int nn);
\ingroup SDPRoutines
\brief Compute the dual matrix S
\param sdpcone semidefinite cone object
\param blockj block number
\param cc the multiple of the matrix C (or A_0)
\param y an array of containing the variables y
\param nvars the length of the array and the number of variables y
\param r the multiple of the identity matrix to add
\param n the dimension of the block
\param s array to where the matrix S will be copied.
\param nn length of the array s 
*/
int SDPConeComputeS(SDPCone sdpcone, int blockj, double cc,double y[], int nvars, double r, int n, double s[], int nn){
  int i,info;
  char UPLQ;
  DSDPVMat T;
  DSDPVec Y=sdpcone->Work;
  DSDPFunctionBegin;
  info=SDPConeCheckN(sdpcone,blockj,n);DSDPCHKBLOCKERR(blockj,info);
  info=SDPConeCheckM(sdpcone,nvars);DSDPCHKERR(info);
  if (n<1){DSDPFunctionReturn(0);}
  info=DSDPVecSetC(Y,-1.0*cc);
  info=DSDPVecSetR(Y,-r);
  for (i=0;i<nvars;i++){info=DSDPVecSetElement(Y,i+1,y[i]);}
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPMakeVMatWithArray(UPLQ,s,nn,n,&T);DSDPCHKBLOCKERR(blockj,info);
  info=SDPConeComputeSS(sdpcone,blockj,Y,T);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVMatDestroy(&T);DSDPCHKBLOCKERR(blockj,info);
  DSDPFunctionReturn(0); 
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeAddADotX"
/*!
\fn int SDPConeAddADotX(SDPCone sdpcone, int blockj, double alpha,double x[], int nn, double adotx[], int m);
\ingroup SDPRoutines
\brief Compute the inner products of a dense matrix X with the data matrices.
\param sdpcone semidefinite cone
\param blockj block number
\param alpha multiply inner product be this multiple
\param x the dense array representing a matrix X.
\param nn length of the array
\param adotx array.
\param m the length of the array and the number of variables y, plus two.
*/
int SDPConeAddADotX(SDPCone sdpcone, int blockj, double alpha,double x[], int nn, double adotx[], int m){
  int info,n;
  char UPLQ;
  SDPblk *blk=sdpcone->blk;
  double scl=blk[blockj].ADATA.scl;
  DSDPVec ADOTX,YW2;
  DSDPVMat T;
  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=SDPConeCheckM(sdpcone,m-2);DSDPCHKERR(info);
  YW2=sdpcone->Work2;
  info=DSDPVecSet(alpha,YW2);DSDPCHKBLOCKERR(blockj,info);
  info=SDPConeGetBlockSize(sdpcone,blockj,&n);DSDPCHKBLOCKERR(blockj,info);
  if (n<1){DSDPFunctionReturn(0);}
  info=DSDPVecCreateWArray(&ADOTX,adotx,m);
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPMakeVMatWithArray(UPLQ,x,nn,n,&T);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPBlockADot(&blk[blockj].ADATA,1.0/scl,YW2,T,ADOTX);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVMatDestroy(&T);DSDPCHKBLOCKERR(blockj,info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeComputeXDot"
/*!
\fn int SDPConeComputeXDot(SDPCone sdpcone, int blockj, DSDPVec Y,DSDPVMat X, DSDPVec AX, double* xtrace,double *xnorm, double *tracexs);
\brief Compute inner product of X with the Data, S, and norm of X.
\param sdpcone semidefinite cone
\param blockj block number
\param Y dual solution
\param X dense array matrix
\param AX inner products
\param xtrace trace of X
\param xnorm norm of X
\param tracexs inner product of X and S
 */
int SDPConeComputeXDot(SDPCone sdpcone, int blockj, DSDPVec Y,DSDPVMat X, DSDPVec AX, double* xtrace,double *xnorm, double *tracexs){
  int info;
  SDPblk *blk=sdpcone->blk;
  DSDPVec YW2=sdpcone->Work2;
  double one=1.0,scl=blk[blockj].ADATA.scl;
  DSDPFunctionBegin;
  info=DSDPVecZero(YW2);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPBlockADot(&blk[blockj].ADATA,-1.0/scl,Y,X,YW2);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVecGetR(YW2,xtrace);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVecSum(YW2,tracexs);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVMatNormF2(X,xnorm);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVecSet(one,YW2);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPBlockADot(&blk[blockj].ADATA,1.0/scl,YW2,X,AX);DSDPCHKBLOCKERR(blockj,info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeComputeX3"
/*!
\fn int SDPConeComputeX3(SDPCone sdpcone, int blockj, double mu, DSDPVec Y,DSDPVec DY,DSDPVMat X);
\brief Compute the matrix X with the given information.
\param sdpcone semidefinite cone
\param blockj block number
\param mu barrier parameter
\param Y dual solution
\param DY Newton direction
\param X destination
\sa DSDPComputeX()
 */
int SDPConeComputeX3(SDPCone sdpcone, int blockj, double mu, DSDPVec Y,DSDPVec DY,DSDPVMat X){
  int info;
  double xshift=1e-12,xscale=1e-12;
  SDPblk *blk=sdpcone->blk;
  DSDPTruth psdefinite1=DSDP_FALSE,psdefinite2=DSDP_FALSE,full;
  DSDPDualMat SS;
  
  DSDPFunctionBegin;
  SS=blk[blockj].SS;  
  info=SDPConeComputeSS(sdpcone,blockj,Y,X);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPDualMatSetArray(SS,X); DSDPCHKBLOCKERR(blockj,info);
  info=DSDPDualMatCholeskyFactor(SS,&psdefinite1); DSDPCHKBLOCKERR(blockj,info);
  if (psdefinite1 == DSDP_FALSE){
    DSDPLogInfo(0,2,"Primal SDP Block %2.0f not PSD\n",blockj);
  }
  info=DSDPDualMatInvert(SS);DSDPCHKBLOCKERR(blockj,info);
  info=SDPConeComputeXX(sdpcone,blockj,DY,mu,SS,X);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPDualMatIsFull(SS,&full);DSDPCHKBLOCKERR(blockj,info);
  psdefinite2=DSDP_FALSE;
  while (full==DSDP_TRUE && psdefinite2==DSDP_FALSE && xscale<2e-1){
    info=DSDPVMatShiftDiagonal(X,xshift); DSDPCHKBLOCKERR(blockj,info);
    info=DSDPVMatScaleDiagonal(X,1.0+xscale); DSDPCHKBLOCKERR(blockj,info);
    DSDPLogInfo(0,10,"VMat: shift diagonal: %4.2e, scale diagonal: %4.2e.\n",xshift,1+xscale);
    info=DSDPDualMatSetArray(SS,X); DSDPCHKBLOCKERR(blockj,info);
    info=DSDPDualMatCholeskyFactor(SS,&psdefinite2); DSDPCHKBLOCKERR(blockj,info);
    xshift*=10;xscale*=10;
  }
  if (full==DSDP_FALSE){
    xshift=1e-12,xscale=1e-10;
    info=DSDPVMatShiftDiagonal(X,xshift); DSDPCHKBLOCKERR(blockj,info);
    info=DSDPVMatScaleDiagonal(X,1.0+xscale); DSDPCHKBLOCKERR(blockj,info);
    DSDPLogInfo(0,10,"XMat: shift diagonal: %4.2e, scale diagonal: %4.2e.\n",xshift,1+xscale);
  }
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "SDPConeComputeX"
/*!
\fn int SDPConeComputeX(SDPCone sdpcone, int blockj, int n, double x[], int nn);
\brief Compute the matrix X
\ingroup SDPRoutines
\param sdpcone semidefinite cone
\param n di
\param blockj block number
\param n the dimension of the block
\param x array to where the matrix X will be copied.
\param nn length of the array x 
\sa DSDPComputeX()
 */
int SDPConeComputeX(SDPCone sdpcone, int blockj, int n, double x[], int nn){
  int info;
  double mu=sdpcone->xmakermu;
  double xnorm,xtrace,trxs;
  char UPLQ;
  DSDPVec DY=sdpcone->DYX,Y=sdpcone->YX,AX=sdpcone->Work;
  DSDPVMat T;

  DSDPFunctionBegin;
  info=SDPConeCheckN(sdpcone,blockj,n);DSDPCHKBLOCKERR(blockj,info);
  if (n<1){DSDPFunctionReturn(0);}
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPMakeVMatWithArray(UPLQ,x,nn,n,&T);DSDPCHKBLOCKERR(blockj,info);
  info=SDPConeComputeX3(sdpcone,blockj,mu,Y,DY,T);DSDPCHKBLOCKERR(blockj,info);
  info=SDPConeComputeXDot(sdpcone,blockj,Y,T,AX,&xtrace,&xnorm,&trxs);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVMatDestroy(&T);DSDPCHKBLOCKERR(blockj,info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeViewX"
/*!
\fn int SDPConeViewX(SDPCone sdpcone, int blockj, int n, double x[], int nn);
\brief Print a dense array X to the screen
\ingroup SDPBasic
\param sdpcone semidefinite cone
\param blockj block number
\param n the dimension of the block
\param x dense matrix array.
\param nn length of the array x 
\sa SDPConeGetXArray()
 */
int SDPConeViewX(SDPCone sdpcone, int blockj, int n, double x[], int nn){
  int info;
  char UPLQ;
  DSDPVMat T;

  DSDPFunctionBegin;
  info=SDPConeCheckN(sdpcone,blockj,n);DSDPCHKBLOCKERR(blockj,info);
  info=SDPConeGetStorageFormat(sdpcone,blockj,&UPLQ);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPMakeVMatWithArray(UPLQ,x,nn,n,&T);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVMatView(T);DSDPCHKBLOCKERR(blockj,info);
  info=DSDPVMatDestroy(&T);DSDPCHKBLOCKERR(blockj,info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeXVMultiply"
/*!
\fn int SDPConeXVMultiply(SDPCone sdpcone, int blockj, double vin[], double vout[], int n);
\brief Multiply an array by a factor V such that \f$ V V^T= X \f$.
\ingroup SDPRoutines
\param sdpcone the SDP cone
\param blockj block number
\param vin input array 
\param vout the product V vin
\param n dimension of the block and the length of the two arrays.
\sa SDPConeComputeXV()
\sa SDPConeAddXVAV()
*/
int SDPConeXVMultiply(SDPCone sdpcone, int blockj, double vin[], double vout[], int n){
  int info;
  SDPblk *blk=sdpcone->blk;
  SDPConeVec V1,V2,V3,V4;
  DSDPDualMat S,SS;

  DSDPFunctionBegin;
  info=SDPConeCheckN(sdpcone,blockj,n);DSDPCHKBLOCKERR(blockj,info);
  if (sdpcone->blk[blockj].n>1){
    S=blk[blockj].S;  SS=blk[blockj].SS;
    V2=blk[blockj].W; V3=blk[blockj].W2;
    info=SDPConeVecCreateWArray(&V1,vin,n);
    info=SDPConeVecCreateWArray(&V4,vout,n);
    if (0){
      info=DSDPDualMatCholeskySolveForward(S,V1,V3);DSDPCHKERR(info);
      info=SDPConeVecScale(sqrt(sdpcone->xmakermu),V3);DSDPCHKERR(info);
      info=DSDPDualMatCholeskySolveBackward(S,V3,V2);DSDPCHKERR(info);
      info=DSDPDualMatCholeskyBackwardMultiply(SS,V2,V1);DSDPCHKERR(info);
    }
    info=DSDPDualMatCholeskyForwardMultiply(SS,V1,V2);DSDPCHKERR(info);
    info=DSDPDualMatCholeskySolveForward(S,V2,V3);DSDPCHKERR(info);
    info=SDPConeVecScale(sqrt(sdpcone->xmakermu),V3);DSDPCHKERR(info);
    info=DSDPDualMatCholeskySolveBackward(S,V3,V4);DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeAddXVAV"
/*!
\fn int SDPConeAddXVAV(SDPCone sdpcone, int blockj, double vin[], int n, double sum[], int mm);
\brief Compute \f$ v^T A_{i,j} v \f$ for i = 0 through m.  
\ingroup SDPRoutines
\param sdpcone the SDP cone
\param blockj block number
\param vin array 
\param n dimension of the block and the length of the two arrays.
\param sum input array 
\param mm the number of variables plus 2
\sa SDPConeAddADotX()
 */
int SDPConeAddXVAV(SDPCone sdpcone, int blockj, double vin[], int n, double sum[], int mm){
  int info;
  SDPblk *blk=sdpcone->blk;
  SDPConeVec V2;
  DSDPVec W,Wout;
  DSDPFunctionBegin;
  info=SDPConeCheckN(sdpcone,blockj,n);DSDPCHKBLOCKERR(blockj,info);
  W=sdpcone->Work;
  info=DSDPVecSet(1.0,sdpcone->Work);DSDPCHKBLOCKERR(blockj,info);
  if (sdpcone->blk[blockj].n>1){
    info=SDPConeVecCreateWArray(&V2,vin,n);DSDPCHKERR(info);
    info=DSDPVecCreateWArray(&Wout,sum,mm);DSDPCHKERR(info);
    info=DSDPBlockvAv(&blk[blockj].ADATA,1.0,sdpcone->Work,V2,Wout);DSDPCHKBLOCKERR(blockj,info);
  }
  DSDPFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "SDPConeComputeXV"
/*!
\fn int SDPConeComputeXV(SDPCone sdpcone, int blockj, int *derror);
\ingroup SDPRoutines
\brief Compute a factor V such that \f$ V V^T= X \f$.
\param sdpcone the SDP cone
\param blockj block number
\param derror nonzero if no such factor could be computed.
\sa SDPConeXVMultiply()

This routine is helpful in semidefinite relaxations of combinatorial problems.

*/
int SDPConeComputeXV(SDPCone sdpcone, int blockj, int *derror){

  int info; double rr;
  DSDPVec Y,DY,W;
  SDPblk *blk=sdpcone->blk;
  DSDPTruth psdefinite1=DSDP_FALSE,psdefinite2=DSDP_FALSE;
  DSDPDualMat S,SS;
  DSDPVMat T;

  DSDPFunctionBegin;
  *derror=0;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKBLOCKERR(blockj,info);
  if (sdpcone->blk[blockj].n>1){
    Y=sdpcone->YX; DY=sdpcone->DYX; W=sdpcone->Work;
    T=blk[blockj].T;  S=blk[blockj].S;  SS=blk[blockj].SS;      
    info=DSDPVecWAXPY(W,-1.0,DY,Y);DSDPCHKBLOCKERR(blockj,info);

    while (psdefinite1==DSDP_FALSE){
      info=DSDPVecGetR(W,&rr);
      info=DSDPVecSetR(W,10*rr-1e-12);
      info=SDPConeComputeSS(sdpcone,blockj,W,T);DSDPCHKBLOCKERR(blockj,info);
      info=DSDPDualMatSetArray(SS,T); DSDPCHKBLOCKERR(blockj,info);
      info=DSDPDualMatCholeskyFactor(SS,&psdefinite1); DSDPCHKBLOCKERR(blockj,info);
    }

    while (psdefinite2==DSDP_FALSE){
      info=SDPConeComputeSS(sdpcone,blockj,Y,T);DSDPCHKBLOCKERR(blockj,info);
      info=DSDPDualMatSetArray(S,T); DSDPCHKBLOCKERR(blockj,info);
      info=DSDPDualMatCholeskyFactor(S,&psdefinite2); DSDPCHKBLOCKERR(blockj,info);    
      if (psdefinite2==DSDP_FALSE){
	info=DSDPVecGetR(Y,&rr);
	info=DSDPVecSetR(Y,10*rr-1e-15);
      }
    }
    if (psdefinite1==DSDP_FALSE || psdefinite2==DSDP_FALSE) *derror=1;
  }
  DSDPFunctionReturn(0);
}
