#include "dsdpsdp.h"
#include "dsdpsys.h"
/*!
\file sdpcompute.c
\brief Compute the gradient vector and Hessian matrix. Also compute X matrices
 */

/*
static int tt1=0,tt2=0;
      DSDPEventLogBegin(tt1);
      DSDPEventLogBegin(tt2);
      DSDPEventLogEnd(tt2);
      DSDPEventLogEnd(tt1);
  if (tt1==0){DSDPEventLogRegister("SDP T1 Event",&tt1);}
  if (tt2==0){DSDPEventLogRegister("SDP T2 Event",&tt2);}
*/
#undef __FUNCT__  
#define __FUNCT__ "SDPConeComputeHessian"
/*!
\fn int SDPConeComputeHessian( SDPCone sdpcone, double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2);

\brief Compute the Hessian to the barrier term.
\param sdpcone cone
\param mu barrier parameter
\param M Schur matrix to insert elements.
\param vrhs1 dual objectvive gradient.
\param vrhs2 barrier gradient

*/
int SDPConeComputeHessian( SDPCone sdpcone, double mu, DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2){

  int i,k,kt,kk,m,n,rank,info;
  int ncols,ii,idA;
  double rtemp,ack,ggamma,bmu,scl;
  double rhs1i,rhs2i;
  DSDPTruth method1;
  SDPConeVec W,W2;
  DSDPVec MRowI=sdpcone->Work, Select=sdpcone->Work2;
  DSDPDataMat AA;
  DSDPDualMat S;
  DSDPVMat T;
  DSDPDataTranspose ATranspose=sdpcone->ATR;
  SDPblk *blk=sdpcone->blk;
  DSDPIndex IS;

  /* Evaluate M */
  DSDPFunctionBegin;
  SDPConeValid(sdpcone);
  info=DSDPVecGetSize(vrhs1,&m);DSDPCHKERR(info);

  for (i=0; i<m; i++){  /* One row at a time */

    /* Which Coluns */
    rhs1i=0;rhs2i=0;
    info=DSDPVecZero(MRowI);DSDPCHKERR(info);
    info=DSDPSchurMatRowColumnScaling(M,i,Select,&ncols); DSDPCHKERR(info); 
    if (ncols==0){continue;}
 
    for (kt=0; kt<ATranspose.nnzblocks[i]; kt++){ /* Loop over all blocks */
      kk=ATranspose.nzblocks[i][kt];
      idA=ATranspose.idA[i][kt];
      info=DSDPBlockGetMatrix(&blk[kk].ADATA,idA,&ii,&scl,&AA);DSDPCHKBLOCKERR(kk,info);
      if (ii!=i){DSDPSETERR2(8,"Data Transpose Error: var %d does not equal %d.\n",i,ii);}
      info = DSDPDataMatGetRank(AA,&rank,blk[kk].n);DSDPCHKBLOCKERR(kk,info);
      if (rank==0) continue;

      T=blk[kk].T; S=blk[kk].S; W=blk[kk].W; W2=blk[kk].W2;
      n=blk[kk].n; ggamma=blk[kk].gammamu; bmu=blk[kk].bmu; IS=blk[kk].IS;

      method1=DSDP_TRUE;  /* Simple heuristic */
      if (rank==1) method1=DSDP_FALSE;
      if (rank==2 && ncols<=n) method1=DSDP_FALSE;
      if (rank*rank*ncols<=n+1)method1=DSDP_FALSE;
      if (ncols*blk[kk].nnz < n*n/10) method1=DSDP_FALSE;
      if (ncols==1 && i==m-1)method1=DSDP_FALSE;
      if (n<5) method1=DSDP_TRUE;
      if (0==1) method1=DSDP_FALSE;
      if (method1==DSDP_TRUE){info=DSDPVMatZeroEntries(T);DSDPCHKBLOCKERR(kk,info);}
      for (k=0; k<rank; k++){
	
	info=DSDPDataMatGetEig(AA,k,W,IS,&ack); DSDPCHKBLOCKERR(kk,info);
	if (ack==0.0) continue;
	ack*=scl;
	info=DSDPDualMatInverseMultiply(S,IS,W,W2);DSDPCHKBLOCKERR(kk,info);

	/* RHS terms */
	info = SDPConeVecDot(W,W2,&rtemp); DSDPCHKBLOCKERR(kk,info);
	if (rtemp==0.0) continue;
	rhs1i+=rtemp*ack*bmu; rhs2i+=rtemp*ack*ggamma*mu;
	ack*=(ggamma+bmu);

	if (method1==DSDP_TRUE){
	  info=DSDPVMatAddOuterProduct(T,ack*mu,W2);DSDPCHKBLOCKERR(kk,info);
	} else {
	  info=DSDPBlockvAv(&blk[kk].ADATA,ack*mu,Select,W2,MRowI);DSDPCHKBLOCKERR(kk,info);
	} /* End row computations for rank kk of block kk */
 
      }   /* End row computations for all of block kk     */

      if (method1==DSDP_TRUE){
	info=DSDPBlockADot(&blk[kk].ADATA,1.0,Select,T,MRowI);DSDPCHKBLOCKERR(kk,info);
      }   /* End row computations for all of block ll     */
    }     /* End row computations for all blocks          */
    info=DSDPVecAddElement(vrhs1,i,rhs1i);DSDPCHKERR(info);
    info=DSDPVecAddElement(vrhs2,i,rhs2i);DSDPCHKERR(info);
    info=DSDPSchurMatAddRow(M,i,1.0,MRowI);DSDPCHKERR(info);
  }

  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeComputeRHS"
/*!
\fn int SDPConeComputeRHS( SDPCone sdpcone, int blockj, double mu, DSDPVec vrow,  DSDPVec vrhs1, DSDPVec vrhs2);

\brief Compute the gradient to the barrier term.
\param sdpcone semidefinite cone
\param blockj block of the cone.
\param mu barrier parameter
\param vrow scalar to multiply each element of gradient.
\param vrhs1 dual objectvive gradient.
\param vrhs2 barrier gradient
*/
int SDPConeComputeRHS( SDPCone sdpcone, int blockj, double mu, DSDPVec vrow, DSDPVec vrhs1, DSDPVec vrhs2){

  int info,i,ii,k,rank,nnzmats;
  double dtmp,dyiscale=1,ack,scl,rtemp;
  SDPblk *sdp=&sdpcone->blk[blockj];
  SDPConeVec W=sdp->W,W2=sdp->W2;
  DSDPDataMat AA;
  DSDPVMat T=sdp->T;
  DSDPDualMat S=sdp->S;
  DSDPTruth method1;
  DSDPIndex IS=sdp->IS;

  DSDPFunctionBegin;

  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  method1=DSDP_TRUE;
  method1=DSDP_FALSE;

  if (method1==DSDP_TRUE){
    info=DSDPBlockCountNonzeroMatrices(&sdp->ADATA,&nnzmats);DSDPCHKERR(info);
    for (i=0;i<nnzmats;i++){
      info=DSDPBlockGetMatrix(&sdp->ADATA,i,&ii,&scl,&AA);DSDPCHKERR(info);
      info=DSDPVecGetElement(vrow,ii,&dyiscale);DSDPCHKVARERR(ii,info);
      if (dyiscale==0) continue;
      info=DSDPDataMatGetRank(AA,&rank,sdp->n);DSDPCHKVARERR(ii,info);
      for (k=0; k<rank; k++){
	info=DSDPDataMatGetEig(AA,k,W,IS,&ack); DSDPCHKVARERR(ii,info);
	if (ack==0) continue;
	info=DSDPDualMatInverseMultiply(S,IS,W,W2);DSDPCHKVARERR(ii,info);
	info=SDPConeVecDot(W,W2,&rtemp); DSDPCHKVARERR(ii,info);
	dtmp=rtemp*ack*mu*dyiscale*scl;
	info=DSDPVecAddElement(vrhs2,ii,dtmp);DSDPCHKVARERR(ii,info);
      }
    }
    
  } else {
    info=DSDPVMatZeroEntries(T);DSDPCHKERR(info);
    info=DSDPDualMatInverseAdd(S,mu,T);DSDPCHKERR(info);
    info=DSDPBlockADot(&sdp->ADATA,1.0,vrow,T,vrhs2);DSDPCHKERR(info);
  }
  
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeMultiply"
/*!
\fn int SDPConeMultiply( SDPCone sdpcone, int blockj, double mu, DSDPVec vrow,  DSDPVec vin, DSDPVec vout);

\brief Compute the gradient to the barrier term.
\param sdpcone semidefinite cone
\param blockj block of cone
\param mu barrier parameter
\param vrow scalar to multiply each element of the product
\param vin in vector.
\param vout product
*/
int SDPConeMultiply(SDPCone sdpcone, int blockj, double mu, DSDPVec vrow, DSDPVec vin, DSDPVec vout){

  int info,i,ii,k,rank,nnzmats;
  double dd2,dtmp,dyiscale,ack,scl,vv;
  SDPblk *sdp=&sdpcone->blk[blockj];
  SDPConeVec W=sdp->W,W2=sdp->W2;
  DSDPDataMat AA;
  DSDPVMat T=sdp->T;
  DSDPDSMat DS=sdp->DS;
  DSDPIndex IS=sdp->IS;
  DSDPDualMat S=sdp->S;

  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  info=DSDPVMatZeroEntries(T); DSDPCHKERR(info);
  info=DSDPBlockASum(&sdp->ADATA,-1.0,vin,T); DSDPCHKERR(info);
  info=DSDPDSMatSetArray(DS,T); DSDPCHKERR(info);
  info=DSDPBlockCountNonzeroMatrices(&sdp->ADATA,&nnzmats);DSDPCHKERR(info);
  for (i=0;i<nnzmats;i++){
    info=DSDPBlockGetMatrix(&sdp->ADATA,i,&ii,&scl,&AA);DSDPCHKERR(info);
    info=DSDPVecGetElement(vrow,ii,&dyiscale);DSDPCHKVARERR(ii,info);
    if (dyiscale==0) continue;

    info=DSDPDataMatGetRank(AA,&rank,sdp->n);DSDPCHKVARERR(ii,info);
    for (dd2=0,k=0; k<rank; k++){
      info=DSDPDataMatGetEig(AA,k,W,IS,&ack); DSDPCHKVARERR(ii,info);
      if (ack==0) continue;
      info=DSDPDualMatInverseMultiply(S,IS,W,W2);DSDPCHKVARERR(ii,info);
      info=DSDPDSMatVecVec(DS,W2,&vv);DSDPCHKVARERR(ii,info);
      dd2+=vv*ack;
    }
    dtmp = dd2 * dyiscale *mu *scl;
    info=DSDPVecAddElement(vout,ii,dtmp);DSDPCHKVARERR(ii,info);
  }

  DSDPFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "SDPConeComputeXX"
/*!
\fn int SDPConeComputeXX( SDPCone sdpcone, int blockj, DSDPVec DY, double mu, DSDPDualMat S, DSDPVMat X);

\brief Compute X
\param sdpcone cone
\param blockj block number in SDP cone.
\param DY step direction vector
\param mu barrier parameter
\param S dual matrix alread inverted.
\param X the result.
*/
int SDPConeComputeXX(SDPCone sdpcone, int blockj, DSDPVec DY, double mu,DSDPDualMat S, DSDPVMat X){
  
  int info, i, ii,k, rank, n, nnzmats;
  double dtmp,dyiscale,ack,scl;
  SDPblk *sdp=&sdpcone->blk[blockj];
  SDPConeVec W=sdp->W,W2=sdp->W2;
  DSDPDataMat AA;
  DSDPIndex IS=sdp->IS;

  DSDPFunctionBegin;
  info=SDPConeCheckJ(sdpcone,blockj);DSDPCHKERR(info);
  mu=mu*sdp->gammamu; n=sdp->n;
  info=DSDPVMatZeroEntries(X);DSDPCHKERR(info);
  info=DSDPBlockCountNonzeroMatrices(&sdp->ADATA,&nnzmats);DSDPCHKERR(info);
  for (i=0;i<nnzmats;i++){
    info=DSDPBlockGetMatrix(&sdp->ADATA,i,&ii,&scl,&AA);DSDPCHKVARERR(ii,info);
    info=DSDPVecGetElement(DY,ii,&dyiscale);DSDPCHKVARERR(ii,info);
    if (dyiscale==0) continue;
    info=DSDPDataMatGetRank(AA,&rank,n);DSDPCHKVARERR(ii,info);
    for (k=0; k<rank; k++){
      info = DSDPDataMatGetEig(AA,k,W,IS,&ack);DSDPCHKVARERR(ii,info);
      if (ack==0) continue;
      info=DSDPDualMatInverseMultiply(S,IS,W,W2);DSDPCHKVARERR(ii,info);
      dtmp = ack * dyiscale * mu * scl;
      info=DSDPVMatAddOuterProduct(X,dtmp,W2);DSDPCHKVARERR(ii,info);
    }
  }

  info=DSDPDualMatInverseAdd(S,mu,X);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

