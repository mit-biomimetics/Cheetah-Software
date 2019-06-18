#include "dsdpsdp.h"
#include "dsdpsys.h"
/*!
\file sdpconesetup.c
\brief Setup the internal data structures needed by the SDPCone object 
*/

#undef __FUNCT__  
#define __FUNCT__ "DSDPDataTransposeInitialize"
/*!
\fn int DSDPDataTransposeInitialize(DSDPDataTranspose *ATranspose);
\brief Initialize transpose structure for data.
\param ATranspose transpose structure for data.
*/
int DSDPDataTransposeInitialize(DSDPDataTranspose *ATranspose){
  DSDPFunctionBegin;
  ATranspose->nnzblocks=0;
  ATranspose->nzblocks=0;
  ATranspose->idA=0;
  ATranspose->idAP=0;
  ATranspose->ttnzmat=0;
  ATranspose->nnzblocks=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DSDPDataTransposeSetup"
/*!
\fn int DSDPDataTransposeSetup(DSDPDataTranspose *ATranspose, SDPblk *blk, int nblocks, int m);
\brief Set up transpose structure for data.
\param ATranspose transpose structure for data.
\param blk semidefinite blocks
\param nblocks number of blocks
\param m dimension of Y vector.
*/
int DSDPDataTransposeSetup(DSDPDataTranspose *ATranspose, SDPblk *blk, int nblocks, int m){
  
  int i,ii,kk,vvar,info;
  int nnzmats,tnzmats=0;
  DSDPFunctionBegin;

  info=DSDPDataTransposeTakeDown(ATranspose);DSDPCHKERR(info);
  /* Determine sparsity pattern of SDP Data Matrices */

  DSDPCALLOC2(&ATranspose->nnzblocks,int,(m),&info);DSDPCHKERR(info);
  DSDPCALLOC2(&ATranspose->nzblocks,int*,(m),&info);DSDPCHKERR(info);
  DSDPCALLOC2(&ATranspose->idA,int*,(m),&info);DSDPCHKERR(info);
  ATranspose->m=m;
  for (i=0;i<m;i++){ ATranspose->nnzblocks[i]=0; }
  for (kk=0; kk<nblocks; kk++){
    info=DSDPBlockDataMarkNonzeroMatrices(&blk[kk].ADATA,ATranspose->nnzblocks);DSDPCHKERR(info);
  }
  for (tnzmats=0,i=0;i<m;i++){ tnzmats += ATranspose->nnzblocks[i];}

  DSDPCALLOC2(&ATranspose->ttnzmat,int,tnzmats,&info);DSDPCHKERR(info);
  ATranspose->nzblocks[0]=ATranspose->ttnzmat;
  for (i=1;i<m;i++){
    ATranspose->nzblocks[i]=ATranspose->nzblocks[i-1]+ATranspose->nnzblocks[i-1];
  }

  DSDPCALLOC2(&ATranspose->idAP,int,tnzmats,&info);DSDPCHKERR(info);
  ATranspose->idA[0]=ATranspose->idAP;
  for (i=1;i<m;i++){
    ATranspose->idA[i]=ATranspose->idA[i-1]+ATranspose->nnzblocks[i-1];
  }

  for (i=0;i<m;i++){ATranspose->nnzblocks[i]=0;}
  for (kk=0; kk<nblocks; kk++){
    info=DSDPBlockCountNonzeroMatrices(&blk[kk].ADATA,&nnzmats);DSDPCHKERR(info);
    for (i=0;i<nnzmats;i++){
      info=DSDPBlockGetMatrix(&blk[kk].ADATA,i,&ii,0,0);DSDPCHKERR(info);
      vvar=ATranspose->nnzblocks[ii];
      ATranspose->nzblocks[ii][vvar]=kk;
      ATranspose->idA[ii][vvar]=i;
      ATranspose->nnzblocks[ii]++;
    }
  }

  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPDataTransposeTakeDown"
/*!
\fn int DSDPDataTransposeTakeDown(DSDPDataTranspose *ATranspose);
\brief Free transpose structure for data.
\param ATranspose transpose structure for data.
*/
int DSDPDataTransposeTakeDown(DSDPDataTranspose *ATranspose){
  int info;
  DSDPFunctionBegin;
  DSDPFREE(&ATranspose->ttnzmat,&info);DSDPCHKERR(info);
  DSDPFREE(&ATranspose->idAP,&info);DSDPCHKERR(info);
  DSDPFREE(&ATranspose->nzblocks,&info);DSDPCHKERR(info);
  DSDPFREE(&ATranspose->nnzblocks,&info);DSDPCHKERR(info);
  DSDPFREE(&ATranspose->idA,&info);DSDPCHKERR(info);
  info=DSDPDataTransposeInitialize(ATranspose);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPCreateSDPCone"
/*!
\fn int DSDPCreateSDPCone(DSDP dsdp, in nblocks, SDPCone *sdpcone);
\brief Create a semidefinite matrix cone with one or more blocks.
\ingroup SDPBasic
\param dsdp the solver
\param nblocks number of blocks
\param sdpcone the new semidefinite cone with multiple blocks
\sa DSDPCreate()
\sa SDPConeSetASparseVecMat()
*/
int DSDPCreateSDPCone(DSDP dsdp, int blocks, SDPCone* dspcone){
  int i,info;
  SDPCone sdpcone;

  DSDPFunctionBegin;
  DSDPCALLOC1(&sdpcone,struct SDPCone_C,&info);DSDPCHKERR(info);
  *dspcone=sdpcone;
  sdpcone->keyid=SDPCONEKEY;
  info=DSDPAddSDP(dsdp,sdpcone);DSDPCHKERR(info);

  info=DSDPGetNumberOfVariables(dsdp,&sdpcone->m);DSDPCHKERR(info);
  DSDPCALLOC2(&sdpcone->blk,SDPblk,blocks,&info); DSDPCHKERR(info);
  for (i=0;i<blocks; i++){
    info=DSDPBlockInitialize(&sdpcone->blk[i]); DSDPCHKBLOCKERR(i,info);
  }

  sdpcone->nblocks=blocks;
  sdpcone->optype=3;
  info=DSDPUseDefaultDualMatrix(sdpcone); DSDPCHKERR(info);

  sdpcone->nn=0;
  sdpcone->dsdp=dsdp;
  info=DSDPDataTransposeInitialize(&sdpcone->ATR); DSDPCHKERR(info);
  info=DSDPBlockEventZero();DSDPCHKERR(info);
  info=DSDPDualMatEventZero();DSDPCHKERR(info);
  info=DSDPVMatEventZero();DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


int DSDPCreateS(DSDPBlockData*,char,int,DSDPVec,DSDPVMat,SDPConeVec,SDPConeVec,DSDPDualMat*, DSDPDualMat*, DSDPDSMat*, void*);

#undef __FUNCT__ 
#define __FUNCT__ "DSDPBlockSetup"
/*!
\fn int DSDPBlockSetup(SDPblk *blk, int blockj, DSDPVec WY);
\brief Allocate data structures of one block the cone.
\param blk block in semidefinite cone
\param blockj block number
\param WY sample variable vector.
*/
int DSDPBlockSetup(SDPblk *blk, int blockj, DSDPVec WY){ 
  int n,info,trank,flag;
  DSDPFunctionBegin;
  /*  
  info=DSDPBlockTakeDown(blk);DSDPCHKERR(info);
  */
  n=blk->n;
  info=DSDPVMatExist(blk->T,&flag);DSDPCHKERR(info);
  if (flag==0){
    info=DSDPMakeVMat(blk->format,n,&blk->T);DSDPCHKERR(info);
  }

  info = DSDPIndexCreate(blk->n,&blk->IS);DSDPCHKERR(info);
  info = SDPConeVecCreate(blk->n,&blk->W);DSDPCHKERR(info);
  info = SDPConeVecDuplicate(blk->W,&blk->W2);DSDPCHKERR(info);

  /* Build Lanczos structures */
  info=DSDPSetMaximumLanczosIterations(&blk->Lanczos,20); DSDPCHKERR(info);
  if (n>30){info=DSDPSetMaximumLanczosIterations(&blk->Lanczos,20); DSDPCHKERR(info);}
  if (n>300){info=DSDPSetMaximumLanczosIterations(&blk->Lanczos,40); DSDPCHKERR(info);}
  if (n>1000){info=DSDPSetMaximumLanczosIterations(&blk->Lanczos,50); DSDPCHKERR(info);}
  if (1){
    info=DSDPFastLanczosSetup(&blk->Lanczos,blk->W);DSDPCHKERR(info);
    DSDPLogInfo(0,19,"SDP Block %d using Fast Lanczos\n",blockj);
  } else {
    info=DSDPRobustLanczosSetup(&blk->Lanczos,blk->W);DSDPCHKERR(info);
    DSDPLogInfo(0,19,"SDP Block %d using Full Lanczos\n",blockj);
  }

  /* Eigenvalues and Eigenvectors */
  info=DSDPBlockFactorData(&blk->ADATA,blk->T,blk->W);DSDPCHKERR(info);
  info=DSDPBlockDataRank(&blk->ADATA,&trank,n);DSDPCHKERR(info);

  info=DSDPCreateS(&blk->ADATA,blk->format,trank,WY,blk->T,blk->W,blk->W2,&blk->S,&blk->SS,&blk->DS,0);DSDPCHKERR(info);
  
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeBlockNNZ"
int SDPConeBlockNNZ(SDPblk *blk,int m){
  int i,ii,n,info,nnz,nnzmats,tnnzmats,tnnz=0;
  double scl;
  DSDPDataMat AA;
  DSDPFunctionBegin;
  nnzmats=blk->ADATA.nnzmats;tnnzmats=nnzmats;
  n=blk->n;

  for (i=0;i<nnzmats;i++){
    info=DSDPBlockGetMatrix(&blk->ADATA,i,&ii,&scl,&AA);DSDPCHKERR(info);
    if (ii==0){tnnzmats--; continue;}
    if (ii==m-1){continue;}
    info = DSDPDataMatCountNonzeros(AA,&nnz,n); DSDPCHKERR(info);
    tnnz+= (nnz*(tnnzmats-i));
  }
  if (tnnzmats>1){ tnnz=tnnz/((tnnzmats)*(tnnzmats+1)/2); }
  if (tnnz<1)  tnnz = 1; 
  blk->nnz=tnnz;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeSetup2"
/*!
\fn int SDPConeSetup2(SDPCone sdpcone, DSDPVec yy0, DSDPSchurMat M);
\brief Allocate data structure of the cone.
\param sdpcone semidefinite cone
\param yy0 variable vector
\param M Schur matrix.
*/
int SDPConeSetup2(SDPCone sdpcone, DSDPVec yy0, DSDPSchurMat M){  
  int kk,n,m,info;
  double nn=0;
  SDPblk *blk;
  DSDPFunctionBegin;
  info=DSDPVecGetSize(yy0,&m);DSDPCHKERR(info);
  for (kk=0; kk<sdpcone->nblocks; kk++){
    blk=&sdpcone->blk[kk];
    n=blk->n;
    info=SDPConeBlockNNZ(blk,m);DSDPCHKERR(info);
    info=DSDPBlockSetup(blk,kk,sdpcone->Work);DSDPCHKERR(info);
    nn+=n*blk->gammamu;
  }
  sdpcone->nn=(int)nn;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeSetup"
/*!
\fn int SDPConeSetup(SDPCone sdpcone, DSDPVec yy0);
\brief Allocate data structure of the cone.
\param sdpcone semidefinite cone
\param yy0 variable vector
*/
int SDPConeSetup(SDPCone sdpcone, DSDPVec yy0){
  int kk,n,m,info;
  DSDPFunctionBegin;

  info = DSDPVecGetSize(yy0,&m);DSDPCHKERR(info);
  if (m!=sdpcone->m+2){DSDPSETERR(8,"CHECK DIMENSION\n");}
  info = DSDPVecDuplicate(yy0,&sdpcone->Work);DSDPCHKERR(info);
  info = DSDPVecDuplicate(yy0,&sdpcone->Work2);DSDPCHKERR(info);
  info = DSDPVecDuplicate(yy0,&sdpcone->YY);DSDPCHKERR(info);
  info = DSDPVecDuplicate(yy0,&sdpcone->YX);DSDPCHKERR(info);
  info = DSDPVecDuplicate(yy0,&sdpcone->DYX);DSDPCHKERR(info);
  for (kk=0; kk<sdpcone->nblocks; kk++){
    n=sdpcone->blk[kk].n;
    info=SDPConeSetRIdentity(sdpcone,kk,n,1.0);DSDPCHKERR(info);
  }
  
  info=DSDPDataTransposeSetup(&sdpcone->ATR,sdpcone->blk,sdpcone->nblocks,m); DSDPCHKERR(info);
  info=DSDPBlockEventInitialize();DSDPCHKERR(info);
  info=DSDPDualMatEventInitialize();DSDPCHKERR(info);
  info=DSDPVMatEventInitialize();DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__ 
#define __FUNCT__ "DSDPBlockInitialize"
/*!
\fn int DSDPBlockInitialize(SDPblk *blk);
\brief Initialize data structures in one block of the cone.
\param blk block of semidefinite cone
*/
int DSDPBlockInitialize(SDPblk *blk){
  int info;
  DSDPFunctionBegin;
  blk->n=0;
  blk->format='N';
  blk->gammamu=1.0;
  blk->bmu=0.0;
  blk->nnz=10000000;

  info = DSDPDualMatInitialize(&blk->S); DSDPCHKERR(info);
  info = DSDPDualMatInitialize(&blk->SS); DSDPCHKERR(info);
  info = DSDPDSMatInitialize(&blk->DS); DSDPCHKERR(info);
  info = DSDPVMatInitialize(&blk->T); DSDPCHKERR(info);
  info = DSDPLanczosInitialize(&blk->Lanczos); DSDPCHKERR(info);
  info = DSDPBlockDataInitialize(&blk->ADATA); DSDPCHKERR(info);
  info = DSDPIndexInitialize(&blk->IS); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}
  
#undef __FUNCT__  
#define __FUNCT__ "DSDPBlockTakeDown"
/*!
\fn int DSDPBlockTakeDown(SDPblk *blk);
\brief Free data structures in one block of the cone.
\param blk block of semidefinite cone
*/
int DSDPBlockTakeDown(SDPblk *blk){
  int info;
  DSDPFunctionBegin;
  if (!blk){DSDPFunctionReturn(0);}
  info=DSDPBlockTakeDownData(&blk->ADATA);DSDPCHKERR(info);
  info=SDPConeVecDestroy(&blk->W);DSDPCHKERR(info);
  info=SDPConeVecDestroy(&blk->W2);DSDPCHKERR(info);
  info=DSDPIndexDestroy(&blk->IS);DSDPCHKERR(info);
  info=DSDPLanczosDestroy(&blk->Lanczos);DSDPCHKERR(info);
  info=DSDPDualMatDestroy(&blk->SS);DSDPCHKERR(info);
  info=DSDPDualMatDestroy(&blk->S);DSDPCHKERR(info);
  info=DSDPDSMatDestroy(&blk->DS);DSDPCHKERR(info);
  info=DSDPVMatDestroy(&blk->T);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPConeTakeDown"
/*!
\fn int DSDPConeTakeDown(SDPCone);
\brief Free data structure of the cone.
\param sdpcone semidefinite cone
*/
int DSDPConeTakeDown(SDPCone sdpcone){
  int blockj,info;
  DSDPFunctionBegin;
  for (blockj=0; blockj<sdpcone->nblocks; blockj++){
    info=DSDPBlockTakeDown(&sdpcone->blk[blockj]);DSDPCHKERR(info);
  }
  info=DSDPVecDestroy(&sdpcone->Work);DSDPCHKERR(info);
  info=DSDPVecDestroy(&sdpcone->Work2);DSDPCHKERR(info);
  info=DSDPVecDestroy(&sdpcone->YY);DSDPCHKERR(info);
  info=DSDPVecDestroy(&sdpcone->YX);DSDPCHKERR(info);
  info=DSDPVecDestroy(&sdpcone->DYX);DSDPCHKERR(info);
  info=DSDPDataTransposeTakeDown(&sdpcone->ATR);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "SDPConeDestroy"
/*!
\fn int SDPConeDestroy(SDPCone sdpcone);
\brief Free data structure of the cone.
\param sdpcone semidefinite cone
*/
int SDPConeDestroy(SDPCone sdpcone){
  int blockj,info;
  DSDPFunctionBegin;
  info=DSDPConeTakeDown(sdpcone);DSDPCHKERR(info);
  for (blockj=0; blockj<sdpcone->nblocks; blockj++){
    info=DSDPBlockDataDestroy(&sdpcone->blk[blockj].ADATA);DSDPCHKERR(info);
  }
  DSDPFREE(&sdpcone->blk,&info);DSDPCHKERR(info);
  DSDPFREE(&sdpcone,&info);DSDPCHKERR(info);
  info=DSDPBlockEventZero();DSDPCHKERR(info);
  info=DSDPDualMatEventZero();DSDPCHKERR(info);
  info=DSDPVMatEventZero();DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

