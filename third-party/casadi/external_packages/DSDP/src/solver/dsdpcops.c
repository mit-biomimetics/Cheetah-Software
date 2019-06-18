#include "dsdp.h"
#include "dsdpsys.h"
/*!
\file dsdpcops.c
\brief Applies conic operations to each cone in the solver.
*/

#define DSDPCHKCONEERR(a,b);  { if (b){ DSDPSETERR1(b,"Cone Number: %d,\n",a);} }

static int ConeSetup=0,ConeComputeS=0,ConeComputeSS=0,ConeComputeH=0,ConeHMultiplyAdd=0,ConeMaxPStep=0,ConeMaxDStep=0,ConePotential=0,ConeComputeX=0,ConeView=0,ConeDestroy=0,ConeXEigs=0,ConeRHS=0,ConeInvertS=0;
static int DSDPRegisterConeEvents(void);
int DSDPSetSchurMatOps(DSDP,struct DSDPSchurMat_Ops*, void*);
int DSDPAddSchurRow(DSDP,int, DSDPVec);
/*
int DSDPIdentifySchurColumns(DSDP,int, int*, int*, int);
*/

#undef __FUNCT__  
#define __FUNCT__ "DSDPZeroConeEvents"
static int DSDPZeroConeEvents(){
  DSDPFunctionBegin;
  ConeSetup=0;ConeComputeS=0;ConeComputeSS=0;ConeComputeH=0;ConeHMultiplyAdd=0;ConeMaxPStep=0;ConeMaxDStep=0;ConePotential=0;ConeComputeX=0;ConeView=0;ConeDestroy=0;ConeXEigs=0;ConeRHS=0;ConeInvertS=0;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPRegisterConeEvents"
static int DSDPRegisterConeEvents(){

  DSDPFunctionBegin;
  DSDPEventLogRegister("Cone Setup 1&2",&ConeSetup);
  DSDPEventLogRegister("Cone Invert S",&ConeInvertS);
  DSDPEventLogRegister("Cone RHS",&ConeRHS);
  DSDPEventLogRegister("Cone Compute Newton Eq.",&ConeComputeH);
  DSDPEventLogRegister("Cone Newton Multiply-Add",&ConeHMultiplyAdd);
  DSDPEventLogRegister("Cone Max P Step Length",&ConeMaxPStep);
  DSDPEventLogRegister("Cone Compute and Factor SP",&ConeComputeSS);
  DSDPEventLogRegister("Cone Max D Step Length",&ConeMaxDStep);
  DSDPEventLogRegister("Cone Compute and Factor S",&ConeComputeS);
  DSDPEventLogRegister("Cone Potential",&ConePotential);
  DSDPEventLogRegister("Cone View",&ConeView);
  DSDPEventLogRegister("Cone Compute X",&ConeComputeX);
  DSDPEventLogRegister("Cone X Residuals",&ConeXEigs);
  DSDPEventLogRegister("Cone Destroy",&ConeDestroy);

  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetUpCones( DSDP dsdp);

\brief Each cone should factor data or allocate internal data structures.
\param dsdp the solver

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetUpCones"
int DSDPSetUpCones(DSDP dsdp){
  int info,kk;
  DSDPVec yy0=dsdp->y;
  DSDPFunctionBegin;
  info=DSDPRegisterConeEvents();
  DSDPEventLogBegin(ConeSetup);
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeSetUp(dsdp->K[kk].cone,yy0);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  DSDPEventLogEnd(ConeSetup);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPSetUpCones2( DSDP dsdp, DSDPVec yy0, DSDPSchurMat M);

\brief Each cone should allocate its data structures .
\param dsdp the solver
\param yy0 variable vector
\param M Shur Matrix

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPSetUpCones2"
int DSDPSetUpCones2(DSDP dsdp, DSDPVec yy0, DSDPSchurMat M){
  int info,kk;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeSetup);
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeSetUp2(dsdp->K[kk].cone,yy0,M);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  DSDPEventLogEnd(ConeSetup);
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPDestroyCones( DSDP dsdp);

\brief Each cone shoudl free its data structures.
\param dsdp the solver

*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPDestroyCones"
int DSDPDestroyCones(DSDP dsdp){
  int info,kk,ncones=dsdp->ncones;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeDestroy);
  for (kk=ncones-1;kk>=0; kk--){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeDestroy(&dsdp->K[kk].cone);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
    info=DSDPConeInitialize(&dsdp->K[kk].cone);DSDPCHKCONEERR(kk,info);
    dsdp->ncones--;
  }
  if (dsdp->maxcones>0){
    DSDPFREE(&dsdp->K,&info);DSDPCHKERR(info);
    dsdp->K=0;
    dsdp->maxcones=0;
  }
  DSDPEventLogEnd(ConeDestroy);
  info=DSDPZeroConeEvents();DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPComputeHessian( DSDP dsdp , DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2);

\brief Compute the Schur complement, or Gram, matrix for each cone.
\param dsdp the solver
\param M matrix
\param vrhs1 gradient of objective (b)
\param vrhs2 gradient of barrier


*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeHessian"
int DSDPComputeHessian( DSDP dsdp , DSDPSchurMat M,  DSDPVec vrhs1, DSDPVec vrhs2){
  int info,kk; double r;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeComputeH);
  dsdp->schurmu=dsdp->mutarget;
  info=DSDPVecGetR(dsdp->y,&r);DSDPCHKERR(info);
  info=DSDPSchurMatSetR(dsdp->M,r);DSDPCHKERR(info);
  info=DSDPSchurMatZeroEntries(M);DSDPCHKERR(info);
  info=DSDPVecZero(vrhs1);DSDPCHKERR(info);
  info=DSDPVecZero(vrhs2);DSDPCHKERR(info);
  info=DSDPVecZero(M.schur->rhs3);DSDPCHKERR(info);
  info=DSDPObjectiveGH(dsdp,M,vrhs1); DSDPCHKERR(info);
  for (kk=dsdp->ncones-1;kk>=0;kk--){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeComputeHessian(dsdp->K[kk].cone,dsdp->schurmu,M,vrhs1,vrhs2);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  info=DSDPSchurMatAssemble(M);DSDPCHKERR(info);
  /*    DSDPSchurMatView(M); */
  info=DSDPSchurMatReducePVec(M,vrhs1);DSDPCHKERR(info);
  info=DSDPSchurMatReducePVec(M,vrhs2);DSDPCHKERR(info);
  info=DSDPSchurMatReducePVec(M,M.schur->rhs3);DSDPCHKERR(info);
  if (0 && dsdp->UsePenalty==DSDPNever){
    info=DSDPVecAXPY(1.0,M.schur->rhs3,vrhs2);DSDPCHKERR(info);
    info=DSDPVecZero(M.schur->rhs3);DSDPCHKERR(info);
    info=DSDPVecZero(M.schur->dy3);DSDPCHKERR(info);
    info=DSDPVecSetR(vrhs1,0);DSDPCHKERR(info);
    info=DSDPVecSetR(vrhs2,r);DSDPCHKERR(info);
  }
  DSDPEventLogEnd(ConeComputeH);
  DSDPFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DSDPHessianMultiplyAdd"
/*!
\fn int DSDPHessianMultiplyAdd( DSDP dsdp , DSDPVec v, DSDPVec vv);

\brief Add the product of Schur matrix with v.
\param dsdp the solver
\param v input vector.
\param vv product gradient of barrier


*/
int DSDPHessianMultiplyAdd( DSDP dsdp , DSDPVec v, DSDPVec vv){
  int info,kk;
  DSDPVec vrow=dsdp->sles->BR;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeHMultiplyAdd);

  info=DSDPSchurMatRowScaling(dsdp->M,vrow);DSDPCHKERR(info);
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeMultiplyAdd(dsdp->K[kk].cone,dsdp->schurmu,vrow,v,vv);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  info=DSDPSchurMatReducePVec(dsdp->M,vv);DSDPCHKERR(info);
  DSDPEventLogEnd(ConeHMultiplyAdd);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeG"
/*!
\fn int DSDPComputeG( DSDP dsdp , DSDPVec vt, DSDPVec vrhs1, DSDPVec vrhs2);
\param dsdp the solver
\param vt scaling for each element in the next two vectors
\param vrhs1 scaled gradient of the objective function
\param vrhs2 scaled gradient of the barrier function
\brief Compute the gradient of the barrier for each cone.
 */
int DSDPComputeG( DSDP dsdp , DSDPVec vt, DSDPVec vrhs1, DSDPVec vrhs2){
  int info,kk; double r;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeRHS);
  info=DSDPVecZero(vrhs1);DSDPCHKERR(info);
  info=DSDPVecZero(vrhs2);DSDPCHKERR(info);
  info=DSDPVecGetR(dsdp->y,&r);DSDPCHKERR(info);
  info=DSDPSchurMatSetR(dsdp->M,r);DSDPCHKERR(info);
  info=DSDPSchurMatRowScaling(dsdp->M,vt);DSDPCHKERR(info);
  info=DSDPObjectiveGH(dsdp,dsdp->M,vrhs1); DSDPCHKERR(info);
  if (0 && r==0){info=DSDPVecSetR(vrhs1,0);info=DSDPVecSetR(vrhs2,0);}
  /*  info=DSDPVecScale(1.0/dsdp->schurmu,vrhs1); DSDPCHKERR(info); */
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeComputeRHS(dsdp->K[kk].cone,dsdp->schurmu,vt,vrhs1,vrhs2);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  DSDPEventLogEnd(ConeRHS);
  info=DSDPSchurMatReducePVec(dsdp->M,vrhs1);DSDPCHKERR(info);
  info=DSDPSchurMatReducePVec(dsdp->M,vrhs2);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeANorm2"
/*!
\fn int DSDPComputeANorm2( DSDP dsdp , DSDPVec Anorm2);
\param dsdp the solver
\param Anorm2 norm of data corresponding to each variable y.
\brief Compute norm of A and C.
 */
int DSDPComputeANorm2( DSDP dsdp , DSDPVec Anorm2){
  int info,kk;
  DSDPFunctionBegin;
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeANorm2(dsdp->K[kk].cone,Anorm2);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPComputeSS(DSDP dsdp, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite);

\brief Compute the dual variables S in each cone.

\param dsdp the solver
\param Y variables
\param flag primal or dual structure
\param ispsdefinite DSDP_TRUE if a member of the cone, DSDP_FALSE otherwise.


*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeSS"
int DSDPComputeSS(DSDP dsdp, DSDPVec Y, DSDPDualFactorMatrix flag, DSDPTruth *ispsdefinite){
  int info,kk;
  DSDPTruth psd=DSDP_TRUE;
  DSDPFunctionBegin;
  if (flag==DUAL_FACTOR){
    DSDPEventLogBegin(ConeComputeS);
  } else if (flag==PRIMAL_FACTOR){
    DSDPEventLogBegin(ConeComputeSS);
  }
  for (kk=dsdp->ncones-1; kk>=0 && psd==DSDP_TRUE;kk--){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeComputeS(dsdp->K[kk].cone,Y,flag,&psd); DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  *ispsdefinite=psd;
  if (flag==DUAL_FACTOR){
    DSDPEventLogEnd(ConeComputeS);
  } else if (flag==PRIMAL_FACTOR){
    DSDPEventLogEnd(ConeComputeSS);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPInvertS(DSDP dsdp);

\brief Invert the S variables in each cone.

\param dsdp the solver

\sa DSDPComputeSS()

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPInvertS"
int DSDPInvertS(DSDP dsdp){
  int info,kk;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeInvertS);
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeInvertS(dsdp->K[kk].cone); DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  DSDPEventLogEnd(ConeInvertS);
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputeMaxStepLength(DSDP dsdp, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength);

\brief Compute the maximum step length for the given step direction.

\param dsdp the solver
\param DY step direction
\param flag primal or dual structure
\param maxsteplength the minumum of maximums on each cone.

\sa DSDPComputeSS()


*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeMaxStepLength"
int DSDPComputeMaxStepLength(DSDP dsdp, DSDPVec DY, DSDPDualFactorMatrix flag, double *maxsteplength){
  int info,kk;
  double msteplength=1.0e30,conesteplength;
  DSDPFunctionBegin;
  if (flag==DUAL_FACTOR){
    DSDPEventLogBegin(ConeMaxDStep);
  } else if (flag==PRIMAL_FACTOR){
    DSDPEventLogBegin(ConeMaxPStep);
  }
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    conesteplength=1.0e20;
    info=DSDPConeComputeMaxStepLength(dsdp->K[kk].cone,DY,flag,&conesteplength);DSDPCHKCONEERR(kk,info);
    msteplength=DSDPMin(msteplength,conesteplength);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  *maxsteplength=msteplength;
  if (flag==DUAL_FACTOR){
    DSDPEventLogEnd(ConeMaxDStep);
  } else if (flag==PRIMAL_FACTOR){
    DSDPEventLogEnd(ConeMaxPStep);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPPassXVectors(DSDP dsdp, double mu, DSDPVec Y, DSDPVec DY);

\brief Pass the information needed to compute the variables X in each
cone but do not compute X.

\param dsdp the solver
\param mu barrier parameter
\param Y input y variables
\param DY input step direction


\sa DSDPComputeXVariables()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPPassXVectors"
int DSDPPassXVectors(DSDP dsdp, double mu, DSDPVec Y, DSDPVec DY){
  int info,kk;
  DSDPFunctionBegin;
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeSetXMaker(dsdp->K[kk].cone,mu,Y,DY);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPGetConicDimension(DSDP dsdp, double *n);
\brief Get the total dimension of the cones.

\param dsdp the solver
\param n dimension


*/
#undef __FUNCT__
#define __FUNCT__ "DSDPGetConicDimension"
int DSDPGetConicDimension(DSDP dsdp, double *n){
  int info,kk;
  double nn,nnn=0;
  DSDPFunctionBegin;
  for (kk=0;kk<dsdp->ncones;kk++){
    nn=0;
    info=DSDPConeGetDimension(dsdp->K[kk].cone,&nn);DSDPCHKCONEERR(kk,info);
    nnn+=nn;
  }
  *n=nnn;
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPViewCones( DSDP dsdp);

\brief Each cone should print its state.
\param dsdp the solver

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPViewCones"
int DSDPViewCones(DSDP dsdp){
  int info,kk;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeView);
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeView(dsdp->K[kk].cone);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  DSDPEventLogEnd(ConeView);
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPMonitorCones( DSDP dsdp, int tag);

\brief This routine is called once per iteration.
\param dsdp the solver
\param tag allow for multiple monitors

The cone can print statistics, visualize data, terminate solver, or 
whatever it wants.
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPMonitorCones"
int DSDPMonitorCones(DSDP dsdp,int tag){
  int info,kk;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeView);
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    info=DSDPConeMonitor(dsdp->K[kk].cone,tag);DSDPCHKCONEERR(kk,info);
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  DSDPEventLogEnd(ConeView);
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPSchurSparsity( DSDP dsdp, int row, int rnnz[], int m);
\brief Each cone should print its state.
\param dsdp the solver
\param row corresponding to the variable y.
\param rnnz nonzeros indicate a nonzero in the Shur matrix at that column.
\param m size of Schur matrix and the arrow.
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSparsityInSchurMat"
int DSDPSchurSparsity(DSDP dsdp, int row, int rnnz[], int m){
  int info,kk;
  DSDPFunctionBegin;
  for (kk=0;kk<dsdp->ncones;kk++){
    /*    DSDPEventLogBegin(dsdp->K[kk].coneid); */
    info=DSDPConeSparsityInSchurMat(dsdp->K[kk].cone,row,rnnz,m+2);DSDPCHKCONEERR(kk,info);
    /*    DSDPEventLogEnd(dsdp->K[kk].coneid); */
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPComputeLogSDeterminant(DSDP dsdp, double *logdet);
\brief Compute the logarithmic barrier function for the dual varialbe S.
\param dsdp the solver
\param logdet evaluated barrier function 

\sa DSDPComputeSS()
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPComputeLogSDeterminant"
int DSDPComputeLogSDeterminant(DSDP dsdp, double *logdet){
  int info,kk;
  double coneobjective,conepotential,llogdet=0;
  DSDPFunctionBegin;
  DSDPEventLogBegin(ConePotential);
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    coneobjective=0;conepotential=0;
    info=DSDPConeComputeLogSDeterminant(dsdp->K[kk].cone,&coneobjective,&conepotential);DSDPCHKCONEERR(kk,info);
    llogdet+=conepotential;
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  *logdet=llogdet;
  DSDPEventLogEnd(ConePotential);
  DSDPFunctionReturn(0);
}


/*!
\fn int DSDPSetCone( DSDP dsdp, DSDPCone tcone);
\brief Pass a cone to the DSDP solver.
\param dsdp the solver
\param tcone a cone object.

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSetCone"
int DSDPSetCone(DSDP dsdp,DSDPCone tcone){
  int info,i,tc;
  char conename[100];
  DCone *ccones;
  DSDPFunctionBegin;
  if (dsdp->ncones>=dsdp->maxcones){
    tc=2*dsdp->maxcones+4;

    DSDPCALLOC2(&ccones,DCone,tc,&info);DSDPCHKERR(info);
    for (i=0;i<dsdp->ncones;i++){ccones[i].cone=dsdp->K[i].cone;}
    for (i=0;i<dsdp->ncones;i++){ccones[i].coneid=dsdp->K[i].coneid;}
    DSDPFREE(&dsdp->K,&info);DSDPCHKERR(info);
    dsdp->K=ccones;
    dsdp->maxcones=tc;
  }
  info=DSDPGetConeName(tcone,conename,100);DSDPCHKERR(info);
  DSDPEventLogRegister(conename,&tc);
  dsdp->K[dsdp->ncones].cone=tcone;
  dsdp->K[dsdp->ncones].coneid=tc;
  dsdp->ncones++;
  DSDPFunctionReturn(0); 
}

/*!
\fn int DSDPAddCone(DSDP dsdp,struct DSDPCone_Ops* dsdpops, void* dsdpcone);
\param dsdp the solver
\param dsdpops address of a structure with function pointers
\param dsdpcone address of a cone structure
\brief Apply DSDP to a conic structure.

DSDP operates on 
cones such as the semidefinite cone and nonnegative orthant.  Given
variables y from the solver, each cone implements operations such
as computing S, maximum step length, computing the Newton matrix,
and computing the Hessian.  Each operation is well defined by
the dual-scaling algorithm.  A cone that implements these operations
can be added to the DSDP solver.

\sa DSDPCreateSDPCone()
\sa DSDPCreateLPCone()
\sa DSDPCreateBCone()

\todo Add SOCP cone and application-specific cones.

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPAddCone"
int DSDPAddCone(DSDP dsdp,struct DSDPCone_Ops* dsdpops, void* dsdpcone){
  int info;
  DSDPCone K;
  DSDPFunctionBegin;
  info=DSDPConeInitialize(&K); DSDPCHKERR(info);
  info=DSDPConeSetData(&K,dsdpops,dsdpcone); DSDPCHKERR(info);
  info=DSDPSetCone(dsdp,K); DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}


/*!
\fn int DSDPSetSchurMatOps(DSDP dsdp,struct DSDPSchurMat_Ops* sops, void* mdata);

\param dsdp the solver
\param sops address of a structure with function pointers
\param mdata address of a matrix object
\brief Set the Schur complement matrix

The step direction in DSDP is the solution to a set of linear equations.  The
cones used by DSDP compute the elements of the matrix and the right-hand
side vectors.   Any matrix that implements the Schur complement matrix
interface can be used by DSDP.  In addition to factoring a matrix
and solving it, this interface also provides matrix assembly routines
for the cones.

\sa DSDPAddCone()

\todo Use SCALAPACK to assemble, factor, and solve the matrix in parallel.

*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSetSchurMatOps"
int DSDPSetSchurMatOps(DSDP dsdp,struct DSDPSchurMat_Ops* sops, void* mdata){
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatSetData(&dsdp->M,sops,mdata);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}


/*!
\fn int DSDPAddSchurRow( DSDP dsdp, int row, DSDPVec R);

\brief Add a row to the Schur matrix.
\param dsdp the solver
\param row corresponding to which variable y.
\param R the elements of the row.

This routine is called by the conic routines that compute the Hessian matrix.
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSetSchurRow"
int DSDPAddSchurRow(DSDP dsdp,int row, DSDPVec R){
  int info;
  DSDPFunctionBegin;
  info=DSDPSchurMatAddRow(dsdp->M,row,1.0,R);DSDPCHKERR(info);
  DSDPFunctionReturn(0); 
}
/*
#undef __FUNCT__
#define __FUNCT__ "DSDPIdentifySchurColumns"
int DSDPIdentifySchurColumns(DSDP dsdp,int row, int *mcol, int *ncols, int m){
  DSDPFunctionBegin;
  int info;
  DSDPVec V;
  info=DSDPSchurMatRowColumnScaling(dsdp->M,row,V,ncols); DSDPCHKERR(info);
  DSDPFunctionReturn(1); 
}
*/

/*!
\fn int DSDPComputeXVariables(DSDP dsdp, double xmakermu, DSDPVec xmakery, DSDPVec xmakerdy, DSDPVec AX, double *tracexs);

\brief Compute the X variables in each cone.

\param dsdp the solver
\param xmakermu barrier parameter
\param xmakery input y variables
\param xmakerdy input step direction
\param AX output product of X and the data
\param tracexs ouput inner product of X and S.
*/
#undef __FUNCT__  
#define __FUNCT__ "DSDPComputeXVariables"
int DSDPComputeXVariables(DSDP dsdp, double xmakermu, DSDPVec xmakery, DSDPVec xmakerdy, DSDPVec AX, double *tracexs){
  int kk,info;
  double ttracexs=0,tttracexs=0,tracex;

  DSDPFunctionBegin;
  DSDPEventLogBegin(ConeComputeX);
  info=DSDPVecZero(AX);DSDPCHKERR(info);
  for (kk=0;kk<dsdp->ncones;kk++){
    DSDPEventLogBegin(dsdp->K[kk].coneid);
    tttracexs=0;
    info=DSDPConeComputeX(dsdp->K[kk].cone,xmakermu,xmakery,xmakerdy,AX,&tttracexs);DSDPCHKCONEERR(kk,info);
    ttracexs+=tttracexs;
    DSDPEventLogEnd(dsdp->K[kk].coneid);
  }
  info=DSDPVecGetR(AX,&tracex); DSDPCHKERR(info);
  DSDPLogInfo(0,2,"Trace(X): %4.2e\n",dsdp->tracex);
  info=DSDPVecAXPY(-1.0,dsdp->b,AX); DSDPCHKERR(info);
  info=DSDPComputeFixedYX(dsdp->M,AX); DSDPCHKERR(info);
  *tracexs=ttracexs;
  info=DSDPVecSetR(AX,tracex); DSDPCHKERR(info);
  DSDPEventLogEnd(ConeComputeX);
  DSDPFunctionReturn(0);
}


