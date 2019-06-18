/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2013, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
 This is the implementation file for the preconditioned conjugate 
 gradient solver in SUNDIALS.
 --------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_pcg.h>
#include <sundials/sundials_math.h>


/*---------------------------------------------------------------
 private constants
 --------------------------------------------------------------*/
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)


/*---------------------------------------------------------------
 Function : PcgMalloc
 --------------------------------------------------------------*/
PcgMem PcgMalloc(int l_max, N_Vector vec_tmpl)
{
  PcgMem mem;
  N_Vector r, p, z, Ap;

  /* Check the input parameters */
  if (l_max <= 0) return(NULL);

  /* Create temporary arrays */
  r = N_VClone(vec_tmpl);
  if (r == NULL) {
    return(NULL);
  }

  p = N_VClone(vec_tmpl);
  if (p == NULL) {
    N_VDestroy(r);
    return(NULL);
  }

  z = N_VClone(vec_tmpl);
  if (z == NULL) {
    N_VDestroy(r);
    N_VDestroy(p);
    return(NULL);
  }

  Ap = N_VClone(vec_tmpl);
  if (Ap == NULL) {
    N_VDestroy(r);
    N_VDestroy(p);
    N_VDestroy(z);
    return(NULL);
  }

  /* Get memory for an PcgMemRec containing PCG vectors */
  mem = NULL;
  mem = (PcgMem) malloc(sizeof(PcgMemRec));
  if (mem == NULL) {
    N_VDestroy(r);
    N_VDestroy(p);
    N_VDestroy(z);
    N_VDestroy(Ap);
    return(NULL);
  }

  /* Set the structure fields */
  mem->l_max = l_max;
  mem->r     = r;
  mem->p     = p;
  mem->z     = z;
  mem->Ap    = Ap;

  /* Return the pointer to PCG memory */
  return(mem);
}


/*---------------------------------------------------------------
 Function : PcgSolve
 --------------------------------------------------------------*/
int PcgSolve(PcgMem mem, void *A_data, N_Vector x, N_Vector b,
	     int pretype, realtype delta, void *P_data,
	     N_Vector w, ATimesFn atimes, PSolveFn psolve,
	     realtype *res_norm, int *nli, int *nps)
{
  realtype alpha, beta, r0_norm, rho, rz, rz_old;
  N_Vector r, p, z, Ap;
  booleantype UsePrec, converged;
  int l, l_max, ier;

  if (mem == NULL)  return(PCG_MEM_NULL);

  /* Make local copies of mem variables */
  l_max  = mem->l_max;
  r      = mem->r;
  p      = mem->p;
  z      = mem->z;
  Ap     = mem->Ap;

  /* Initialize counters and converged flag */
  *nli = *nps = 0;
  converged = FALSE;

  /* Set preconditioning flag */
  UsePrec = ((pretype == PREC_BOTH) || (pretype == PREC_LEFT) || (pretype == PREC_RIGHT));

  /* Set r to initial residual r_0 = b - A*x_0 */
  if (N_VDotProd(x, x) == ZERO)  N_VScale(ONE, b, r);
  else {
    ier = atimes(A_data, x, r);
    if (ier != 0)
      return((ier < 0) ? PCG_ATIMES_FAIL_UNREC : PCG_ATIMES_FAIL_REC);
    N_VLinearSum(ONE, b, -ONE, r, r);
  }

  /* Set rho to L2 norm of r, and return if small */
  *res_norm = r0_norm = rho = N_VWrmsNorm(r,w);
  if (rho <= delta) return(PCG_SUCCESS);

  /* Apply preconditioner and b-scaling to r = r_0 */
  if (UsePrec) {
    ier = psolve(P_data, r, z, PREC_LEFT);   /* z = P^{-1}r */
    (*nps)++;
    if (ier != 0) return((ier < 0) ? PCG_PSOLVE_FAIL_UNREC : PCG_PSOLVE_FAIL_REC);
  }
  else N_VScale(ONE, r, z);

  /* Initialize rz to <r,z> */
  rz = N_VDotProd(r, z);

  /* Copy z to p */
  N_VScale(ONE, z, p);

  /* Begin main iteration loop */
  for(l=0; l<l_max; l++) {

    /* increment counter */
    (*nli)++;

    /* Generate Ap = A*p */
    ier = atimes(A_data, p, Ap );
    if (ier != 0)
      return((ier < 0) ? PCG_ATIMES_FAIL_UNREC : PCG_ATIMES_FAIL_REC);

    /* Calculate alpha = <r,z> / <Ap,p> */
    alpha = rz / N_VDotProd(Ap, p);

    /* Update x = x + alpha*p */
    N_VLinearSum(ONE, x, alpha, p, x);

    /* Update r = r - alpha*Ap */
    N_VLinearSum(ONE, r, -alpha, Ap, r);

    /* Set rho and check convergence */
    *res_norm = rho = N_VWrmsNorm(r, w);
    if (rho <= delta) {
      converged = TRUE;
      break;
    }

    /* Apply preconditioner:  z = P^{-1}*r */
    if (UsePrec) {
      ier = psolve(P_data, r, z, PREC_LEFT);
      (*nps)++;
      if (ier != 0) return((ier < 0) ? PCG_PSOLVE_FAIL_UNREC : PCG_PSOLVE_FAIL_REC);
    }
    else N_VScale(ONE, r, z);

    /* update rz */
    rz_old = rz;
    rz = N_VDotProd(r, z);
    
    /* Calculate beta = <r,z> / <r_old,z_old> */
    beta = rz / rz_old;

    /* Update p = z + beta*p */
    N_VLinearSum(ONE, z, beta, p, p);

  }

  /* Main loop finished, return with result */
  if (converged == TRUE)  return(PCG_SUCCESS);
  if (rho < r0_norm)      return(PCG_RES_REDUCED);
  return(PCG_CONV_FAIL);
}


/*---------------------------------------------------------------
 Function : PcgFree
 --------------------------------------------------------------*/
void PcgFree(PcgMem mem)
{
  if (mem == NULL) return;

  N_VDestroy(mem->r);
  N_VDestroy(mem->p);
  N_VDestroy(mem->z);
  N_VDestroy(mem->Ap);

  free(mem); mem = NULL;
}


/*===============================================================
   EOF
===============================================================*/
