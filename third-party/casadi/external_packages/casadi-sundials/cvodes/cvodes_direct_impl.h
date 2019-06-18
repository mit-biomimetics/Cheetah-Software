/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Common implementation header file for the CVDLS linear solvers.
 * -----------------------------------------------------------------
 */

#ifndef _CVSDLS_IMPL_H
#define _CVSDLS_IMPL_H

#include <cvodes/cvodes_direct.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 * C V S D I R E C T    I N T E R N A L    C O N S T A N T S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * CVDLS solver constants
 * -----------------------------------------------------------------
 * CVD_MSBJ   maximum number of steps between Jacobian evaluations
 * CVD_DGMAX  maximum change in gamma between Jacobian evaluations
 * -----------------------------------------------------------------
 */

#define CVD_MSBJ  50
#define CVD_DGMAX RCONST(0.2)

/*
 * =================================================================
 * PART I:  F O R W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types: CVDlsMemRec, CVDlsMem                             
 * -----------------------------------------------------------------
 * CVDlsMem is pointer to a CVDlsMemRec structure.
 * -----------------------------------------------------------------
 */

typedef struct CVDlsMemRec {

  int d_type;             /* SUNDIALS_DENSE or SUNDIALS_BAND              */

  long int d_n;           /* problem dimension                            */

  long int d_ml;          /* lower bandwidth of Jacobian                  */
  long int d_mu;          /* upper bandwidth of Jacobian                  */ 
  long int d_smu;         /* upper bandwith of M = MIN(N-1,d_mu+d_ml)     */

  booleantype d_jacDQ;    /* TRUE if using internal DQ Jacobian approx.   */
  CVDlsDenseJacFn d_djac; /* dense Jacobian routine to be called          */
  CVDlsBandJacFn d_bjac;  /* band Jacobian routine to be called           */
  void *d_J_data;         /* data pointer passed to djac or bjac          */

  DlsMat d_M;             /* M = I - gamma * df/dy                        */
  DlsMat d_savedJ;        /* savedJ = old Jacobian                        */

  int *d_pivots;          /* pivots = int pivot array for PM = LU         */
  long int *d_lpivots;    /* lpivots = long int pivot array for PM = LU   */
  
  long int  d_nstlj;      /* nstlj = nst at last Jacobian eval.           */

  long int d_nje;         /* nje = no. of calls to jac                    */

  long int d_nfeDQ;       /* no. of calls to f due to DQ Jacobian approx. */

  long int d_last_flag;   /* last error return flag                       */
  
} *CVDlsMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

int cvDlsDenseDQJac(long int N, realtype t,
		    N_Vector y, N_Vector fy, 
		    DlsMat Jac, void *data,
		    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
int cvDlsBandDQJac(long int N, long int mupper, long int mlower,
		   realtype t, N_Vector y, N_Vector fy, 
		   DlsMat Jac, void *data,
		   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


/*
 * =================================================================
 * PART II:  B A C K W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types : CVDlsMemRecB, CVDlsMemB       
 * -----------------------------------------------------------------
 * A CVDLS linear solver's specification function attaches such
 * a structure to the lmemB filed of CVodeBMem
 * -----------------------------------------------------------------
 */

typedef struct CVDlsMemRecB {

  int d_typeB;

  CVDlsDenseJacFnB d_djacB;
  CVDlsDenseJacFnBS d_djacBS;
  CVDlsBandJacFnB d_bjacB;
  CVDlsBandJacFnBS d_bjacBS;

} *CVDlsMemB;


/*
 * =================================================================
 * E R R O R   M E S S A G E S
 * =================================================================
 */

#define MSGD_CVMEM_NULL "Integrator memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_MEM_FAIL "A memory request failed."
#define MSGD_LMEM_NULL "Linear solver memory is NULL."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."


#define MSGD_NO_ADJ "Illegal attempt to call before calling CVodeAdjMalloc."
#define MSGD_BAD_WHICH "Illegal value for which."
#define MSGD_LMEMB_NULL "Linear solver memory is NULL for the backward integration."
#define MSGD_BAD_TINTERP "Bad t for interpolation."

#ifdef __cplusplus
}
#endif

#endif
