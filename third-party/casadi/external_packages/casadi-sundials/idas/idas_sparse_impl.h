/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
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
 * Implementation header file for the IDASLS linear solvers.
 * -----------------------------------------------------------------
 */

#ifndef _IDASSPARSE_IMPL_H
#define _IDASSPARSE_IMPL_H

#include "idas/idas_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 * I D A S S P A R S E    I N T E R N A L    C O N S T A N T S
 * =================================================================
 */

/*
 * =================================================================
 * PART I:  F O R W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types : IDASlsMemRec, IDASlsMem                             
 * -----------------------------------------------------------------
 * IDASlsMem is pointer to a IDASlsMemRec structure.
 * -----------------------------------------------------------------
 */

typedef struct IDASlsMemRec {

  IDASlsSparseJacFn s_jaceval; /* user Jacobian evaluation routine 
				  to be called */
  void *s_jacdata;           /* J_data is passed to djac or bjac */

  long int s_nje;           /* nje = no. of calls to jac */

  long int s_last_flag;     /* last error return flag */

  int s_first_factorize;    /* flag telling whether the first 
			       factorization needs to happen */
  SlsMat s_JacMat;          /* J = dF/dy + cj*dF/dy' */

  void *s_solver_data;      /* structure for solver-specific data */
  

} *IDASlsMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */
  

/*
 * =================================================================
 * PART II:  B A C K W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types : IDASlsMemRecB, IDASlsMemB       
 * -----------------------------------------------------------------
 * An IDASLS linear solver's specification function attaches such
 * a structure to the lmemB filed of IDABMem
 * -----------------------------------------------------------------
 */

typedef struct IDASlsMemRecB {

  IDASlsSparseJacFnB s_djacB;
  IDASlsSparseJacFnBS s_djacBS;

} *IDASlsMemB;


/*
 * =================================================================
 * E R R O R   M E S S A G E S
 * =================================================================
 */

#define MSGSP_IDAMEM_NULL "Integrator memory is NULL."
#define MSGSP_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGSP_MEM_FAIL "A memory request failed."
#define MSGSP_LMEM_NULL "Linear solver memory is NULL."
#define MSGSP_ILL_INPUT "Invalid input detected."
#define MSGSP_JAC_NOSET "Jacobian evaluation function has not been set."
#define MSGSP_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."
#define MSGSP_PACKAGE_FAIL "A call to an external package failed."

#define MSGSP_CAMEM_NULL "idaadj_mem = NULL illegal."
#define MSGSP_LMEMB_NULL "Linear solver memory is NULL for the backward integration."
#define MSGSP_BAD_T "Bad t for interpolation."
#define MSGSP_BAD_WHICH "Illegal value for which."
#define MSGSP_NO_ADJ "Illegal attempt to call before calling IDAAdjInit."

#ifdef __cplusplus
}
#endif

#endif
