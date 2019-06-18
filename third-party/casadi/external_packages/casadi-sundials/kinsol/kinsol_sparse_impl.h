/*
 * -----------------------------------------------------------------
 * $Revision: 4357 $
 * $Date: 2015-02-09 13:22:31 -0800 (Mon, 09 Feb 2015) $
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
 * Implementation header file for the KINSLS linear solvers.
 * -----------------------------------------------------------------
 */

#ifndef _KINSPARSE_IMPL_H
#define _KINSPARSE_IMPL_H

#include "kinsol/kinsol_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 * K I N S P A R S E    I N T E R N A L    C O N S T A N T S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types : KINSlsMemRec, KINSlsMem                             
 * -----------------------------------------------------------------
 * KINSlsMem is pointer to a KINSlsMemRec structure.
 * -----------------------------------------------------------------
 */

typedef struct KINSlsMemRec {

  KINSlsSparseJacFn s_jaceval; /* user Jacobian evaluation routine 
				  to be called */
  void *s_jacdata;           /* J_data is passed to eval. routine */

  long int s_nje;           /* nje = no. of calls to jac */

  long int s_last_flag;     /* last error return flag */

  int s_first_factorize;    /* flag telling whether the first 
			       factorization needs to happen */
  SlsMat s_JacMat;          /* J = dF/du */

  void *s_solver_data;      /* structure for solver-specific data */
  

} *KINSlsMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */
  

/*
 * =================================================================
 * E R R O R   M E S S A G E S
 * =================================================================
 */

#define MSGSP_KINMEM_NULL "Solver memory is NULL."
#define MSGSP_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGSP_MEM_FAIL "A memory request failed."
#define MSGSP_LMEM_NULL "Linear solver memory is NULL."
#define MSGSP_JAC_NOSET "Jacobian evaluation function has not been set."
#define MSGSP_ILL_INPUT "Invalid input detected."
#define MSGSP_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."
#define MSGSP_PACKAGE_FAIL "A call to an external package failed."

#ifdef __cplusplus
}
#endif

#endif
