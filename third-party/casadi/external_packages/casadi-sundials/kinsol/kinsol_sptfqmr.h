/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
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
 * This is the public header file for the KINSOL scaled preconditioned
 * TFQMR linear solver module, KINSPTFQMR.
 * -----------------------------------------------------------------
 */

#ifndef _KINSPTFQMR_H
#define _KINSPTFQMR_H

#include <kinsol/kinsol_spils.h>
#include <sundials/sundials_sptfqmr.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmr
 * -----------------------------------------------------------------
 * KINSptfqmr links the main KINSOL solver module with the SPTFQMR
 * linear solver module. The routine establishes the inter-module
 * interface by setting the generic KINSOL pointers linit, lsetup,
 * lsolve, and lfree to KINSptfqmrInit, KINSptfqmrSetup, KINSptfqmrSolve,
 * and KINSptfqmrFree, respectively.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  maxl  maximum allowable dimension of Krylov subspace (passing
 *        a value of 0 (zero) will cause the default value
 *        KINSPTFQMR_MAXL (predefined constant) to be used)
 *
 * If successful, KINSptfqmr returns KINSPTFQMR_SUCCESS. If an error
 * occurs, then KINSptfqmr returns an error code (negative integer
 * value).
 *
 * -----------------------------------------------------------------
 * KINSptfqmr Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSptfqmr subroutine are the
 * following:
 *
 * KINSPTFQMR_SUCCESS : means the KINSPTFQMR linear solver module
 *                      (implementation of the TFQMR method) was
 *                      successfully initialized - allocated system
 *                      memory and set shared variables to default
 *                      values [0]
 *
 * KINSPTFQMR_MEM_NULL : means a NULL KINSOL memory block pointer
 *                       was given (must call the KINCreate and
 *                       KINMalloc memory allocation subroutines
 *                       prior to calling KINSptfqmr) [-1]
 *
 * KINSPTFQMR_MEM_FAIL : means either insufficient system resources
 *                       were available to allocate memory for the
 *                       main KINSPTFQMR data structure (type
 *                       KINSptfqmrMemRec), or the SptfqmrMalloc
 *                       subroutine failed (unable to allocate enough
 *                       system memory for vector storate and/or the
 *                       main SPTFQMR data structure
 *                       (type SptfqmrMemRec)) [-4]
 *
 * KINSPTFQMR_ILL_INPUT : means either a supplied parameter was invalid,
 *                        or the NVECTOR implementation is NOT
 *                        compatible [-3]
 *
 * The above constants are defined in kinsol_spils.h
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSptfqmr(void *kinmem, int maxl);


#ifdef __cplusplus
}
#endif

#endif
