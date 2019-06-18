/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
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
 * This is the public header file for the IDAS scaled preconditioned
 * TFQMR linear solver module, IDASPTFQMR.
 *
 * Part I contains function prototypes for using IDASPTFQMR on forward 
 * problems (DAE integration and/or FSA)
 *
 * Part II contains function prototypes for using IDASPTFQMR on adjoint 
 * (backward) problems
 * -----------------------------------------------------------------
 */

#ifndef _IDASSPTFQMR_H
#define _IDASSPTFQMR_H

#include <idas/idas_spils.h>
#include <sundials/sundials_sptfqmr.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* 
 * -----------------------------------------------------------------
 * PART I - forward problems
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : IDASptfqmr
 * -----------------------------------------------------------------
 * A call to the IDASptfqmr function links the main integrator with
 * the IDASPTFQMR linear solver module. Its parameters are as
 * follows:
 *
 * IDA_mem  is the pointer to memory block returned by IDACreate.
 *
 * maxl     is the maximum Krylov subspace dimension, an
 *          optional input. Pass 0 to use the default value.
 *          Otherwise pass a positive integer.
 *
 * The return values of IDASptfqmr are:
 *    IDASPILS_SUCCESS    if successful
 *    IDASPILS_MEM_NULL   if the IDAS memory was NULL
 *    IDASPILS_MEM_FAIL   if there was a memory allocation failure
 *    IDASPILS_ILL_INPUT  if there was illegal input.
 * The above constants are defined in idas_spils.h
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASptfqmr(void *ida_mem, int maxl);

/* 
 * -----------------------------------------------------------------
 * PART II - backward problems
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASptfqmrB(void *ida_mem, int which, int maxlB);

#ifdef __cplusplus
}
#endif

#endif
