/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
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
 * This is the header file for the KINSOL dense linear solver module, 
 * KINDENSE.
 * -----------------------------------------------------------------
 */

#include <kinsol/kinsol_direct.h>
#include <sundials/sundials_dense.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _KINDENSE_H
#define _KINDENSE_H

/*
 * -----------------------------------------------------------------
 * Function : KINDense
 * -----------------------------------------------------------------
 * A call to the KINDense function links the main solver with the
 * KINDENSE linear solver. Its arguments are as follows:
 *
 * kinmem - pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 * N      - problem size
 *
 * The return value of KINDense is one of:
 *    KINDLS_SUCCESS   if successful
 *    KINDLS_MEM_NULL  if the kinsol memory was NULL
 *    KINDLS_MEM_FAIL  if there was a memory allocation failure
 *    KINDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINDense(void *kinmem, long int N);

#endif

#ifdef __cplusplus
}
#endif
