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
 * This is the header file for the CVSBANDPRE module, which
 * provides a banded difference quotient Jacobian-based
 * preconditioner and solver routines for use with CVSPGMR,
 * CVSPBCG, or CVSPTFQMR.
 *
 * Part I contains type definitions and function prototypes for using
 * CVSBANDPRE on forward problems (IVP integration and/or FSA)
 *
 * Part II contains type definitions and function prototypes for using
 * CVSBANDPRE on adjopint (backward) problems
 * -----------------------------------------------------------------
 */

#ifndef _CVSBANDPRE_H
#define _CVSBANDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* 
 * =================================================================
 * PART I - forward problems
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 *
 * SUMMARY
 * 
 * These routines provide a band matrix preconditioner based on
 * difference quotients of the ODE right-hand side function f.
 * The user supplies parameters
 *   mu = upper half-bandwidth (number of super-diagonals)
 *   ml = lower half-bandwidth (number of sub-diagonals)
 * The routines generate a band matrix of bandwidth ml + mu + 1
 * and use this to form a preconditioner for use with the Krylov
 * linear solver in CVSP*. Although this matrix is intended to
 * approximate the Jacobian df/dy, it may be a very crude
 * approximation. The true Jacobian need not be banded, or its
 * true bandwith may be larger than ml + mu + 1, as long as the
 * banded approximation generated here is sufficiently accurate
 * to speed convergence as a preconditioner.
 *
 * Usage:
 *   The following is a summary of the usage of this module.
 *   Details of the calls to CVodeCreate, CVodeInit, CVSp*,
 *   and CVode are available in the User Guide.
 *   To use these routines, the sequence of calls in the user
 *   main program should be as follows:
 *
 *   #include <cvodes/cvodes_bandpre.h>
 *   #include <nvector_serial.h>
 *   ...
 *   Set y0
 *   ...
 *   cvode_mem = CVodeCreate(...);
 *   ier = CVodeInit(...);
 *   ...
 *   flag = CVSptfqmr(cvode_mem, pretype, maxl);
 *     -or-
 *   flag = CVSpgmr(cvode_mem, pretype, maxl);
 *     -or-
 *   flag = CVSpbcg(cvode_mem, pretype, maxl);
 *   ...
 *   flag = CVBandPrecInit(cvode_mem, N, mu, ml);
 *   ...
 *   flag = CVode(...);
 *   ...
 *   Free y0
 *   ...
 *   CVodeFree(&cvode_mem);
 *
 * Notes:
 * (1) Include this file for the CVBandPrecData type definition.
 * (2) In the CVBandPrecInit call, the arguments N is the
 *     problem dimension.
 * (3) In the CVBPSp* call, the user is free to specify
 *     the input pretype and the optional input maxl.
 * -----------------------------------------------------------------
 */


/*
 * -----------------------------------------------------------------
 * Function : CVBandPrecInit
 * -----------------------------------------------------------------
 * CVBandPrecInit allocates and initializes the BANDPRE preconditioner
 * module. This functino must be called AFTER one of the SPILS linear
 * solver modules has been attached to the CVODE integrator.
 *
 * The parameters of CVBandPrecInit are as follows:
 *
 * cvode_mem is the pointer to CVODE memory returned by CVodeCreate.
 *
 * N is the problem size.
 *
 * mu is the upper half bandwidth.
 *
 * ml is the lower half bandwidth.
 *
 * The return value of CVBandPrecInit is one of:
 *   CVSPILS_SUCCESS if no errors occurred
 *   CVSPILS_MEM_NULL if the integrator memory is NULL
 *   CVSPILS_LMEM_NULL if the linear solver memory is NULL
 *   CVSPILS_ILL_INPUT if an input has an illegal value
 *   CVSPILS_MEM_FAIL if a memory allocation request failed
 *
 * NOTE: The band preconditioner assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVBandPrecInit will
 *       first test for a compatible N_Vector internal
 *       representation by checking for required functions.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVBandPrecInit(void *cvode_mem, long int N, long int mu, long int ml);

/*
 * -----------------------------------------------------------------
 * Optional output functions : CVBandPrecGet*
 * -----------------------------------------------------------------
 * CVBandPrecGetWorkSpace returns the real and integer work space used
 *                        by CVBANDPRE.
 * CVBandPrecGetNumRhsEvals returns the number of calls made from
 *                          CVBANDPRE to the user's right-hand side
 *                          routine f.
 *
 * The return value of CVBandPrecGet* is one of:
 *   CVSPILS_SUCCESS if no errors occurred
 *   CVSPILS_MEM_NULL if the integrator memory is NULL
 *   CVSPILS_LMEM_NULL if the linear solver memory is NULL
 *   CVSPILS_PMEM_NULL if the preconditioner memory is NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVBandPrecGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CVBandPrecGetNumRhsEvals(void *cvode_mem, long int *nfevalsBP);

/* 
 * =================================================================
 * PART II - backward problems
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Functions: CVBandPrecInitB, CVBPSp*B
 * -----------------------------------------------------------------
 * Interface functions for the CVBANDPRE preconditioner to be used
 * on the backward phase.
 *
 * CVBandPrecInitB interfaces to the CVBANDPRE preconditioner
 * for the backward integration.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVBandPrecInitB(void *cvode_mem, int which,
                                    long int nB, long int muB, long int mlB);

#ifdef __cplusplus
}
#endif

#endif
