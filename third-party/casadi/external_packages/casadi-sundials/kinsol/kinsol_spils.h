/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott Cohen, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
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
 * This is the common header file for the Scaled Preconditioned
 * Iterative Linear Solvers in KINSOL.
 * -----------------------------------------------------------------
 */

#ifndef _KINSPILS_H
#define _KINSPILS_H

#include <sundials/sundials_iterative.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * KINSPILS return values
 * -----------------------------------------------------------------
 */

#define KINSPILS_SUCCESS    0

#define KINSPILS_MEM_NULL  -1
#define KINSPILS_LMEM_NULL -2
#define KINSPILS_ILL_INPUT -3
#define KINSPILS_MEM_FAIL  -4
#define KINSPILS_PMEM_NULL -5

/*
 * -----------------------------------------------------------------
 * KINSPILS solver constant
 * -----------------------------------------------------------------
 * KINSPILS_MAXL : default maximum dimension of Krylov subspace
 * -----------------------------------------------------------------
 */

#define KINSPILS_MAXL 10

/*
 * -----------------------------------------------------------------
 * Type : KINSpilsPrecSetupFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner setup subroutine should
 * compute the right-preconditioner matrix P (stored in memory
 * block referenced by P_data pointer) used to form the
 * scaled preconditioned linear system:
 *
 *  (Df*J(uu)*(P^-1)*(Du^-1)) * (Du*P*x) = Df*(-F(uu))
 *
 * where Du and Df denote the diagonal scaling matrices whose
 * diagonal elements are stored in the vectors uscale and
 * fscale, repsectively.
 *
 * The preconditioner setup routine (referenced by iterative linear
 * solver modules via pset (type KINSpilsPrecSetupFn)) will not be
 * called prior to every call made to the psolve function, but will
 * instead be called only as often as necessary to achieve convergence
 * of the Newton iteration.
 *
 * Note: If the psolve routine requires no preparation, then a
 * preconditioner setup function need not be given.
 *
 *  uu  current iterate (unscaled) [input]
 *
 *  uscale  vector (type N_Vector) containing diagonal elements
 *          of scaling matrix for vector uu [input]
 *
 *  fval  vector (type N_Vector) containing result of nonliear
 *        system function evaluated at current iterate:
 *        fval = F(uu) [input]
 *
 *  fscale  vector (type N_Vector) containing diagonal elements
 *          of scaling matrix for fval [input]
 *
 *  user_data  pointer to user-allocated data memory block
 *
 *  vtemp1/vtemp2  available scratch vectors (temporary storage)
 *
 * If successful, the function should return 0 (zero). If an error
 * occurs, then the routine should return a non-zero integer value.
 * -----------------------------------------------------------------
 */

typedef int (*KINSpilsPrecSetupFn)(N_Vector uu, N_Vector uscale,
                                   N_Vector fval, N_Vector fscale,
                                   void *user_data, N_Vector vtemp1,
				   N_Vector vtemp2);

/*
 * -----------------------------------------------------------------
 * Type : KINSpilsPrecSolveFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner solve subroutine (referenced
 * by iterative linear solver modules via psolve (type
 * KINSpilsPrecSolveFn)) should solve a (scaled) preconditioned
 * linear system of the generic form P*z = r, where P denotes the
 * right-preconditioner matrix computed by the pset routine.
 *
 *  uu  current iterate (unscaled) [input]
 *
 *  uscale  vector (type N_Vector) containing diagonal elements
 *          of scaling matrix for vector uu [input]
 *
 *  fval  vector (type N_Vector) containing result of nonliear
 *        system function evaluated at current iterate:
 *        fval = F(uu) [input]
 *
 *  fscale  vector (type N_Vector) containing diagonal elements
 *          of scaling matrix for fval [input]
 *
 *  vv  vector initially set to the right-hand side vector r, but
 *      which upon return contains a solution of the linear system
 *      P*z = r [input/output]
 *
 *  user_data  pointer to user-allocated data memory block
 *
 *  vtemp  available scratch vector (volatile storage)
 *
 * If successful, the function should return 0 (zero). If a
 * recoverable error occurs, then the subroutine should return
 * a positive integer value (in this case, KINSOL attempts to
 * correct by calling the preconditioner setup function if the 
 * preconditioner information is out of date). If an unrecoverable 
 * error occurs, then the preconditioner solve function should return 
 * a negative integer value.
 * -----------------------------------------------------------------
 */

typedef int (*KINSpilsPrecSolveFn)(N_Vector uu, N_Vector uscale, 
                                   N_Vector fval, N_Vector fscale, 
                                   N_Vector vv, void *user_data,
                                   N_Vector vtemp);

/*
 * -----------------------------------------------------------------
 * Type : KINSpilsJacTimesVecFn
 * -----------------------------------------------------------------
 * The (optional) user-supplied matrix-vector product subroutine
 * (referenced internally via jtimes (type KINSpilsJacTimesVecFn))
 * is used to compute Jv = J(uu)*v (system Jacobian applied to a
 * given vector). If a user-defined routine is not given, then the
 * private routine is used.
 *
 *  v  unscaled variant of vector to be multiplied by J(uu) [input]
 *
 *  Jv  vector containing result of matrix-vector product J(uu)*v
 *      [output]
 *
 *  uu  current iterate (unscaled) [input]
 *
 *  new_uu  flag (reset by user) indicating if the iterate uu
 *          has been updated in the interim - Jacobian needs
 *          to be updated/reevaluated, if appropriate, unless
 *          new_uu = FALSE [input/output]
 *
 *  user_data  pointer to user data, the same as the user_data
 *          parameter passed to the KINSetUserData function.
 *
 * If successful, the function should return 0 (zero). If an error
 * occurs, then the routine should return a non-zero integer value.
 * -----------------------------------------------------------------
 */

typedef int (*KINSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv,
                                     N_Vector uu, booleantype *new_uu, 
                                     void *J_data);





/*
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs:
 *
 *       Function Name       |   Optional Input  [Default Value]
 *                           |
 * -----------------------------------------------------------------
 *                           |
 * KINSpilsSetMaxRestarts    | maximum number of times the SPGMR
 *                           | (scaled preconditioned GMRES) linear
 *                           | solver can be restarted
 *                           | [0]
 *                           |
 * KINSpilsSetPreconditioner | used to set the following:
 *                           |   (a) name of user-supplied routine
 *                           |       used to compute a preconditioner
 *                           |       matrix for the given linear
 *                           |       system (pset)
 *                           |       [NULL]
 *                           |   (b) name of user-supplied routine
 *                           |       used to apply preconditioner to
 *                           |       linear system (psolve)
 *                           |       [NULL]
 *                           |
 * KINSpilsSetJacTimesVecFn  | used to set the following the name
 *                           | of user-supplied subroutine used to 
 *                           | compute the matrix-vector product J(u)*v,
 *                           | where J denotes the system Jacobian.
 *                           | [KINSpilsDQJtimes]
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSpilsSetMaxRestarts(void *kinmem, int maxrs);
SUNDIALS_EXPORT int KINSpilsSetPreconditioner(void *kinmem,
					      KINSpilsPrecSetupFn pset,
					      KINSpilsPrecSolveFn psolve);
SUNDIALS_EXPORT int KINSpilsSetJacTimesVecFn(void *kinmem,
                                             KINSpilsJacTimesVecFn jtv);

/*
 * -----------------------------------------------------------------
 * KINSpilsSet* Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSpilsSet* subroutines
 * are the following:
 *
 * KINSPILS_SUCCESS : means the associated parameter was successfully
 *                    set [0]
 *
 * KINSPILS_ILL_INPUT : means the supplied parameter was invalid
 *                      (check error message) [-3]
 *
 * KINSPILS_MEM_NULL : means a NULL KINSOL memory block pointer
 *                     was given [-1]
 *
 * KINSPILS_LMEM_NULL : means system memory has not yet been
 *                      allocated for the linear solver 
 *                      (lmem == NULL) [-2]
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Optional Output Extraction Functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistical information related to the KINSPILS linear
 * solvers:
 *
 *        Function Name       |      Returned Value
 *                            |
 * -----------------------------------------------------------------
 *                            |
 * KINSpilsGetWorkSpace       | returns both integer workspace size
 *                            | (total number of long int-sized blocks
 *                            | of memory allocated  for
 *                            | vector storage), and real workspace
 *                            | size (total number of realtype-sized
 *                            | blocks of memory allocated
 *                            | for vector storage)
 *                            |
 * KINSpilsGetNumPrecEvals    | total number of preconditioner
 *                            | evaluations (number of calls made
 *                            | to the user-defined pset routine)
 *                            |
 * KINSpilsGetNumPrecSolves   | total number of times preconditioner
 *                            | was applied to linear system (number
 *                            | of calls made to the user-supplied
 *                            | psolve function)
 *                            |
 * KINSpilsGetNumLinIters     | total number of linear iterations
 *                            | performed
 *                            |
 * KINSpilsGetNumConvFails    | total number of linear convergence
 *                            | failures
 *                            |
 * KINSpilsGetNumJtimesEvals  | total number of times the matrix-
 *                            | vector product J(u)*v was computed
 *                            | (number of calls made to the jtimes
 *                            | subroutine)
 *                            |
 * KINSpilsGetNumFuncEvals    | total number of evaluations of the
 *                            | system function F(u) (number of
 *                            | calls made to the user-supplied
 *                            | func routine by the linear solver
 *                            | module member subroutines)
 *                            |
 * KINSpilsGetLastFlag        | returns the last flag returned by
 *                            | the linear solver
 *                            |
 * KINSpilsGetReturnFlagName  | returns the name of the constant
 *                            | associated with a KINSPILS return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSpilsGetWorkSpace(void *kinmem, long int *lenrwSG, long int *leniwSG);
SUNDIALS_EXPORT int KINSpilsGetNumPrecEvals(void *kinmem, long int *npevals);
SUNDIALS_EXPORT int KINSpilsGetNumPrecSolves(void *kinmem, long int *npsolves);
SUNDIALS_EXPORT int KINSpilsGetNumLinIters(void *kinmem, long int *nliters);
SUNDIALS_EXPORT int KINSpilsGetNumConvFails(void *kinmem, long int *nlcfails);
SUNDIALS_EXPORT int KINSpilsGetNumJtimesEvals(void *kinmem, long int *njvevals);
SUNDIALS_EXPORT int KINSpilsGetNumFuncEvals(void *kinmem, long int *nfevalsS); 
SUNDIALS_EXPORT int KINSpilsGetLastFlag(void *kinmem, long int *flag);
SUNDIALS_EXPORT char *KINSpilsGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
