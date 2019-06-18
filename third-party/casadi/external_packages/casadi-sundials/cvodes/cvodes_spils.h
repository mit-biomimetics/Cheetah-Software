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
 * This is the common header file for the Scaled, Preconditioned
 * Iterative Linear Solvers in CVODES.
 *
 * Part I contains type definitions and functions for using the 
 * iterative linear solvers on forward problems 
 * (IVP integration and/or FSA)
 *
 * Part II contains type definitions and functions for using the 
 * iterative linear solvers on adjoint (backward) problems
 * -----------------------------------------------------------------
 */

#ifndef _CVSSPILS_H
#define _CVSSPILS_H

#include <sundials/sundials_iterative.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * CVSPILS return values
 * -----------------------------------------------------------------
 */

#define CVSPILS_SUCCESS          0
#define CVSPILS_MEM_NULL        -1
#define CVSPILS_LMEM_NULL       -2
#define CVSPILS_ILL_INPUT       -3
#define CVSPILS_MEM_FAIL        -4
#define CVSPILS_PMEM_NULL       -5

/* Return values for the adjoint module */

#define CVSPILS_NO_ADJ          -101
#define CVSPILS_LMEMB_NULL      -102

/*
 * -----------------------------------------------------------------
 * CVSPILS solver constants
 * -----------------------------------------------------------------
 * CVSPILS_MAXL   : default value for the maximum Krylov
 *                  dimension
 *
 * CVSPILS_MSBPRE : maximum number of steps between
 *                  preconditioner evaluations
 *
 * CVSPILS_DGMAX  : maximum change in gamma between
 *                  preconditioner evaluations
 *
 * CVSPILS_EPLIN  : default value for factor by which the
 *                  tolerance on the nonlinear iteration is
 *                  multiplied to get a tolerance on the linear
 *                  iteration
 * -----------------------------------------------------------------
 */

#define CVSPILS_MAXL   5
#define CVSPILS_MSBPRE 50
#define CVSPILS_DGMAX  RCONST(0.2)
#define CVSPILS_EPLIN  RCONST(0.05)

/* 
 * -----------------------------------------------------------------
 * PART I - forward problems
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Type : CVSpilsPrecSetupFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner setup function PrecSetup and
 * the user-supplied preconditioner solve function PrecSolve
 * together must define left and right preconditoner matrices
 * P1 and P2 (either of which may be trivial), such that the
 * product P1*P2 is an approximation to the Newton matrix
 * M = I - gamma*J.  Here J is the system Jacobian J = df/dy,
 * and gamma is a scalar proportional to the integration step
 * size h.  The solution of systems P z = r, with P = P1 or P2,
 * is to be carried out by the PrecSolve function, and PrecSetup
 * is to do any necessary setup operations.
 *
 * The user-supplied preconditioner setup function PrecSetup
 * is to evaluate and preprocess any Jacobian-related data
 * needed by the preconditioner solve function PrecSolve.
 * This might include forming a crude approximate Jacobian,
 * and performing an LU factorization on the resulting
 * approximation to M.  This function will not be called in
 * advance of every call to PrecSolve, but instead will be called
 * only as often as necessary to achieve convergence within the
 * Newton iteration.  If the PrecSolve function needs no
 * preparation, the PrecSetup function can be NULL.
 *
 * For greater efficiency, the PrecSetup function may save
 * Jacobian-related data and reuse it, rather than generating it
 * from scratch.  In this case, it should use the input flag jok
 * to decide whether to recompute the data, and set the output
 * flag *jcurPtr accordingly.
 *
 * Each call to the PrecSetup function is preceded by a call to
 * the RhsFn f with the same (t,y) arguments.  Thus the PrecSetup
 * function can use any auxiliary data that is computed and
 * saved by the f function and made accessible to PrecSetup.
 *
 * A function PrecSetup must have the prototype given below.
 * Its parameters are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *          namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == FALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == TRUE  means that Jacobian data, if saved from
 *                  the previous PrecSetup call, can be reused
 *                  (with the current value of gamma).
 *         A Precset call with jok == TRUE can only occur after
 *         a call with jok == FALSE.
 *
 * jcurPtr is a pointer to an output integer flag which is
 *         to be set by PrecSetup as follows:
 *         Set *jcurPtr = TRUE if Jacobian data was recomputed.
 *         Set *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                        but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * user_data  is a pointer to user data - the same as the user_data
 *         parameter passed to the CVodeSetUserData function.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *                      for N_Vectors which can be used by
 *                      CVSpilsPrecSetupFn as temporary storage or
 *                      work space.
 *
 * NOTE: If the user's preconditioner needs other quantities,
 *       they are accessible as follows: hcur (the current stepsize)
 *       and ewt (the error weight vector) are accessible through
 *       CVodeGetCurrentStep and CVodeGetErrWeights, respectively).
 *       The unit roundoff is available as UNIT_ROUNDOFF defined in
 *       sundials_types.h.
 *
 * Returned value:
 * The value to be returned by the PrecSetup function is a flag
 * indicating whether it was successful.  This value should be
 *   0   if successful,
 *   > 0 for a recoverable error (step will be retried),
 *   < 0 for an unrecoverable error (integration is halted).
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
				  booleantype jok, booleantype *jcurPtr,
				  realtype gamma, void *user_data,
				  N_Vector tmp1, N_Vector tmp2,
				  N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * Type : CVSpilsPrecSolveFn
 * -----------------------------------------------------------------
 * The user-supplied preconditioner solve function PrecSolve
 * is to solve a linear system P z = r in which the matrix P is
 * one of the preconditioner matrices P1 or P2, depending on the
 * type of preconditioning chosen.
 *
 * A function PrecSolve must have the prototype given below.
 * Its parameters are as follows:
 *
 * t      is the current value of the independent variable.
 *
 * y      is the current value of the dependent variable vector.
 *
 * fy     is the vector f(t,y).
 *
 * r      is the right-hand side vector of the linear system.
 *
 * z      is the output vector computed by PrecSolve.
 *
 * gamma  is the scalar appearing in the Newton matrix.
 *
 * delta  is an input tolerance for use by PSolve if it uses
 *        an iterative method in its solution.  In that case,
 *        the residual vector Res = r - P z of the system
 *        should be made less than delta in weighted L2 norm,
 *        i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta.
 *        Note: the error weight vector ewt can be obtained
 *        through a call to the routine CVodeGetErrWeights.
 *
 * lr     is an input flag indicating whether PrecSolve is to use
 *        the left preconditioner P1 or right preconditioner
 *        P2: lr = 1 means use P1, and lr = 2 means use P2.
 *
 * user_data is a pointer to user data - the same as the user_data
 *        parameter passed to the CVodeSetUserData function.
 *
 * tmp    is a pointer to memory allocated for an N_Vector
 *        which can be used by PSolve for work space.
 *
 * Returned value:
 * The value to be returned by the PrecSolve function is a flag
 * indicating whether it was successful.  This value should be
 *   0 if successful,
 *   positive for a recoverable error (step will be retried),
 *   negative for an unrecoverable error (integration is halted).
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
				  N_Vector r, N_Vector z,
				  realtype gamma, realtype delta,
				  int lr, void *user_data, N_Vector tmp);

/*
 * -----------------------------------------------------------------
 * Type : CVSpilsJacTimesVecFn
 * -----------------------------------------------------------------
 * The user-supplied function jtimes is to generate the product
 * J*v for given v, where J is the Jacobian df/dy, or an
 * approximation to it, and v is a given vector. It should return
 * 0 if successful a positive value for a recoverable error or 
 * a negative value for an unrecoverable failure.
 *
 * A function jtimes must have the prototype given below. Its
 * parameters are as follows:
 *
 *   v        is the N_Vector to be multiplied by J.
 *
 *   Jv       is the output N_Vector containing J*v.
 *
 *   t        is the current value of the independent variable.
 *
 *   y        is the current value of the dependent variable
 *            vector.
 *
 *   fy       is the vector f(t,y).
 *
 *   user_data   is a pointer to user data, the same as the user_data
 *            parameter passed to the CVodeSetUserData function. 
 *
 *   tmp      is a pointer to memory allocated for an N_Vector
 *            which can be used by Jtimes for work space.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t,
				    N_Vector y, N_Vector fy,
				    void *user_data, N_Vector tmp);


/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVSPILS linear solver
 * -----------------------------------------------------------------
 *
 * CVSpilsSetPrecType resets the type of preconditioner, pretype,
 *                from the value previously set.
 *                This must be one of PREC_NONE, PREC_LEFT, 
 *                PREC_RIGHT, or PREC_BOTH.
 *
 * CVSpilsSetGSType specifies the type of Gram-Schmidt
 *                orthogonalization to be used. This must be one of
 *                the two enumeration constants MODIFIED_GS or
 *                CLASSICAL_GS defined in iterative.h. These correspond
 *                to using modified Gram-Schmidt and classical
 *                Gram-Schmidt, respectively.
 *                Default value is MODIFIED_GS.
 *
 * CVSpilsSetMaxl resets the maximum Krylov subspace size, maxl,
 *                from the value previously set.
 *                An input value <= 0, gives the default value.
 *
 * CVSpilsSetEpsLin specifies the factor by which the tolerance on
 *                the nonlinear iteration is multiplied to get a
 *                tolerance on the linear iteration.
 *                Default value is 0.05.
 *
 * CVSpilsSetPreconditioner specifies the PrecSetup and PrecSolve functions.
 *                Default is NULL for both arguments (no preconditioning).
 *
 * CVSpilsSetJacTimesVecFn specifies the jtimes function. Default is to use 
 *                an internal finite difference approximation routine.
 *
 * The return value of CVSpilsSet* is one of:
 *    CVSPILS_SUCCESS   if successful
 *    CVSPILS_MEM_NULL  if the cvode memory was NULL
 *    CVSPILS_LMEM_NULL if the linear solver memory was NULL
 *    CVSPILS_ILL_INPUT if an input has an illegal value
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVSpilsSetPrecType(void *cvode_mem, int pretype);
SUNDIALS_EXPORT int CVSpilsSetGSType(void *cvode_mem, int gstype);
SUNDIALS_EXPORT int CVSpilsSetMaxl(void *cvode_mem, int maxl);
SUNDIALS_EXPORT int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac);
SUNDIALS_EXPORT int CVSpilsSetPreconditioner(void *cvode_mem,
                                             CVSpilsPrecSetupFn pset, 
					     CVSpilsPrecSolveFn psolve);
SUNDIALS_EXPORT int CVSpilsSetJacTimesVecFn(void *cvode_mem,
                                            CVSpilsJacTimesVecFn jtv);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVSPILS linear solver
 * -----------------------------------------------------------------
 * CVSpilsGetWorkSpace returns the real and integer workspace used
 *                by the SPILS module.
 *
 * CVSpilsGetNumPrecEvals returns the number of preconditioner
 *                 evaluations, i.e. the number of calls made
 *                 to PrecSetup with jok==FALSE.
 *
 * CVSpilsGetNumPrecSolves returns the number of calls made to
 *                 PrecSolve.
 *
 * CVSpilsGetNumLinIters returns the number of linear iterations.
 *
 * CVSpilsGetNumConvFails returns the number of linear
 *                 convergence failures.
 *
 * CVSpilsGetNumJtimesEvals returns the number of calls to jtimes.
 *
 * CVSpilsGetNumRhsEvals returns the number of calls to the user
 *                 f routine due to finite difference Jacobian
 *                 times vector evaluation.
 *
 * CVSpilsGetLastFlag returns the last error flag set by any of
 *                 the CVSPILS interface functions.
 *
 * The return value of CVSpilsGet* is one of:
 *    CVSPILS_SUCCESS   if successful
 *    CVSPILS_MEM_NULL  if the cvode memory was NULL
 *    CVSPILS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVSpilsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals);
SUNDIALS_EXPORT int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves);
SUNDIALS_EXPORT int CVSpilsGetNumLinIters(void *cvode_mem, long int *nliters);
SUNDIALS_EXPORT int CVSpilsGetNumConvFails(void *cvode_mem, long int *nlcfails);
SUNDIALS_EXPORT int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals);
SUNDIALS_EXPORT int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS); 
SUNDIALS_EXPORT int CVSpilsGetLastFlag(void *cvode_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a CVSPILS return flag
 * -----------------------------------------------------------------
 */
  
SUNDIALS_EXPORT char *CVSpilsGetReturnFlagName(long int flag);


/* 
 * -----------------------------------------------------------------
 * PART II - backward problems
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Type : CVSpilsPrecSetupFnB
 * -----------------------------------------------------------------
 * A function PrecSetupB for the adjoint (backward) problem must have 
 * the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsPrecSetupFnB)(realtype t, N_Vector y,
                                   N_Vector yB, N_Vector fyB,
                                   booleantype jokB,
                                   booleantype *jcurPtrB, realtype gammaB,
                                   void *user_dataB,
                                   N_Vector tmp1B, N_Vector tmp2B,
                                   N_Vector tmp3B);


/*
 * -----------------------------------------------------------------
 * Type : CVSpilsPrecSetupFnBS
 * -----------------------------------------------------------------
 * A function PrecSetupBS for the adjoint (backward) problem must have 
 * the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsPrecSetupFnBS)(realtype t, N_Vector y, N_Vector *yS,
                                    N_Vector yB, N_Vector fyB,
                                    booleantype jokB,
                                    booleantype *jcurPtrB, realtype gammaB,
                                    void *user_dataB,
                                    N_Vector tmp1B, N_Vector tmp2B,
                                    N_Vector tmp3B);


/*
 * -----------------------------------------------------------------
 * Type : CVSpilsPrecSolveFnB
 * -----------------------------------------------------------------
 * A function PrecSolveB for the adjoint (backward) problem  must 
 * have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsPrecSolveFnB)(realtype t, N_Vector y,
                                   N_Vector yB, N_Vector fyB,
                                   N_Vector rB, N_Vector zB,
                                   realtype gammaB, realtype deltaB,
                                   int lrB, void *user_dataB, N_Vector tmpB);

/*
 * -----------------------------------------------------------------
 * Type : CVSpilsPrecSolveFnBS
 * -----------------------------------------------------------------
 * A function PrecSolveBS for the adjoint (backward) problem  must 
 * have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsPrecSolveFnBS)(realtype t, N_Vector y, N_Vector *yS,
                                    N_Vector yB, N_Vector fyB,
                                    N_Vector rB, N_Vector zB,
                                    realtype gammaB, realtype deltaB,
                                    int lrB, void *user_dataB, N_Vector tmpB);

/*
 * -----------------------------------------------------------------
 * Type : CVSpilsJacTimesVecFnB
 * -----------------------------------------------------------------
 * A function jtimesB for the adjoint (backward) problem must have 
 * the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsJacTimesVecFnB)(N_Vector vB, N_Vector JvB, realtype t,
                                     N_Vector y, N_Vector yB, N_Vector fyB,
                                     void *jac_dataB, N_Vector tmpB);

/*
 * -----------------------------------------------------------------
 * Type : CVSpilsJacTimesVecFnBS
 * -----------------------------------------------------------------
 * A function jtimesBS for the adjoint (backward) problem must have 
 * the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpilsJacTimesVecFnBS)(N_Vector vB, N_Vector JvB,
                                      realtype t, N_Vector y, N_Vector *yS,
                                      N_Vector yB, N_Vector fyB,
                                      void *jac_dataB, N_Vector tmpB);

/*
 * -----------------------------------------------------------------
 * Functions
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVSpilsSetPrecTypeB(void *cvode_mem, int which, int pretypeB);
SUNDIALS_EXPORT int CVSpilsSetGSTypeB(void *cvode_mem, int which, int gstypeB);
SUNDIALS_EXPORT int CVSpilsSetEpsLinB(void *cvode_mem, int which, realtype eplifacB);
SUNDIALS_EXPORT int CVSpilsSetMaxlB(void *cvode_mem, int which, int maxlB);

SUNDIALS_EXPORT int CVSpilsSetPreconditionerB(void *cvode_mem, int which, 
                                              CVSpilsPrecSetupFnB psetB,
					      CVSpilsPrecSolveFnB psolveB);
SUNDIALS_EXPORT int CVSpilsSetPreconditionerBS(void *cvode_mem, int which, 
                                               CVSpilsPrecSetupFnBS psetBS,
					       CVSpilsPrecSolveFnBS psolveBS);

SUNDIALS_EXPORT int CVSpilsSetJacTimesVecFnB(void *cvode_mem, int which, 
                                             CVSpilsJacTimesVecFnB jtvB);
SUNDIALS_EXPORT int CVSpilsSetJacTimesVecFnBS(void *cvode_mem, int which, 
                                              CVSpilsJacTimesVecFnBS jtvBS);


#ifdef __cplusplus
}
#endif

#endif
