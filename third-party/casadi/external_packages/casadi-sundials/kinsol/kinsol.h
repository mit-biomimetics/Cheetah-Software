/*
 * -----------------------------------------------------------------
 * $Revision: 4423 $
 * $Date: 2015-03-08 17:23:10 -0700 (Sun, 08 Mar 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
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
 * KINSOL solver module header file
 * -----------------------------------------------------------------
 */

#ifndef _KINSOL_H
#define _KINSOL_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *              K I N S O L     C O N S T A N T S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * KINSOL return flags
 * -----------------------------------------------------------------
 */

#define KIN_SUCCESS             0
#define KIN_INITIAL_GUESS_OK    1
#define KIN_STEP_LT_STPTOL      2

#define KIN_WARNING             99

#define KIN_MEM_NULL            -1
#define KIN_ILL_INPUT           -2
#define KIN_NO_MALLOC           -3
#define KIN_MEM_FAIL            -4
#define KIN_LINESEARCH_NONCONV  -5
#define KIN_MAXITER_REACHED     -6
#define KIN_MXNEWT_5X_EXCEEDED  -7
#define KIN_LINESEARCH_BCFAIL   -8
#define KIN_LINSOLV_NO_RECOVERY -9
#define KIN_LINIT_FAIL          -10
#define KIN_LSETUP_FAIL         -11
#define KIN_LSOLVE_FAIL         -12

#define KIN_SYSFUNC_FAIL        -13
#define KIN_FIRST_SYSFUNC_ERR   -14
#define KIN_REPTD_SYSFUNC_ERR   -15


/*
 * -----------------------------------------------------------------
 * Enumeration for inputs to KINSetEtaForm (eta choice)
 * -----------------------------------------------------------------
 * KIN_ETACONSTANT : use constant value for eta (default value is
 *                   0.1 but a different value can be specified via
 *                   a call to KINSetEtaConstValue)
 *
 * KIN_ETACHOICE1 : use choice #1 as given in Eisenstat and Walker's
 *                  paper of SIAM J.Sci.Comput.,17 (1996), pp 16-32,
 *                  wherein eta is defined to be:
 *
 *              eta(k+1) = ABS(||F(u_k+1)||_L2-||F(u_k)+J(u_k)*p_k||_L2)
 *                       ---------------------------------------------
 *                                       ||F(u_k)||_L2
 *
 *                                                      1+sqrt(5)
 *              eta_safe = eta(k)^ealpha where ealpha = ---------
 *                                                          2
 *
 * KIN_ETACHOICE2 : use choice #2 as given in Eisenstat and Walker's
 *                  paper wherein eta is defined to be:
 *
 *                                  [ ||F(u_k+1)||_L2 ]^ealpha
 *              eta(k+1) = egamma * [ --------------- ]
 *                                  [  ||F(u_k)||_L2  ]
 *
 *                  where egamma = [0,1] and ealpha = (1,2]
 *
 *              eta_safe = egamma*(eta(k)^ealpha)
 *
 *                  Note: The default values of the scalar
 *                  coefficients egamma and ealpha (both required)
 *                  are egamma = 0.9 and ealpha = 2.0, but the
 *                  routine KINSetEtaParams can be used to specify
 *                  different values.
 *
 * When using either KIN_ETACHOICE1 or KIN_ETACHOICE2, if
 * eta_safe > 0.1 then the following safeguard is applied:
 *
 *  eta(k+1) = MAX {eta(k+1), eta_safe}
 *
 * The following safeguards are always applied when using either
 * KIN_ETACHOICE1 or KIN_ETACHOICE2 so that eta_min <= eta <= eta_max:
 *
 *  eta(k+1) = MAX {eta(k+1), eta_min}
 *  eta(k+1) = MIN {eta(k+1), eta_max}
 *
 * where eta_min = 1.0e-4 and eta_max = 0.9 (see KINForcingTerm).
 * -----------------------------------------------------------------
 */

#define KIN_ETACHOICE1  1
#define KIN_ETACHOICE2  2
#define KIN_ETACONSTANT 3
  
/*
 * -----------------------------------------------------------------
 * Enumeration for global strategy
 * -----------------------------------------------------------------
 */
  
#define KIN_NONE       0
#define KIN_LINESEARCH 1
#define KIN_PICARD     2
#define KIN_FP         3

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type : KINSysFn
 * -----------------------------------------------------------------
 * The user-supplied subroutine implementing the nonlinear system
 * function (vector-valued function) F must take as input the
 * dependent variable vector uu (type N_Vector), and set fval (type
 * N_Vector) equal to F(uu) before returning. Additional workspace
 * is allocated by the user and referenced by the user_data memory
 * pointer.
 * 
 * Note: The user-defined routine (internally referenced by a
 * a pointer (type KINSysFn) named func) should have an 'int' return
 * value type. However, the return value is currently ignored.
 * -----------------------------------------------------------------
 */

typedef int (*KINSysFn)(N_Vector uu, N_Vector fval, void *user_data );


/*
 * -----------------------------------------------------------------
 * Type : KINErrHandlerFn
 * -----------------------------------------------------------------
 * A function eh, which handles error messages, must have type
 * KINErrHandlerFn.
 * The function eh takes as input the error code, the name of the
 * module reporting the error, the error message, and a pointer to
 * user data, the same as that passed to KINSetUserData.
 * 
 * All error codes are negative, except KIN_WARNING which indicates 
 * a warning (the solver continues).
 *
 * A KINErrHandlerFn has no return value.
 * -----------------------------------------------------------------
 */

typedef void (*KINErrHandlerFn)(int error_code, 
				const char *module, const char *function, 
				char *msg, void *user_data); 


/*
 * -----------------------------------------------------------------
 * Type : KINInfoHandlerFn
 * -----------------------------------------------------------------
 * A function ih, which handles info messages, must have type
 * KINInfoHandlerFn.
 * The function ih takes as input the name of the module and of the
 * function reporting the info message and a pointer to
 * user data, the same as that passed to KINSetfdata.
 * 
 * A KINInfoHandlerFn has no return value.
 * -----------------------------------------------------------------
 */

typedef void (*KINInfoHandlerFn)(const char *module, const char *function, 
				 char *msg, void *user_data); 

/*
 * ================================================================
 *          U S E R - C A L L A B L E   R O U T I N E S           
 * ================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : KINCreate
 * -----------------------------------------------------------------
 * KINCreate allocates and initializes an internal memory block for
 * the KINSOL solver module.
 *
 * If successful, KINCreate returns a pointer to the initialized
 * memory block which should be passed to KINInit. If an
 * error occurs, then KINCreate returns a NULL pointer.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void *KINCreate(void);

/*
 * -----------------------------------------------------------------
 * Optional Input Specification Functions (KINSOL)
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs:
 *
 *     Function Name      |    Optional Input  [Default Value]
 *                        |
 * -----------------------------------------------------------------
 *                        |
 * KINSetErrHandlerFn     | user-provided ErrHandler function.
 *                        | [internal]
 *                        |
 * KINSetErrFile          | pointer (type FILE) indicating where all
 *                        | warning/error messages should be sent
 *                        | if the default internal error handler 
 *                        | is used
 *                        | [stderr]
 *                        |
 * KINSetPrintLevel       | level of verbosity of output:
 *                        |
 *                        |  0  no statistical information is
 *                        |     displayed (default level)
 *                        |
 *                        |  1  for each nonlinear iteration display
 *                        |     the following information: the scaled
 *                        |     norm (L2) of the system function
 *                        |     evaluated at the current iterate, the
 *                        |     scaled norm of the Newton step (only if
 *                        |     using KIN_NONE), and the
 *                        |     number of function evaluations performed
 *                        |     thus far
 *                        |
 *                        |  2  display level 1 output and the
 *                        |     following values for each iteration:
 *                        |
 *                        |       fnorm (L2) = ||fscale*func(u)||_L2
 *                        |       (only for KIN_NONE)
 *                        |
 *                        |       scaled fnorm (for stopping) =
 *                        |       ||fscale*ABS(func(u))||_L-infinity
 *                        |       (for KIN_NONE and
 *                        |       KIN_LINESEARCH)
 *                        |
 *                        |  3  display level 2 output plus additional
 *                        |     values used by the global strategy
 *                        |     (only if using KIN_LINESEARCH), and
 *                        |     statistical information for the linear
 *                        |     solver
 *                        | [0]
 *                        |
 * KINSetInfoHandlerFn    | user-provided InfoHandler function.
 *                        | [internal]
 *                        |
 * KINSetInfoFile         | pointer (type FILE) specifying where
 *                        | informative (non-error) messages should
 *                        | be sent if the default internal info
 *                        | handler is used
 *                        | [stdout]
 *                        |
 * KINSetUserData         | pointer to user-allocated memory that is
 *                        | passed to the user-supplied subroutine
 *                        | implementing the nonlinear system function
 *                        | F(u)
 *                        | [NULL]
 *                        |
 * KINSetNumMaxIters      | maximum number of nonlinear iterations
 *                        | [MXITER_DEFAULT] (defined in kinsol_impl.h)
 *                        |
 * KINSetNoInitSetup      | flag controlling whether or not the
 *                        | KINSol routine makes an initial call
 *                        | to the linear solver setup routine (lsetup)
 *                        | (possible values are TRUE and FALSE)
 *                        | [FALSE]
 *                        |
 * KINSetNoResMon         | flag controlling whether or not the nonlinear
 *                        | residual monitoring scheme is used to control
 *                        | Jacobian updating (possible values are TRUE
 *                        | and FALSE)
 *                        | [FALSE if using direct linear solver]
 *                        | [TRUE if using inexact linear solver]
 *                        |
 * KINSetMaxSetupCalls    | mbset, number of nonlinear iteraions, such 
 *                        | that a call to the linear solver setup routine
 *                        | (lsetup) is forced every mbset iterations.
 *                        | If mbset=1, lsetup s called at every iteration.
 *                        | [MSBSET_DEFAULT] (defined in kinsol_impl.h)
 *                        |
 * KINSetMaxSubSetupCalls | mbsetsub is the number of nonlinear iterations
 *                        | between checks by the nonlinear residual
 *                        | monitoring algorithm (specifies length of
 *                        | subinterval)
 *                        | NOTE: mbset should be a multiple of mbsetsub
 *                        | [MSBSET_SUB_DEFAULT] (defined in kinsol_impl.h)
 *                        |
 * KINSetEtaForm          | flag indicating which method to use to
 *                        | compute the value of the eta coefficient
 *                        | used in the calculation of the linear
 *                        | solver convergence tolerance:
 *                        |
 *                        |  eps = (eta+uround)*||fscale*func(u)||_L2
 *                        |
 *                        | the linear solver tests for convergence by
 *                        | checking if the following inequality has
 *                        | been satisfied:
 *                        |
 *                        |  ||fscale*(func(u)+J(u)*p)||_L2 <= eps
 *                        |
 *                        | where J(u) is the system Jacobian
 *                        | evaluated at the current iterate, and p
 *                        | denotes the Newton step
 *                        |
 *                        | choices for computing eta are as follows:
 *                        |
 *                        |  KIN_ETACHOICE1  (refer to KINForcingTerm)
 *                        |
 *                        |  eta = ABS(||F(u_k+1)||_L2-||F(u_k)+J(u_k)*p_k||_L2)
 *                        |        ---------------------------------------------
 *                        |                        ||F(u_k)||_L2
 *                        | 
 *                        |  KIN_ETACHOICE2  (refer to KINForcingTerm)
 *                        |
 *                        |                [ ||F(u_k+1)||_L2 ]^alpha
 *                        |  eta = gamma * [ --------------- ]
 *                        |                [  ||F(u_k)||_L2  ]
 *                        |
 *                        |  where gamma = [0,1] and alpha = (1,2]
 *                        |
 *                        |  KIN_ETACONSTANT  use a constant value for eta
 *                        | [KIN_ETACHOICE1]
 *                        |
 * KINSetEtaConstValue    | constant value of eta - use with
 *                        | KIN_ETACONSTANT option
 *                        | [0.1]
 *                        |
 * KINSetEtaParams        | values of eta_gamma (egamma) and eta_alpha
 *                        | (ealpha) coefficients - use with KIN_ETACHOICE2
 *                        | option
 *                        | [0.9 and 2.0]
 *                        |
 * KINSetResMonParams     | values of omega_min and omega_max scalars
 *                        | used by nonlinear residual monitoring
 *                        | algorithm (see KINStop)
 *                        | [0.00001 and 0.9]
 *                        |
 * KINSetResMonConstValue | constant value used by residual monitoring
 *                        | algorithm. If omega=0, then it is estimated
 *                        | using omega_min and omega_max.
 *                        | [0.0]
 *                        |
 * KINSetNoMinEps         | flag controlling whether or not the value
 *                        | of eps is bounded below by 0.01*fnormtol
 *                        | (see KINSetFuncNormTol)
 *                        |
 *                        |  FALSE  constrain value of eps by setting
 *                        |         to the following:
 *                        |
 *                        |          eps = MAX{0.01*fnormtol, eps}
 *                        |
 *                        |  TRUE  do not constrain value of eps
 *                        | [FALSE]
 *                        |
 * KINSetMaxNewtonStep    | maximum scaled length of Newton step
 *                        | (reset to value of one if user-supplied
 *                        | value is less than one)
 *                        | [1000*||uscale*u_0||_L2]
 *                        |
 * KINSetMaxBetaFails     | maximum number of beta condition failures
 *                        | in the line search algorithm.
 *                        | [MXNBCF_DEFAULT] (defined in kinsol_impl.h)
 *                        |
 * KINSetRelErrFunc       | real scalar equal to realative error in
 *                        | computing F(u) (used in difference-
 *                        | quotient approximation of matrix-vector
 *                        | product J(u)*v)
 *                        | [(uround)^1/2]
 *                        |
 * KINSetFuncNormTol      | real scalar used as stopping tolerance on
 *                        | ||fscale*ABS(func(u))||_L-infinity (see
 *                        | KINStop and KINInitialStop)
 *                        | [(uround)^1/3]
 *                        |
 * KINSetScaledStepTol    | real scalar used as stopping tolerance on
 *                        | the maximum scaled step length:
 *                        |
 *                        |  ||    u_k+1 - u_k    ||
 *                        |  || ----------------- ||_L-infinity
 *                        |  || ABS(u_k+1)+uscale ||
 *                        |
 *                        | (see KINStop)
 *                        | [(uround)^2/3]
 *                        |
 * KINSetConstraints      | pointer to an array (type N_Vector) of
 *                        | constraints on the solution vector u
 *                        | 
 *                        | if constraints[i] =
 *                        |
 *                        |   0  u[i] not constrained
 *                        |
 *                        |  +1  u[i] constrained to be >= 0
 *                        |  -1  u[i] constrained to be <= 0
 *                        |
 *                        |  +2  u[i] constrained to be > 0
 *                        |  -2  u[i] constrained to be < 0
 *                        |
 *                        | if a NULL pointer is given, then no
 *                        | constraints are applied to vector u
 *                        | [NULL]
 *                        |
 * KINSetSysFunc          | set the user-provided routine which
 *                        | defines the nonlinear problem to be 
 *                        | solved
 *                        | [none]
 * -----------------------------------------------------------------
 * The possible return values for the KINSet* subroutines are the
 * following:
 *
 * KIN_SUCCESS : means the associated variable was successfully
 *               set [0]
 *
 * KIN_MEM_NULL : means a NULL KINSOL memory block pointer was given
 *                (must call the KINCreate and KINInit memory
 *                allocation subroutines prior to calling KINSol) [-1]
 *
 * KIN_ILL_INPUT : means the supplied parameter was invalid (check
 *                 error message) [-2]
 * -----------------------------------------------------------------
 * Note: If successful, these functions return KIN_SUCCESS. If an
 * argument has an illegal value, then an error message is printed
 * to the file specified by errfp and an error code is returned.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSetErrHandlerFn(void *kinmem, KINErrHandlerFn ehfun, void *eh_data);
SUNDIALS_EXPORT int KINSetErrFile(void *kinmem, FILE *errfp);
SUNDIALS_EXPORT int KINSetInfoHandlerFn(void *kinmem, KINInfoHandlerFn ihfun, void *ih_data);
SUNDIALS_EXPORT int KINSetInfoFile(void *kinmem, FILE *infofp);
SUNDIALS_EXPORT int KINSetUserData(void *kinmem, void *user_data);
SUNDIALS_EXPORT int KINSetPrintLevel(void *kinmemm, int printfl);
SUNDIALS_EXPORT int KINSetMAA(void *kinmem, long int maa);
  /* SUNDIALS_EXPORT int KINSetAAStopCrit(void *kinmem, booleantype setstop); */
SUNDIALS_EXPORT int KINSetNumMaxIters(void *kinmem, long int mxiter);
SUNDIALS_EXPORT int KINSetNoInitSetup(void *kinmem, booleantype noInitSetup);
SUNDIALS_EXPORT int KINSetNoResMon(void *kinmem, booleantype noNNIResMon);
SUNDIALS_EXPORT int KINSetMaxSetupCalls(void *kinmem, long int msbset);
SUNDIALS_EXPORT int KINSetMaxSubSetupCalls(void *kinmem, long int msbsetsub);
SUNDIALS_EXPORT int KINSetEtaForm(void *kinmem, int etachoice);
SUNDIALS_EXPORT int KINSetEtaConstValue(void *kinmem, realtype eta);
SUNDIALS_EXPORT int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha);
SUNDIALS_EXPORT int KINSetResMonParams(void *kinmem, realtype omegamin, realtype omegamax);
SUNDIALS_EXPORT int KINSetResMonConstValue(void *kinmem, realtype omegaconst);
SUNDIALS_EXPORT int KINSetNoMinEps(void *kinmem, booleantype noMinEps);
SUNDIALS_EXPORT int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep);
SUNDIALS_EXPORT int KINSetMaxBetaFails(void *kinmem, long int mxnbcf);
SUNDIALS_EXPORT int KINSetRelErrFunc(void *kinmem, realtype relfunc);
SUNDIALS_EXPORT int KINSetFuncNormTol(void *kinmem, realtype fnormtol);
SUNDIALS_EXPORT int KINSetScaledStepTol(void *kinmem, realtype scsteptol);
SUNDIALS_EXPORT int KINSetConstraints(void *kinmem, N_Vector constraints);
SUNDIALS_EXPORT int KINSetSysFunc(void *kinmem, KINSysFn func);

/*
 * -----------------------------------------------------------------
 * Function : KINInit
 * -----------------------------------------------------------------
 * KINInit allocates additional memory for vector storage and
 * sets a couple problem-specific KINSOL variables.
 *
 * Note: Additional vectors must be initialized by the user and
 * passed to the KINSol routine.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  func  name of user-supplied subroutine implementing the
 *        nonlinear function F(u)
 *
 *  tmpl  implementation-specific template vector (type N_Vector)
 *        (created using either N_VNew_Serial or N_VNew_Parallel)
 *
 * KINInit return flags: KIN_SUCCESS, KIN_MEM_NULL, KIN_ILL_INPUT,
 * and KIN_MEM_FAIL (see below). If an error occurs, then KINInit
 * prints an error message.
 *
 * -----------------------------------------------------------------
 * The possible return values for the KINInit subroutine are the
 * following:
 *
 * KIN_SUCCESS : means the necessary system memory was successfully
 *               allocated [0]
 *
 * KIN_MEM_NULL : means a NULL KINSOL memory block pointer was given
 *                (must call the KINCreate routine before calling
 *                KINInit) [-1]
 *
 * KIN_ILL_INPUT : means the name of a user-supplied subroutine
 *                 implementing the nonlinear system function F(u)
 *                 was not given [-2]
 *
 * KIN_MEM_FAIL : means an error occurred during memory allocation
 *                (either insufficient system resources are available
 *                or the vector kernel has not yet been initialized)
 *                [-4]
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINInit(void *kinmem, KINSysFn func, N_Vector tmpl);

/*
 * -----------------------------------------------------------------
 * Function : KINSol
 * -----------------------------------------------------------------
 * KINSol (main KINSOL driver routine) manages the computational
 * process of computing an approximate solution of the nonlinear
 * system. If the initial guess (initial value assigned to vector u)
 * doesn't violate any user-defined constraints, then the subroutine
 * attempts to solve the system F(u) = 0 using a nonlinear Krylov
 * subspace projection method. The Newton-Krylov iterations are
 * stopped if either of the following conditions is satisfied:
 *
 *  ||F(u)||_L-infinity <= 0.01*fnormtol
 *
 *  ||u[i+1] - u[i]||_L-infinity <= scsteptol
 *
 * However, if the current iterate satisfies the second stopping
 * criterion, it doesn't necessarily mean an approximate solution
 * has been found since the algorithm may have stalled, or the
 * user-specified step tolerance (scsteptol) may be too large.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  uu  vector set to initial guess by user before calling KINSol,
 *      but which upon return contains an approximate solution of
 *      the nonlinear system F(u) = 0
 *
 *  strategy  global strategy applied to Newton step if unsatisfactory
 *            (KIN_NONE or KIN_LINESEARCH)
 *
 *  u_scale  vector containing diagonal elements of scaling matrix
 *           for vector u chosen so that the components of
 *           u_scale*u (as a matrix multiplication) all have
 *           about the same magnitude when u is close to a root
 *           of F(u)
 *
 *  f_scale  vector containing diagonal elements of scaling matrix
 *           for F(u) chosen so that the components of
 *           f_scale*F(u) (as a matrix multiplication) all have
 *           roughly the same magnitude when u is not too near a
 *           root of F(u)
 *
 * Note: The components of vectors u_scale and f_scale should be
 * positive.
 *
 * If successful, KINSol returns a vector uu contains an approximate
 * solution of the given nonlinear system. If an error occurs, then
 * an error message is printed and an error code is returned.
 *
 * -----------------------------------------------------------------
 * KINSol Return Values
 * -----------------------------------------------------------------
 *
 * The possible return values for the KINSol subroutine are the
 * following:
 *
 * KIN_SUCCESS : means ||fscale*ABS(func(u))||_L-infinity <= 0.01*fnormtol
 *               and the current iterate uu is probably an approximate
 *               solution of the nonlinear system F(u) = 0 [0]
 *
 * KIN_INITIAL_GUESS_OK : means the initial user-supplied guess
 *                        already satisfies the stopping criterion
 *                        given above [1]
 *
 * KIN_STEP_LT_STPTOL : means the following inequality has been
 *                      satisfied (stopping tolerance on scaled
 *                      step length):
 *
 *                    ||    u_k+1 - u_k    ||
 *                    || ----------------- ||_L-infinity <= scsteptol
 *                    || ABS(u_k+1)+uscale ||
 *
 *                      so the current iterate (denoted above by u_k+1)
 *                      may be an approximate solution of the given
 *                      nonlinear system, but it is also quite possible
 *                      that the algorithm is "stalled" (making
 *                      insufficient progress) near an invalid solution,
 *                      or the real scalar scsteptol is too large [2]
 *
 * KIN_LINESEARCH_NONCONV : means the line search algorithm was unable
 *                          to find an iterate sufficiently distinct
 *                          from the current iterate
 *
 *                          failure to satisfy the sufficient decrease
 *                          condition could mean the current iterate is
 *                          "close" to an approximate solution of the given
 *                          nonlinear system, the finite-difference
 *                          approximation of the matrix-vector product
 *                          J(u)*v is inaccurate, or the real scalar
 *                          scsteptol is too large [-5]
 *
 * KIN_MAXITER_REACHED : means the maximum number of nonlinear iterations
 *                       has been reached [-6]
 *
 * KIN_MXNEWT_5X_EXCEEDED : means five consecutive steps have been taken
 *                          that satisfy the following inequality:
 *
 *                            ||uscale*p||_L2 > 0.99*mxnewtstep
 *
 *                          where p denotes the current step and
 *                          mxnewtstep is a real scalar upper bound
 *                          on the scaled step length
 *
 *                          such a failure may mean ||fscale*func(u)||_L2
 *                          asymptotes from above to a finite value, or
 *                          the real scalar mxnewtstep is too small [-7]
 *
 * KIN_LINESEARCH_BCFAIL : means the line search algorithm (implemented
 *                         in KINLineSearch) was unable to satisfy the
 *                         beta-condition for MXNBCF + 1 nonlinear
 *                         iterations (not necessarily consecutive),
 *                         which may indicate the algorithm is making
 *                         poor progress [-8]
 *
 * KIN_LINSOLV_NO_RECOVERY : means the user-supplied routine psolve
 *                           encountered a recoverable error, but
 *                           the preconditioner is already current [-9]
 *
 * KIN_LINIT_FAIL : means the linear solver initialization routine (linit)
 *                  encountered an error [-10]
 *
 * KIN_LSETUP_FAIL : means the user-supplied routine pset (used to compute
 *                   the preconditioner) encountered an unrecoverable
 *                   error [-11]
 *
 * KIN_LSOLVE_FAIL : means either the user-supplied routine psolve (used to
 *                   to solve the preconditioned linear system) encountered
 *                   an unrecoverable error, or the linear solver routine
 *                   (lsolve) encountered an error condition [-12]
 *
 * KIN_MEM_NULL : means a NULL KINSOL memory block pointer was given
 *                (must call the KINCreate and KINInit memory
 *                allocation subroutines prior to calling KINSol) [-1]
 *
 * KIN_NO_MALLOC : means additional system memory has not yet been
 *                 allocated for vector storage (forgot to call the
 *                 KINInit routine) [-3]
 *
 * KIN_ILL_INPUT : means at least one input parameter was invalid
 *                 (check error output) [-2]
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSol(void *kinmem, N_Vector uu, int strategy,
			   N_Vector u_scale, N_Vector f_scale);

/*
 * -----------------------------------------------------------------
 * Optional Output Extraction Functions (KINSOL)
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistical information related to the KINSOL solver:
 *
 *       Function Name       |      Returned Value
 *                           |
 * -----------------------------------------------------------------
 *                           |
 * KINGetWorkSpace           | returns both integer workspace size
 *                           | (total number of long int-sized blocks
 *                           | of memory allocated by KINSOL for
 *                           | vector storage) and real workspace
 *                           | size (total number of realtype-sized
 *                           | blocks of memory allocated by KINSOL
 *                           | for vector storage)
 *                           |
 * KINGetNumFuncEvals        | total number evaluations of the
 *                           | nonlinear system function F(u)
 *                           | (number of direct calls made to the
 *                           | user-supplied subroutine by KINSOL
 *                           | module member functions)
 *                           |
 * KINGetNumNonlinSolvIters  | total number of nonlinear iterations
 *                           | performed
 *                           |
 * KINGetNumBetaCondFails    | total number of beta-condition
 *                           | failures (see KINLineSearch)
 *                           |
 *                           | KINSOL halts if the number of such
 *                           | failures exceeds the value of the
 *                           | constant MXNBCF (defined in kinsol.c)
 *                           |
 * KINGetNumBacktrackOps     | total number of backtrack operations
 *                           | (step length adjustments) performed
 *                           | by the line search algorithm (see
 *                           | KINLineSearch)
 *                           |
 * KINGetFuncNorm            | scaled norm of the nonlinear system
 *                           | function F(u) evaluated at the
 *                           | current iterate:
 *                           |
 *                           |  ||fscale*func(u)||_L2
 *                           |
 * KINGetStepLength          | scaled norm (or length) of the step
 *                           | used during the previous iteration:
 *                           |
 *                           |  ||uscale*p||_L2
 *                           |
 * KINGetReturnFlagName      | returns the name of the constant
 *                           | associated with a KINSOL return flag
 *                           |
 * -----------------------------------------------------------------
 *
 * The possible return values for the KINSet* subroutines are the
 * following:
 *
 * KIN_SUCCESS : means the information was successfully retrieved [0]
 * 
 * KIN_MEM_NULL : means a NULL KINSOL memory block pointer was given
 *                (must call the KINCreate and KINInit memory
 *                allocation subroutines prior to calling KINSol) [-1]
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw);
SUNDIALS_EXPORT int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters);
SUNDIALS_EXPORT int KINGetNumFuncEvals(void *kinmem, long int *nfevals);
SUNDIALS_EXPORT int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails); 
SUNDIALS_EXPORT int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr);
SUNDIALS_EXPORT int KINGetFuncNorm(void *kinmem, realtype *fnorm);
SUNDIALS_EXPORT int KINGetStepLength(void *kinmem, realtype *steplength);
SUNDIALS_EXPORT char *KINGetReturnFlagName(long int flag);

/*
 * -----------------------------------------------------------------
 * Function : KINFree
 * -----------------------------------------------------------------
 * KINFree frees system memory resources reserved for the KINSOL
 * solver module.
 *
 *  kinmem  pointer to an internal memory block allocated during
 *          prior calls to KINCreate and KINInit
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void KINFree(void **kinmem);


#ifdef __cplusplus
}
#endif

#endif
