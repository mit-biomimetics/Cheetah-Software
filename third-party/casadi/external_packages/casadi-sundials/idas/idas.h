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
 * This is the header (include) file for the main IDAS solver.
 * -----------------------------------------------------------------
 *
 * IDAS is used to solve numerically the initial value problem     
 * for the differential algebraic equation (DAE) system           
 *   F(t,y,y') = 0,                                               
 * given initial conditions                                       
 *   y(t0) = y0,   y'(t0) = yp0.                                  
 * Here y and F are vectors of length N.                          
 *
 * Additionally, IDAS can perform forward or adjoint sensitivity
 * analysis.
 * -----------------------------------------------------------------
 */

#ifndef _IDAS_H
#define _IDAS_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

/* * =================================================================
 *              I D A S     C O N S T A N T S
 * =================================================================
 */

/*
 * ----------------------------------------------------------------
 * Inputs to:
 *  IDAInit, IDAReInit, 
 *  IDASensInit, IDASensReInit, 
 *  IDAQuadInit, IDAQuadReInit,
 *  IDAQuadSensInit, IDAQuadSensReInit,
 *  IDACalcIC, IDASolve,
 *  IDAAdjInit
 * ----------------------------------------------------------------
 */

/* itask */
#define IDA_NORMAL           1
#define IDA_ONE_STEP         2

/* icopt */
#define IDA_YA_YDP_INIT      1 
#define IDA_Y_INIT           2

/* ism */
#define IDA_SIMULTANEOUS     1
#define IDA_STAGGERED        2

/* DQtype */
#define IDA_CENTERED         1
#define IDA_FORWARD          2

/* interp */
#define IDA_HERMITE          1
#define IDA_POLYNOMIAL       2

/*
 * ===============================================================
 * IDAS RETURN VALUES
 * ===============================================================
 */

#define IDA_SUCCESS          0
#define IDA_TSTOP_RETURN     1
#define IDA_ROOT_RETURN      2

#define IDA_WARNING          99

#define IDA_TOO_MUCH_WORK   -1
#define IDA_TOO_MUCH_ACC    -2
#define IDA_ERR_FAIL        -3
#define IDA_CONV_FAIL       -4

#define IDA_LINIT_FAIL      -5
#define IDA_LSETUP_FAIL     -6
#define IDA_LSOLVE_FAIL     -7
#define IDA_RES_FAIL        -8
#define IDA_REP_RES_ERR     -9
#define IDA_RTFUNC_FAIL     -10
#define IDA_CONSTR_FAIL     -11

#define IDA_FIRST_RES_FAIL  -12
#define IDA_LINESEARCH_FAIL -13
#define IDA_NO_RECOVERY     -14

#define IDA_MEM_NULL        -20
#define IDA_MEM_FAIL        -21
#define IDA_ILL_INPUT       -22
#define IDA_NO_MALLOC       -23
#define IDA_BAD_EWT         -24
#define IDA_BAD_K           -25
#define IDA_BAD_T           -26
#define IDA_BAD_DKY         -27

#define IDA_NO_QUAD         -30
#define IDA_QRHS_FAIL       -31
#define IDA_FIRST_QRHS_ERR  -32
#define IDA_REP_QRHS_ERR    -33

#define IDA_NO_SENS         -40
#define IDA_SRES_FAIL       -41
#define IDA_REP_SRES_ERR    -42
#define IDA_BAD_IS          -43

#define IDA_NO_QUADSENS     -50
#define IDA_QSRHS_FAIL      -51
#define IDA_FIRST_QSRHS_ERR -52
#define IDA_REP_QSRHS_ERR   -53

/*
 * -----------------------------------------
 * IDAA return flags
 * -----------------------------------------
 */

#define IDA_NO_ADJ          -101
#define IDA_NO_FWD          -102
#define IDA_NO_BCK          -103
#define IDA_BAD_TB0         -104
#define IDA_REIFWD_FAIL     -105
#define IDA_FWD_FAIL        -106
#define IDA_GETY_BADT       -107

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * ----------------------------------------------------------------
 * Type : IDAResFn                                                   
 * ----------------------------------------------------------------
 * The F function which defines the DAE system   F(t,y,y')=0      
 * must have type IDAResFn.                                          
 * Symbols are as follows: 
 *                  t  <-> t        y <-> yy               
 *                  y' <-> yp       F <-> rr
 * A IDAResFn takes as input the independent variable value t,    
 * the dependent variable vector yy, and the derivative (with     
 * respect to t) of the yy vector, yp.  It stores the result of   
 * F(t,y,y') in the vector rr. The yy, yp, and rr arguments are of 
 * type N_Vector. The user_data parameter is the pointer user_data 
 * passed by the user to the IDASetUserData routine. This user-supplied 
 * pointer is passed to the user's res function every time it is called, 
 * to provide access in res to user data.                                    
 *                                                                
 * A IDAResFn res should return a value of 0 if successful, a positive
 * value if a recoverable error occured (e.g. yy has an illegal value),
 * or a negative value if a nonrecoverable error occured. In the latter
 * case, the program halts. If a recoverable error occured, the integrator
 * will attempt to correct and retry.
 * ----------------------------------------------------------------
 */

typedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp,
			N_Vector rr, void *user_data);

/*
 * -----------------------------------------------------------------
 * Type : IDARootFn
 * -----------------------------------------------------------------
 * A function g, which defines a set of functions g_i(t,y,y') whose
 * roots are sought during the integration, must have type IDARootFn.
 * The function g takes as input the independent variable value t,
 * the dependent variable vector y, and its t-derivative yp (= y').
 * It stores the nrtfn values g_i(t,y,y') in the realtype array gout.
 * (Allocation of memory for gout is handled within IDA.)
 * The user_data parameter is the same as that passed by the user
 * to the IDASetUserData routine.  This user-supplied pointer is
 * passed to the user's g function every time it is called.
 *
 * An IDARootFn should return 0 if successful or a non-zero value
 * if an error occured (in which case the integration will be halted).
 * -----------------------------------------------------------------
 */

typedef int (*IDARootFn)(realtype t, N_Vector y, N_Vector yp,
			 realtype *gout, void *user_data);

/*
 * -----------------------------------------------------------------
 * Type : IDAEwtFn
 * -----------------------------------------------------------------
 * A function e, which sets the error weight vector ewt, must have
 * type IDAEwtFn.
 * The function e takes as input the current dependent variable y.
 * It must set the vector of error weights used in the WRMS norm:
 * 
 *   ||y||_WRMS = sqrt [ 1/N * sum ( ewt_i * y_i)^2 ]
 *
 * Typically, the vector ewt has components:
 * 
 *   ewt_i = 1 / (reltol * |y_i| + abstol_i)
 *
 * The user_data parameter is the same as that passed by the user
 * to the IDASetUserData routine.  This user-supplied pointer is
 * passed to the user's e function every time it is called.
 * An IDAEwtFn e must return 0 if the error weight vector has been
 * successfuly set and a non-zero value otherwise.
 * -----------------------------------------------------------------
 */

typedef int (*IDAEwtFn)(N_Vector y, N_Vector ewt, void *user_data);

/*
 * -----------------------------------------------------------------
 * Type : IDAErrHandlerFn
 * -----------------------------------------------------------------
 * A function eh, which handles error messages, must have type
 * IDAErrHandlerFn.
 * The function eh takes as input the error code, the name of the
 * module reporting the error, the error message, and a pointer to
 * user data, the same as that passed to IDASetUserData.
 * 
 * All error codes are negative, except IDA_WARNING which indicates 
 * a warning (the solver continues).
 *
 * An IDAErrHandlerFn has no return value.
 * -----------------------------------------------------------------
 */

typedef void (*IDAErrHandlerFn)(int error_code, 
				const char *module, const char *function, 
				char *msg, void *user_data); 

/*
 * -----------------------------------------------------------------
 * Type : IDAQuadRhsFn
 * -----------------------------------------------------------------
 * The rhsQ function which defines the right hand side of the
 * quadrature equations yQ' = rhsQ(t,y) must have type IDAQuadRhsFn.
 * rhsQ takes as input the value of the independent variable t,
 * the vector of states y and y' and must store the result of rhsQ in
 * rrQ. (Allocation of memory for rrQ is handled by IDAS).
 *
 * The user_data parameter is the same as the user_data parameter
 * set by the user through the IDASetUserData routine and is
 * passed to the rhsQ function every time it is called.
 *
 * A function of type IDAQuadRhsFn should return 0 if successful,
 * a negative value if an unrecoverable error occured, and a positive
 * value if a recoverable error (e.g. invalid y values) occured. 
 * If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) IDAS will
 * try to correct and retry.
 * -----------------------------------------------------------------
 */
  
typedef int (*IDAQuadRhsFn)(realtype tres, 
			    N_Vector yy, N_Vector yp,
			    N_Vector rrQ,
			    void *user_data);

/*
 * -----------------------------------------------------------------
 * Type : IDASensResFn
 * -----------------------------------------------------------------
 * The resS function which defines the right hand side of the
 * sensitivity DAE systems F_y * s + F_y' * s' + F_p = 0 
 * must have type IDASensResFn.
 * 
 * resS takes as input the number of sensitivities Ns, the
 * independent variable value t, the states yy and yp and the
 * corresponding value of the residual in resval, and the dependent
 * sensitivity vectors yyS and ypS. It stores the residual in 
 * resvalS. (Memory allocation for resvalS is handled within IDAS)
 *
 * The user_data parameter is the same as the user_data parameter
 * set by the user through the IDASetUserData routine and is
 * passed to the resS function every time it is called.
 *
 * A IDASensResFn should return 0 if successful, a negative value if
 * an unrecoverable error occured, and a positive value if a 
 * recoverable error (e.g. invalid y, yp, yyS or ypS values) 
 * occured. If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) IDAS will
 * try to correct and retry.
 * -----------------------------------------------------------------
 */

typedef int (*IDASensResFn)(int Ns, realtype t, 
			    N_Vector yy, N_Vector yp, N_Vector resval,
			    N_Vector *yyS, N_Vector *ypS, 
                            N_Vector *resvalS, void *user_data,
			    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * Type : IDAQuadSensRhsFn
 * -----------------------------------------------------------------
 * The rhsQS function which defines the RHS of the sensitivity DAE 
 * systems for quadratures  must have type IDAQuadSensRhsFn.
 *
 * rhsQS takes as input the number of sensitivities Ns (the same as
 * that passed to IDAQuadSensInit), the independent variable 
 * value t, the states yy, yp and the dependent sensitivity vectors 
 * yyS and ypS, as well as the current value of the quadrature RHS 
 * rrQ. It stores the result of rhsQS in rhsvalQS.
 * (Allocation of memory for resvalQS is handled within IDAS)
 *
 * A IDAQuadSensRhsFn should return 0 if successful, a negative
 * value if an unrecoverable error occured, and a positive value
 * if a recoverable error (e.g. invalid yy, yp, yyS or ypS values) 
 * occured. If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) IDAS
 * will try to correct and retry.
 * -----------------------------------------------------------------
 */

typedef int (*IDAQuadSensRhsFn)(int Ns, realtype t,
                               N_Vector yy, N_Vector yp, 
                               N_Vector *yyS, N_Vector *ypS, 
                               N_Vector rrQ, N_Vector *rhsvalQS,
                               void *user_data,
                               N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS);

/*
 * -----------------------------------------------------------------
 * Types: IDAResFnB and IDAResFnBS
 * -----------------------------------------------------------------
 *    The resB function which defines the right hand side of the
 *    DAE systems to be integrated backwards must have type IDAResFnB.
 *    If the backward problem depends on forward sensitivities, its
 *    RHS function must have type IDAResFnBS.
 * -----------------------------------------------------------------
 * Types: IDAQuadRhsFnB and IDAQuadRhsFnBS
 * -----------------------------------------------------------------
 *    The rhsQB function which defines the quadratures to be integrated
 *    backwards must have type IDAQuadRhsFnB.
 *    If the backward problem depends on forward sensitivities, its
 *    quadrature RHS function must have type IDAQuadRhsFnBS.
 * -----------------------------------------------------------------
 */

typedef int (*IDAResFnB)(realtype tt, 
			 N_Vector yy, N_Vector yp,
			 N_Vector yyB, N_Vector ypB, 
                         N_Vector rrB, void *user_dataB);

typedef int (*IDAResFnBS)(realtype t, 
                          N_Vector yy, N_Vector yp, 
                          N_Vector *yyS, N_Vector *ypS,
                          N_Vector yyB, N_Vector ypB,
                          N_Vector rrBS, void *user_dataB);

typedef int (*IDAQuadRhsFnB)(realtype tt, 
                             N_Vector yy, N_Vector yp, 
                             N_Vector yyB, N_Vector ypB,
                             N_Vector rhsvalBQ, void *user_dataB);

typedef int (*IDAQuadRhsFnBS)(realtype t, 
                              N_Vector yy, N_Vector yp,
                              N_Vector *yyS, N_Vector *ypS,
                              N_Vector yyB, N_Vector ypB,
                              N_Vector rhsvalBQS, void *user_dataB);
/*
 * ================================================================
 *          U S E R - C A L L A B L E   R O U T I N E S           
 * ================================================================
 */

/* 
 * ----------------------------------------------------------------
 * Function : IDACreate                                           
 * ----------------------------------------------------------------
 * IDACreate creates an internal memory block for a problem to    
 * be solved by IDA.                                              
 *                                                                
 * If successful, IDACreate returns a pointer to initialized      
 * problem memory. This pointer should be passed to IDAInit.    
 * If an initialization error occurs, IDACreate prints an error   
 * message to standard err and returns NULL.                      
 *                                                                
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT void *IDACreate(void);

/*
 * ----------------------------------------------------------------
 * Integrator optional input specification functions              
 * ----------------------------------------------------------------
 * The following functions can be called to set optional inputs   
 * to values other than the defaults given below:                 
 *                                                                
 *                      |                                         
 * Function             |  Optional input / [ default value ]     
 *                      |                                          
 * ---------------------------------------------------------------- 
 *                      |                                          
 * IDASetErrHandlerFn   | user-provided ErrHandler function.
 *                      | [internal]
 *                      |
 * IDASetErrFile        | the file pointer for an error file
 *                      | where all IDA warning and error
 *                      | messages will be written if the default
 *                      | internal error handling function is used. 
 *                      | This parameter can be stdout (standard 
 *                      | output), stderr (standard error), or a 
 *                      | file pointer (corresponding to a user 
 *                      | error file opened for writing) returned 
 *                      | by fopen.
 *                      | If not called, then all messages will
 *                      | be written to the standard error stream.
 *                      | [stderr]
 *                      |                                          
 * IDASetUserData       | a pointer to user data that will be     
 *                      | passed to all user-supplied functions.
 *                      | [NULL]                                  
 *                      |         
 * IDASetMaxOrd         | maximum lmm order to be used by the     
 *                      | solver.                                 
 *                      | [5]                                      
 *                      |                                          
 * IDASetMaxNumSteps    | maximum number of internal steps to be  
 *                      | taken by the solver in its attempt to   
 *                      | reach tout.                             
 *                      | [500]                                   
 *                      |                                          
 * IDASetInitStep       | initial step size.                      
 *                      | [estimated by IDA]                       
 *                      |                                          
 * IDASetMaxStep        | maximum absolute value of step size     
 *                      | allowed.                                
 *                      | [infinity]                              
 *                      |                                          
 * IDASetStopTime       | the independent variable value past     
 *                      | which the solution is not to proceed.   
 *                      | [infinity]                              
 *                      |                                          
 * IDASetNonlinConvCoef | Newton convergence test  constant       
 *                      | for use during integration.             
 *                      | [0.33]                                  
 *                      |                                          
 * IDASetMaxErrTestFails| Maximum number of error test failures   
 *                      | in attempting one step.                 
 *                      | [10]                                    
 *                      |                                         
 * IDASetMaxNonlinIters | Maximum number of nonlinear solver      
 *                      | iterations at one solution.             
 *                      | [4]                                     
 *                      |                                         
 * IDASetMaxConvFails   | Maximum number of allowable conv.       
 *                      | failures in attempting one step.        
 *                      | [10]                                    
 *                      |                                         
 * IDASetSuppressAlg    | flag to indicate whether or not to      
 *                      | suppress algebraic variables in the     
 *                      | local error tests:                      
 *                      | FALSE = do not suppress;                 
 *                      | TRUE = do suppress;                     
 *                      | [FALSE]                                 
 *                      | NOTE: if suppressed algebraic variables 
 *                      | is selected, the nvector 'id' must be   
 *                      | supplied for identification of those    
 *                      | algebraic components (see IDASetId)    
 *                      |                                          
 * IDASetId             | an N_Vector, which states a given       
 *                      | element to be either algebraic or       
 *                      | differential.                           
 *                      | A value of 1.0 indicates a differential 
 *                      | variable while a 0.0 indicates an       
 *                      | algebraic variable. 'id' is required    
 *                      | if optional input SUPPRESSALG is set,   
 *                      | or if IDACalcIC is to be called with    
 *                      | icopt = IDA_YA_YDP_INIT.               
 *                      |                                         
 * IDASetConstraints    | an N_Vector defining inequality         
 *                      | constraints for each component of the   
 *                      | solution vector y. If a given element   
 *                      | of this vector has values +2 or -2,     
 *                      | then the corresponding component of y   
 *                      | will be constrained to be > 0.0 or      
 *                      | <0.0, respectively, while if it is +1   
 *                      | or -1, the y component is constrained   
 *                      | to be >= 0.0 or <= 0.0, respectively.   
 *                      | If a component of constraints is 0.0,   
 *                      | then no constraint is imposed on the    
 *                      | corresponding component of y.           
 *                      | The presence of a non-NULL constraints  
 *                      | vector that is not 0.0 (ZERO) in all    
 *                      | components will cause constraint        
 *                      | checking to be performed.               
 *                      |                                         
 * -----------------------------------------------------------------
 *                          |
 * IDASetRootDirection      | Specifies the direction of zero
 *                          | crossings to be monitored
 *                          | [both directions]
 *                          |
 * IDASetNoInactiveRootWarn | disable warning about possible
 *                          | g==0 at beginning of integration
 *                          | 
 * ---------------------------------------------------------------- 
 * Return flag:
 *   IDA_SUCCESS   if successful
 *   IDA_MEM_NULL  if the IDAS memory is NULL
 *   IDA_ILL_INPUT if an argument has an illegal value
 *
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASetErrHandlerFn(void *ida_mem, IDAErrHandlerFn ehfun, void *eh_data);
SUNDIALS_EXPORT int IDASetErrFile(void *ida_mem, FILE *errfp);
SUNDIALS_EXPORT int IDASetUserData(void *ida_mem, void *user_data);
SUNDIALS_EXPORT int IDASetMaxOrd(void *ida_mem, int maxord);
SUNDIALS_EXPORT int IDASetMaxNumSteps(void *ida_mem, long int mxsteps);
SUNDIALS_EXPORT int IDASetInitStep(void *ida_mem, realtype hin);
SUNDIALS_EXPORT int IDASetMaxStep(void *ida_mem, realtype hmax);
SUNDIALS_EXPORT int IDASetStopTime(void *ida_mem, realtype tstop);
SUNDIALS_EXPORT int IDASetNonlinConvCoef(void *ida_mem, realtype epcon);
SUNDIALS_EXPORT int IDASetMaxErrTestFails(void *ida_mem, int maxnef);
SUNDIALS_EXPORT int IDASetMaxNonlinIters(void *ida_mem, int maxcor);
SUNDIALS_EXPORT int IDASetMaxConvFails(void *ida_mem, int maxncf);
SUNDIALS_EXPORT int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg);
SUNDIALS_EXPORT int IDASetId(void *ida_mem, N_Vector id);
SUNDIALS_EXPORT int IDASetConstraints(void *ida_mem, N_Vector constraints);

SUNDIALS_EXPORT int IDASetRootDirection(void *ida_mem, int *rootdir);
SUNDIALS_EXPORT int IDASetNoInactiveRootWarn(void *ida_mem);

/*
 * ----------------------------------------------------------------
 * Function : IDAInit                                           
 * ----------------------------------------------------------------
 * IDAInit allocates and initializes memory for a problem to    
 * to be solved by IDAS.                                           
 *                                                                
 * res     is the residual function F in F(t,y,y') = 0.                     
 *                                                                
 * t0      is the initial value of t, the independent variable.   
 *                                                                
 * yy0     is the initial condition vector y(t0).                 
 *                                                                
 * yp0     is the initial condition vector y'(t0)                 
 *                                                                
 *  IDA_SUCCESS if successful
 *  IDA_MEM_NULL if the IDAS memory was NULL
 *  IDA_MEM_FAIL if a memory allocation failed
 *  IDA_ILL_INPUT f an argument has an illegal value.
 *                                                                
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAInit(void *ida_mem, IDAResFn res,
                            realtype t0, N_Vector yy0, N_Vector yp0);

/*
 * ----------------------------------------------------------------
 * Function : IDAReInit                                           
 * ----------------------------------------------------------------
 * IDAReInit re-initializes IDAS for the solution of a problem,    
 * where a prior call to IDAInit has been made.                 
 * IDAReInit performs the same input checking and initializations 
 * that IDAInit does.                                           
 * But it does no memory allocation, assuming that the existing   
 * internal memory is sufficient for the new problem.             
 *                                                                
 * The use of IDAReInit requires that the maximum method order,   
 * maxord, is no larger for the new problem than for the problem  
 * specified in the last call to IDAInit.  This condition is    
 * automatically fulfilled if the default value for maxord is     
 * specified.                                                     
 *                                                                
 * Following the call to IDAReInit, a call to the linear solver   
 * specification routine is necessary if a different linear solver
 * is chosen, but may not be otherwise.  If the same linear solver
 * is chosen, and there are no changes in its input parameters,   
 * then no call to that routine is needed.                        
 *                                                                
 * The first argument to IDAReInit is:                            
 *                                                                
 * ida_mem = pointer to IDA memory returned by IDACreate.         
 *                                                                
 * All the remaining arguments to IDAReInit have names and        
 * meanings identical to those of IDAInit.                      
 *                                                                
 * The return value of IDAReInit is equal to SUCCESS = 0 if there 
 * were no errors; otherwise it is a negative int equal to:       
 *   IDA_MEM_NULL   indicating ida_mem was NULL, or            
 *   IDA_NO_MALLOC  indicating that ida_mem was not allocated. 
 *   IDA_ILL_INPUT  indicating an input argument was illegal   
 *                  (including an attempt to increase maxord). 
 * In case of an error return, an error message is also printed.  
 * ----------------------------------------------------------------
 */                                                                

SUNDIALS_EXPORT int IDAReInit(void *ida_mem,
			      realtype t0, N_Vector yy0, N_Vector yp0);
 
/*
 * -----------------------------------------------------------------
 * Functions : IDASStolerances
 *             IDASVtolerances
 *             IDAWFtolerances
 * -----------------------------------------------------------------
 *
 * These functions specify the integration tolerances. One of them
 * MUST be called before the first call to IDA.
 *
 * IDASStolerances specifies scalar relative and absolute tolerances.
 * IDASVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance (a potentially different absolute tolerance 
 *   for each vector component).
 * IDAWFtolerances specifies a user-provides function (of type IDAEwtFn)
 *   which will be called to set the error weight vector.
 *
 * The tolerances reltol and abstol define a vector of error weights,
 * ewt, with components
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol)      (in the SS case), or
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (in the SV case).
 * This vector is used in all error and convergence tests, which
 * use a weighted RMS norm on all error-like vectors v:
 *    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
 * where N is the problem dimension.
 *
 * The return value of these functions is equal to IDA_SUCCESS = 0 if
 * there were no errors; otherwise it is a negative int equal to:
 *   IDa_MEM_NULL     indicating ida_mem was NULL (i.e.,
 *                    IDACreate has not been called).
 *   IDA_NO_MALLOC    indicating that ida_mem has not been
 *                    allocated (i.e., IDAInit has not been
 *                    called).
 *   IDA_ILL_INPUT    indicating an input argument was illegal
 *                    (e.g. a negative tolerance)
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASStolerances(void *ida_mem, realtype reltol, realtype abstol);
SUNDIALS_EXPORT int IDASVtolerances(void *ida_mem, realtype reltol, N_Vector abstol);
SUNDIALS_EXPORT int IDAWFtolerances(void *ida_mem, IDAEwtFn efun);

/* ----------------------------------------------------------------
 * Initial Conditions optional input specification functions      
 * ----------------------------------------------------------------
 * The following functions can be called to set optional inputs   
 * to control the initial conditions calculations.                
 *                                                                
 *                        |                                        
 * Function               |  Optional input / [ default value ]   
 *                        |                                        
 * -------------------------------------------------------------- 
 *                        |                                        
 * IDASetNonlinConvCoefIC | positive coeficient in the Newton     
 *                        | convergence test.  This test uses a   
 *                        | weighted RMS norm (with weights       
 *                        | defined by the tolerances, as in      
 *                        | IDASolve).  For new initial value     
 *                        | vectors y and y' to be accepted, the  
 *                        | norm of J-inverse F(t0,y,y') is       
 *                        | required to be less than epiccon,     
 *                        | where J is the system Jacobian.       
 *                        | [0.01 * 0.33]                          
 *                        |                                        
 * IDASetMaxNumStepsIC    | maximum number of values of h allowed 
 *                        | when icopt = IDA_YA_YDP_INIT, where  
 *                        | h appears in the system Jacobian,     
 *                        | J = dF/dy + (1/h)dF/dy'.              
 *                        | [5]                                   
 *                        |                                        
 * IDASetMaxNumJacsIC     | maximum number of values of the       
 *                        | approximate Jacobian or preconditioner
 *                        | allowed, when the Newton iterations   
 *                        | appear to be slowly converging.       
 *                        | [4]                                    
 *                        |                                        
 * IDASetMaxNumItersIC    | maximum number of Newton iterations   
 *                        | allowed in any one attempt to solve   
 *                        | the IC problem.                       
 *                        | [10]                                  
 *                        |                                        
 * IDASetLineSearchOffIC  | a boolean flag to turn off the        
 *                        | linesearch algorithm.                 
 *                        | [FALSE]                               
 *                        |                                        
 * IDASetStepToleranceIC  | positive lower bound on the norm of   
 *                        | a Newton step.                        
 *                        | [(unit roundoff)^(2/3)                
 *                                                                
 * ---------------------------------------------------------------- 
 * Return flag:
 *   IDA_SUCCESS   if successful
 *   IDA_MEM_NULL  if the IDAS memory is NULL
 *   IDA_ILL_INPUT if an argument has an illegal value
 *
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASetNonlinConvCoefIC(void *ida_mem, realtype epiccon);
SUNDIALS_EXPORT int IDASetMaxNumStepsIC(void *ida_mem, int maxnh);
SUNDIALS_EXPORT int IDASetMaxNumJacsIC(void *ida_mem, int maxnj);
SUNDIALS_EXPORT int IDASetMaxNumItersIC(void *ida_mem, int maxnit);
SUNDIALS_EXPORT int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff);
SUNDIALS_EXPORT int IDASetStepToleranceIC(void *ida_mem, realtype steptol);

/*
 * -----------------------------------------------------------------
 * Function : IDARootInit
 * -----------------------------------------------------------------
 * IDARootInit initializes a rootfinding problem to be solved
 * during the integration of the DAE system.  It must be called
 * after IDACreate, and before IDASolve.  The arguments are:
 *
 * ida_mem = pointer to IDA memory returned by IDACreate.
 *
 * nrtfn   = number of functions g_i, an int >= 0.
 *
 * g       = name of user-supplied function, of type IDARootFn,
 *           defining the functions g_i whose roots are sought.
 *
 * If a new problem is to be solved with a call to IDAReInit,
 * where the new problem has no root functions but the prior one
 * did, then call IDARootInit with nrtfn = 0.
 *
 * The return value of IDARootInit is IDA_SUCCESS = 0 if there were
 * no errors; otherwise it is a negative int equal to:
 *   IDA_MEM_NULL     indicating ida_mem was NULL, or
 *   IDA_MEM_FAIL     indicating a memory allocation failed.
 *                    (including an attempt to increase maxord).
 *   IDA_ILL_INPUT    indicating nrtfn > 0 but g = NULL.
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g);

/*
 * -----------------------------------------------------------------
 * Quadrature optional input specification functions
 * -----------------------------------------------------------------
 * The following function can be called to set optional inputs
 * to values other than the defaults given below:
 *
 * Function             |  Optional input / [ default value ]
 * --------------------------------------------------------------
 *                      |
 * IDASetQuadErrCon     | are quadrature variables considered in
 *                      | the error control?
 *                      | If yes, set tolerances for quadrature
 *                      | integration. 
 *                      | [errconQ = FALSE]
 *                      |
 * -----------------------------------------------------------------
 * If successful, the function return IDA_SUCCESS. If an argument
 * has an illegal value, they print an error message to the
 * file specified by errfp and return one of the error flags
 * defined for the IDASet* routines.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASetQuadErrCon(void *ida_mem, booleantype errconQ);

/*
 * ----------------------------------------------------------------
 * Function : IDAQuadInit and IDAQuadReInit                                     
 * ----------------------------------------------------------------
 * IDAQuadInit allocates and initializes memory related to      
 * quadrature integration.                                        
 *
 * IDAQuadReInit re-initializes IDAS's quadrature related         
 * memory for a problem, assuming it has already been allocated   
 * in prior calls to IDAInit and IDAQuadInit. 
 *                                                                
 * ida_mem is a pointer to IDAS memory returned by IDACreate      
 *                                                                
 * rhsQ  is the user-provided integrand routine.                  
 *                                                                
 * yQ0   is a pointer to a vector specification structure         
 *       for N_Vectors containing quadrature variables.           
 *                                                                
 *                                                                
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAQuadInit(void *ida_mem, IDAQuadRhsFn rhsQ, N_Vector yQ0);
SUNDIALS_EXPORT int IDAQuadReInit(void *ida_mem, N_Vector yQ0);

/*
 * -----------------------------------------------------------------
 * Functions : IDAQuadSStolerances
 *             IDAQuadSVtolerances
 * -----------------------------------------------------------------
 *
 * These functions specify the integration tolerances for quadrature
 * variables. One of them MUST be called before the first call to
 * IDA IF error control on the quadrature variables is enabled
 * (see IDASetQuadErrCon).
 *
 * IDASStolerances specifies scalar relative and absolute tolerances.
 * IDASVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance (a potentially different absolute tolerance 
 *   for each vector component).
 *
 * Return values:
 *  IDA_SUCCESS    if successful
 *  IDA_MEM_NULL   if the solver memory was NULL
 *  IDA_NO_QUAD    if quadratures were not initialized
 *  IDA_ILL_INPUT  if an input argument was illegal
 *                 (e.g. a negative tolerance)
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAQuadSStolerances(void *ida_mem, realtype reltolQ, realtype abstolQ);
SUNDIALS_EXPORT int IDAQuadSVtolerances(void *ida_mem, realtype reltolQ, N_Vector abstolQ);

/* 
 * ----------------------------------------------------------------
 * Forward sensitivity optional input specification functions     
 * ----------------------------------------------------------------
 * The following functions can be called to set optional inputs   
 * to other values than the defaults given below:                 
 *                                                                
 * Function                 |  Optional input / [ default value ]     
 *                          |                                          
 * -------------------------------------------------------------- 
 *                          |                                         
 * IDASetSensDQMethod       | controls the selection of finite
 *                          | difference schemes used in evaluating
 *                          | the sensitivity right hand sides:
 *                          | (centered vs. forward and 
 *                          | simultaneous vs. separate)
 *                          | [DQtype=IDA_CENTERED]
 *                          | [DQrhomax=0.0]                                   
 *                          |                                         
 * IDASetSensParams         |   parameter information:
 *                          | p: pointer to problem parameters
 *                          | plist: list of parameters with respect
 *                          |        to which sensitivities are to be
 *                          |        computed.
 *                          | pbar: order of magnitude info. 
 *                          |       Typically, if p[plist[i]] is nonzero, 
 *                          |       pbar[i]=p[plist[i]].
 *                          | [p=NULL]
 *                          | [plist=NULL]
 *                          | [pbar=NULL]                               
 *                          |                                         
 * IDASetSensErrCon         | are sensitivity variables considered in 
 *                          | the error control?                      
 *                          | [TRUE]                                  
 *                          |                                         
 * IDASetSensMaxNonlinIters | Maximum number of nonlinear solver  
 *                          | iterations for sensitivity systems  
 *                          | (staggered)                         
 *                          | [4]                                 
 *                          |                                     
 * -------------------------------------------------------------- 
 * If successful, these functions return IDA_SUCCESS. If an argument  
 * has an illegal value, they return one of the error flags      
 * defined for the IDASet* routines.                              
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASetSensDQMethod(void *ida_mem, int DQtype, realtype DQrhomax);
SUNDIALS_EXPORT int IDASetSensParams(void *ida_mem, realtype *p, realtype *pbar, int *plist);
SUNDIALS_EXPORT int IDASetSensErrCon(void *ida_mem, booleantype errconS);
SUNDIALS_EXPORT int IDASetSensMaxNonlinIters(void *ida_mem, int maxcorS);
  
/*
 * ----------------------------------------------------------------
 * Function : IDASensInit                                       
 * ----------------------------------------------------------------
 * IDASensInit allocates and initializes memory related to      
 * sensitivity computations.                                      
 *                                                                
 * ida_mem is a pointer to IDAS memory returned by IDACreate.     
 *                                                                
 * Ns        is the number of sensitivities to be computed.       
 *                                                                
 * ism       is the type of corrector used in sensitivity         
 *           analysis. The legal values are: SIMULTANEOUS
 *           and STAGGERED (see previous description) 
 *                                                                
 * yS0       is the array of initial condition vectors for        
 *           sensitivity variables.                                
 *                                                                
 * ypS0      is the array of initial condition vectors for        
 *           sensitivity derivatives.                              
 *                                                                
 * If successful, IDASensInit returns SUCCESS. If an            
 * initialization error occurs, IDASensInit returns one of      
 * the error flags defined above.                                 
 *                                                                
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASensInit(void *ida_mem, int Ns, int ism, 
                                IDASensResFn resS,
                                N_Vector *yS0, N_Vector *ypS0);
  
  
/*
 * ----------------------------------------------------------------
 * Function : IDASensReInit                                       
 * ----------------------------------------------------------------
 * IDASensReInit re-initializes the IDAS sensitivity related      
 * memory for a problem, assuming it has already been allocated   
 * in prior calls to IDAInit and IDASensInit.                 
 *                                                                
 * All problem specification inputs are checked for errors.       
 * The number of sensitivities Ns is assumed to be unchanged      
 * since the previous call to IDASensInit.                      
 * If any error occurs during initialization, it is reported to   
 * the file whose file pointer is errfp.                          
 *                                                                
 * IDASensReInit potentially does some minimal memory allocation  
 * (for the sensitivity absolute tolerance).
 *                                                                
 * ----------------------------------------------------------------
 */
SUNDIALS_EXPORT int IDASensReInit(void *ida_mem, int ism, N_Vector *yS0, N_Vector *ypS0);

/*
 * -----------------------------------------------------------------
 * Function : IDASensToggleOff
 * -----------------------------------------------------------------
 * IDASensToggleOff deactivates sensitivity calculations.
 * It does NOT deallocate sensitivity-related memory so that 
 * sensitivity computations can be later toggled ON (through
 * IDASensReInit).
 * 
 * 
 * The return value is equal to IDA_SUCCESS = 0 if there were no
 * errors or IDA_MEM_NULL if ida_mem was NULL
 * -----------------------------------------------------------------
 */
SUNDIALS_EXPORT int IDASensToggleOff(void *ida_mem);

/*
 * -----------------------------------------------------------------
 * Functions : IDASensSStolerances
 *             IDASensSVtolerances
 *             IDASensEEtolerances
 * -----------------------------------------------------------------
 *
 * These functions specify the integration tolerances for sensitivity
 * variables. One of them MUST be called before the first call to IDASolve.
 *
 * IDASensSStolerances specifies scalar relative and absolute tolerances.
 * IDASensSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance for each sensitivity vector (a potentially different
 *   absolute tolerance for each vector component).
 * IDASensEEtolerances specifies that tolerances for sensitivity variables
 *   should be estimated from those provided for the state variables.
 *
 * The return value is equal to IDA_SUCCESS = 0 if there were no
 * errors; otherwise it is a negative int equal to:
 *   IDA_MEM_NULL  indicating ida_mem was NULL, or
 *   IDA_NO_SENS   indicating there was not a prior call to
 *                IDASensInit.
 *   IDA_ILL_INPUT indicating an input argument was illegal
 *                (e.g. negative tolerances)
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */
SUNDIALS_EXPORT int IDASensSStolerances(void *ida_mem, realtype reltolS, realtype *abstolS);
SUNDIALS_EXPORT int IDASensSVtolerances(void *ida_mem, realtype reltolS, N_Vector *abstolS);
SUNDIALS_EXPORT int IDASensEEtolerances(void *ida_mem);


/*
 * -----------------------------------------------------------------
 * Function : IDAQuadSensInit and IDAQuadSensReInit
 * -----------------------------------------------------------------
 * IDAQuadSensInit allocates and initializes memory related to
 * quadrature integration.
 *
 * IDAQuadSensReInit re-initializes IDAS' sensitivity quadrature 
 * related memory for a problem, assuming it has already been 
 * allocated in prior calls to IDAInit and IDAQuadSensInit.
 * The number of quadratures Ns is assumed to be unchanged
 * since the previous call to IDAQuadInit.
 *
 * ida_mem is a pointer to IDAS memory returned by IDACreate
 *
 * resQS     is the sensitivity righ-hand side function
 *        (pass NULL to use the internal DQ approximation)
 *
 * yQS    is an N_Vector with initial values for sensitivities

 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAQuadSensInit(void *ida_mem, IDAQuadSensRhsFn resQS, N_Vector *yQS0);
SUNDIALS_EXPORT int IDAQuadSensReInit(void *ida_mem, N_Vector *yQS0);

/*
 * -----------------------------------------------------------------
 * Functions : IDAQuadSensSStolerances
 *             IDAQuadSensSVtolerances
 *             IDAQuadSensEEtolerances
 * -----------------------------------------------------------------
 *
 * These functions specify the integration tolerances for quadrature
 * sensitivity variables. One of them MUST be called before the first
 * call to IDAS IF these variables are included in the error test.
 *
 * IDAQuadSensSStolerances specifies scalar relative and absolute tolerances.
 * IDAQuadSensSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance for each quadrature sensitivity vector (a potentially
 *   different absolute tolerance for each vector component).
 * IDAQuadSensEEtolerances specifies that tolerances for sensitivity variables
 *   should be estimated from those provided for the quadrature variables.
 *   In this case, tolerances for the quadrature variables must be
 *   specified through a call to one of IDAQuad**tolerances.
 *
 * The return value is equal to IDA_SUCCESS = 0 if there were no
 * errors; otherwise it is a negative int equal to:
 *   IDA_MEM_NULL     if ida_mem was NULL, or
 *   IDA_NO_QUADSENS  if there was not a prior call to
 *                    IDAQuadSensInit.
 *   IDA_ILL_INPUT    if an input argument was illegal
 *                   (e.g. negative tolerances)
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAQuadSensSStolerances(void *ida_mem, realtype reltolQS, realtype *abstolQS);
SUNDIALS_EXPORT int IDAQuadSensSVtolerances(void *ida_mem, realtype reltolQS, N_Vector *abstolQS);
SUNDIALS_EXPORT int IDAQuadSensEEtolerances(void *ida_mem);

/*
 * -----------------------------------------------------------------
 * Function: IDASetQuadSensErrCon 
 * -----------------------------------------------------------------
 * IDASetQuadSensErrCon specifies if quadrature sensitivity variables
 * are considered or not in the error control.
 *
 * If yes, tolerances for quadrature sensitivity variables are 
 * required. The function is optional, by default IDAS does not
 * quadrature sensitivities in error control.
 * 
 * The return value is equal to IDA_SUCCESS = 0 if there were no
 * errors or IDA_MEM_NULL if ida_mem was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASetQuadSensErrCon(void *ida_mem, booleantype errconQS);


/*
 * ----------------------------------------------------------------
 * Function : IDACalcIC                                          
 * ----------------------------------------------------------------
 * IDACalcIC calculates corrected initial conditions for the DAE  
 * system for a class of index-one problems of semi-implicit form.
 * It uses Newton iteration combined with a Linesearch algorithm. 
 * Calling IDACalcIC is optional. It is only necessary when the   
 * initial conditions do not solve the given system.  I.e., if    
 * y0 and yp0 are known to satisfy F(t0, y0, yp0) = 0, then       
 * a call to IDACalcIC is NOT necessary (for index-one problems). 
 *                                                                
 * A call to IDACalcIC must be preceded by a successful call to   
 * IDAInit or IDAReInit for the given DAE problem, and by a     
 * successful call to the linear system solver specification      
 * routine.                                                       
 *
 * The call to IDACalcIC should precede the call(s) to IDASolve   
 * for the given problem.                                         
 *                                                                
 * The arguments to IDACalcIC are as follows:                             
 *                                                                
 * ida_mem is the pointer to IDA memory returned by IDACreate.    
 *                                                                
 * icopt  is the option of IDACalcIC to be used.                  
 *        icopt = IDA_YA_YDP_INIT   directs IDACalcIC to compute 
 *                the algebraic components of y and differential  
 *                components of y', given the differential        
 *                components of y.  This option requires that the 
 *                N_Vector id was set through a call to IDASetId  
 *                specifying the differential and algebraic       
 *                components.                                     
 *        icopt = IDA_Y_INIT   directs IDACalcIC to compute all  
 *                components of y, given y'.  id is not required. 
 *                                                                
 * tout1  is the first value of t at which a soluton will be      
 *        requested (from IDASolve).  (This is needed here to     
 *        determine the direction of integration and rough scale  
 *        in the independent variable t.)                          
 *                                                                
 *                                                                
 * IDACalcIC returns an int flag.  Its symbolic values and their  
 * meanings are as follows.  (The numerical return values are set 
 * above in this file.)  All unsuccessful returns give a negative 
 * return value.  If IFACalcIC failed, y0 and yp0 contain         
 * (possibly) altered values, computed during the attempt.        
 *                                                                
 * IDA_SUCCESS         IDACalcIC was successful.  The corrected   
 *                     initial value vectors were stored internally.
 *                                                                
 * IDA_MEM_NULL        The argument ida_mem was NULL.             
 *                                                                
 * IDA_ILL_INPUT       One of the input arguments was illegal.    
 *                     See printed message.                       
 *                                                                
 * IDA_LINIT_FAIL      The linear solver's init routine failed.   
 *                                                                
 * IDA_BAD_EWT         Some component of the error weight vector  
 *                     is zero (illegal), either for the input    
 *                     value of y0 or a corrected value.          
 *                                                                
 * IDA_RES_FAIL        The user's residual routine returned 
 *                     a non-recoverable error flag.              
 *                                                                
 * IDA_FIRST_RES_FAIL  The user's residual routine returned 
 *                     a recoverable error flag on the first call,
 *                     but IDACalcIC was unable to recover.       
 *                                                                
 * IDA_LSETUP_FAIL     The linear solver's setup routine had a    
 *                     non-recoverable error.                     
 *                                                                
 * IDA_LSOLVE_FAIL     The linear solver's solve routine had a    
 *                     non-recoverable error.                     
 *                                                                
 * IDA_NO_RECOVERY     The user's residual routine, or the linear 
 *                     solver's setup or solve routine had a      
 *                     recoverable error, but IDACalcIC was       
 *                     unable to recover.                         
 *                                                                
 * IDA_CONSTR_FAIL     IDACalcIC was unable to find a solution    
 *                     satisfying the inequality constraints.     
 *                                                                
 * IDA_LINESEARCH_FAIL The Linesearch algorithm failed to find a  
 *                     solution with a step larger than steptol   
 *                     in weighted RMS norm.                      
 *                                                                
 * IDA_CONV_FAIL       IDACalcIC failed to get convergence of the 
 *                     Newton iterations.                         
 *                                                                
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDACalcIC(void *ida_mem, int icopt, realtype tout1); 

/*
 * ----------------------------------------------------------------
 * Function : IDASolve                                            
 * ----------------------------------------------------------------
 * IDASolve integrates the DAE over an interval in t, the         
 * independent variable. If itask is IDA_NORMAL, then the solver      
 * integrates from its current internal t value to a point at or  
 * beyond tout, then interpolates to t = tout and returns y(tret) 
 * in the user-allocated vector yret. In general, tret = tout.    
 * If itask is IDA_ONE_STEP, then the solver takes one internal
 * step of the independent variable and returns in yret the value
 * of y at the new internal independent variable value. In this
 * case, tout is used only during the first call to IDASolve to         
 * determine the direction of integration and the rough scale of  
 * the problem. If tstop is enabled (through a call to IDASetStopTime),
 * then IDASolve returns the solution at tstop. Once the integrator
 * returns at a tstop time, any future testing for tstop is disabled
 * (and can be reenabled only though a new call to IDASetStopTime).
 * The time reached by the solver is placed in (*tret). The
 * user is responsible for allocating the memory for this value.
 *                                                                
 * ida_mem is the pointer (void) to IDA memory returned by        
 *         IDACreate.
 *                                                                
 * tout    is the next independent variable value at which a      
 *         computed solution is desired.                          
 *                                                                
 * tret    is a pointer to a real location.  IDASolve sets (*tret)
 *         to the actual t value reached, corresponding to the
 *         solution vector yret.  In IDA_NORMAL mode, with no
 *         errors and no roots found, (*tret) = tout.
 *
 * yret    is the computed solution vector.  With no errors,
 *         yret = y(tret).                                        
 *                                                                
 * ypret   is the derivative of the computed solution at t = tret.
 *                                                                
 * Note: yret and ypret may be the same N_Vectors as y0 and yp0   
 * in the call to IDAInit or IDAReInit.                         
 *                                                                
 * itask   is IDA_NORMAL or IDA_ONE_STEP. These two modes are described above.
 *
 *
 * The return values for IDASolve are described below.            
 * (The numerical return values are defined above in this file.)  
 * All unsuccessful returns give a negative return value.         
 *                                                                
 * IDA_SUCCESS
 *   IDASolve succeeded and no roots were found.                       
 *
 * IDA_ROOT_RETURN:  IDASolve succeeded, and found one or more roots.
 *   If nrtfn > 1, call IDAGetRootInfo to see which g_i were found
 *   to have a root at (*tret).
 *
 * IDA_TSTOP_RETURN: 
 *   IDASolve returns computed results for the independent variable 
 *   value tstop. That is, tstop was reached.                            
 *                                                                
 * IDA_MEM_NULL: 
 *   The ida_mem argument was NULL.            
 *                                                                
 * IDA_ILL_INPUT: 
 *   One of the inputs to IDASolve is illegal. This includes the 
 *   situation when a component of the error weight vectors 
 *   becomes < 0 during internal stepping.  It also includes the
 *   situation where a root of one of the root functions was found
 *   both at t0 and very near t0.  The ILL_INPUT flag          
 *   will also be returned if the linear solver function IDA---
 *   (called by the user after calling IDACreate) failed to set one 
 *   of the linear solver-related fields in ida_mem or if the linear 
 *   solver's init routine failed. In any case, the user should see 
 *   the printed error message for more details.                
 *                                                                
 * IDA_TOO_MUCH_WORK: 
 *   The solver took mxstep internal steps but could not reach tout. 
 *   The default value for mxstep is MXSTEP_DEFAULT = 500.                
 *                                                                
 * IDA_TOO_MUCH_ACC: 
 *   The solver could not satisfy the accuracy demanded by the user 
 *   for some internal step.   
 *                                                                
 * IDA_ERR_FAIL:
 *   Error test failures occurred too many times (=MXETF = 10) during 
 *   one internal step.  
 *                                                                
 * IDA_CONV_FAIL: 
 *   Convergence test failures occurred too many times (= MXNCF = 10) 
 *   during one internal step.                                          
 *                                                                
 * IDA_LSETUP_FAIL: 
 *   The linear solver's setup routine failed  
 *   in an unrecoverable manner.                    
 *                                                                
 * IDA_LSOLVE_FAIL: 
 *   The linear solver's solve routine failed  
 *   in an unrecoverable manner.                    
 *                                                                
 * IDA_CONSTR_FAIL:
 *    The inequality constraints were violated, 
 *    and the solver was unable to recover.         
 *                                                                
 * IDA_REP_RES_ERR: 
 *    The user's residual function repeatedly returned a recoverable 
 *    error flag, but the solver was unable to recover.                 
 *                                                                
 * IDA_RES_FAIL:
 *    The user's residual function returned a nonrecoverable error 
 *    flag.
 *                                                                
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASolve(void *ida_mem, realtype tout, realtype *tret,
			     N_Vector yret, N_Vector ypret, int itask);

/*
 * ----------------------------------------------------------------
 * Function: IDAGetDky                                     
 * ----------------------------------------------------------------
 *
 * This routine computes the k-th derivative of the interpolating
 * polynomial at the time t and stores the result in the vector dky.
 *
 * The return values are:                                         
 *   IDA_SUCCESS:  succeess.                                  
 *   IDA_BAD_T:    t is not in the interval [tn-hu,tn].                   
 *   IDA_MEM_NULL: The ida_mem argument was NULL.
 *   IDA_BAD_DKY  if the dky vector is NULL.
 *   IDA_BAD_K    if the requested k is not in the range 0,1,...,order used 
 *                                                                
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetDky(void *ida_mem, realtype t, int k, N_Vector dky);

/* ----------------------------------------------------------------
 * Integrator optional output extraction functions                
 * ----------------------------------------------------------------
 *                                                                
 * The following functions can be called to get optional outputs  
 * and statistics related to the main integrator.                 
 * ---------------------------------------------------------------- 
 *                                                                
 * IDAGetWorkSpace returns the IDA real and integer workspace sizes      
 * IDAGetNumSteps returns the cumulative number of internal       
 *       steps taken by the solver                                
 * IDAGetNumResEvals returns the number of calls to the user's    
 *       res function                                             
 * IDAGetNumLinSolvSetups returns the number of calls made to     
 *       the linear solver's setup routine                        
 * IDAGetNumErrTestFails returns the number of local error test   
 *       failures that have occured                               
 * IDAGetNumBacktrackOps returns the number of backtrack          
 *       operations done in the linesearch algorithm in IDACalcIC 
 * IDAGetConsistentIC returns the consistent initial conditions
 *       computed by IDACalcIC
 * IDAGetLastOrder returns the order used during the last         
 *       internal step                                            
 * IDAGetCurentOrder returns the order to be used on the next     
 *       internal step                                            
 * IDAGetActualInitStep returns the actual initial step size      
 *       used by IDA                                              
 * IDAGetLastStep returns the step size for the last internal     
 *       step (if from IDASolve), or the last value of the        
 *       artificial step size h (if from IDACalcIC)               
 * IDAGetCurrentStep returns the step size to be attempted on the 
 *       next internal step                                       
 * IDAGetCurrentTime returns the current internal time reached    
 *       by the solver                                            
 * IDAGetTolScaleFactor returns a suggested factor by which the   
 *       user's tolerances should be scaled when too much         
 *       accuracy has been requested for some internal step       
 * IDAGetErrWeights returns the current state error weight vector.        
 *       The user must allocate space for eweight.
 * IDAGetEstLocalErrors returns the estimated local errors. The user
 *       must allocate space for the vector ele.
 * IDAGetNumGEvals returns the number of calls to the user's
 *       g function (for rootfinding)
 * IDAGetRootInfo returns the indices for which g_i was found to 
 *       have a root. The user must allocate space for rootsfound.
 *       For i = 0 ... nrtfn-1, rootsfound[i] = 1 if g_i has a root,
 *       and rootsfound[i]= 0 if not.
 *                                                                
 * IDAGet* return values:
 *   IDA_SUCCESS   if succesful
 *   IDA_MEM_NULL  if the IDAS memory was NULL
 *   IDA_ILL_INPUT if some input is illegal
 *
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetWorkSpace(void *ida_mem, long int *lenrw, long int *leniw);
SUNDIALS_EXPORT int IDAGetNumSteps(void *ida_mem, long int *nsteps);
SUNDIALS_EXPORT int IDAGetNumResEvals(void *ida_mem, long int *nrevals);
SUNDIALS_EXPORT int IDAGetNumLinSolvSetups(void *ida_mem, long int *nlinsetups);
SUNDIALS_EXPORT int IDAGetNumErrTestFails(void *ida_mem, long int *netfails);
SUNDIALS_EXPORT int IDAGetNumBacktrackOps(void *ida_mem, long int *nbacktr);
SUNDIALS_EXPORT int IDAGetConsistentIC(void *ida_mem, N_Vector yy0_mod, N_Vector yp0_mod);
SUNDIALS_EXPORT int IDAGetLastOrder(void *ida_mem, int *klast);
SUNDIALS_EXPORT int IDAGetCurrentOrder(void *ida_mem, int *kcur);
SUNDIALS_EXPORT int IDAGetActualInitStep(void *ida_mem, realtype *hinused);
SUNDIALS_EXPORT int IDAGetLastStep(void *ida_mem, realtype *hlast);
SUNDIALS_EXPORT int IDAGetCurrentStep(void *ida_mem, realtype *hcur);
SUNDIALS_EXPORT int IDAGetCurrentTime(void *ida_mem, realtype *tcur);
SUNDIALS_EXPORT int IDAGetTolScaleFactor(void *ida_mem, realtype *tolsfact);
SUNDIALS_EXPORT int IDAGetErrWeights(void *ida_mem, N_Vector eweight);
SUNDIALS_EXPORT int IDAGetEstLocalErrors(void *ida_mem, N_Vector ele);
SUNDIALS_EXPORT int IDAGetNumGEvals(void *ida_mem, long int *ngevals);
SUNDIALS_EXPORT int IDAGetRootInfo(void *ida_mem, int *rootsfound);

/*
 * ----------------------------------------------------------------
 * As a convenience, the following function provides the          
 * optional outputs in a group.                                   
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetIntegratorStats(void *ida_mem, long int *nsteps, 
					  long int *nrevals, long int *nlinsetups, 
					  long int *netfails, int *qlast, int *qcur, 
					  realtype *hinused, realtype *hlast, realtype *hcur, 
					  realtype *tcur);
/*
 * ----------------------------------------------------------------
 * Nonlinear solver optional output extraction functions          
 * ----------------------------------------------------------------
 *                                                                
 * The following functions can be called to get optional outputs  
 * and statistics related to the nonlinear solver.                
 * -------------------------------------------------------------- 
 *                                                                
 * IDAGetNumNonlinSolvIters returns the number of nonlinear       
 *       solver iterations performed.                             
 * IDAGetNumNonlinSolvConvFails returns the number of nonlinear   
 *       convergence failures.                                    
 *                                                                
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters);
SUNDIALS_EXPORT int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails);

/*
 * ----------------------------------------------------------------
 * As a convenience, the following function provides the          
 * nonlinear solver optional outputs in a group.                                   
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetNonlinSolvStats(void *ida_mem, long int *nniters, 
					  long int *nncfails);

/*
 * -----------------------------------------------------------------
 * Quadrature integration solution extraction routines
 * -----------------------------------------------------------------
 * The following function can be called to obtain the quadrature
 * variables after a successful integration step.
 * If quadratures were not computed, it returns IDA_NO_QUAD.
 *
 * IDAGetQuad returns the quadrature variables at the same time
 *   as that at which IDASolve returned the solution.
 *
 * IDAGetQuadDky returns the quadrature variables (or their 
 *   derivatives up to the current method order) at any time within
 *   the last integration step (dense output). 
 *
 * The output vectors yQout and dky must be allocated by the user.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetQuad(void *ida_mem, realtype *t, N_Vector yQout);
SUNDIALS_EXPORT int IDAGetQuadDky(void *ida_mem, realtype t, int k, N_Vector dky);
/*
 * -----------------------------------------------------------------
 * Quadrature integration optional output extraction routines
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the integration of quadratures.
 * -----------------------------------------------------------------
 * IDAGetQuadNumRhsEvals returns the number of calls to the
 *                       user function rhsQ defining the right hand
 *                       side of the quadrature variables.
 * IDAGetQuadNumErrTestFails returns the number of local error
 *                           test failures for quadrature variables.
 * IDAGetQuadErrWeights returns the vector of error weights for
 *                      the quadrature variables. The user must
 *                      allocate space for ewtQ.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetQuadNumRhsEvals(void *ida_mem, long int *nrhsQevals);
SUNDIALS_EXPORT int IDAGetQuadNumErrTestFails(void *ida_mem, long int *nQetfails);
SUNDIALS_EXPORT int IDAGetQuadErrWeights(void *ida_mem, N_Vector eQweight);

/*
 * -----------------------------------------------------------------
 * As a convenience, the following function provides the
 * optional outputs in a group.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetQuadStats(void *ida_mem, 
                                    long int *nrhsQevals, long int *nQetfails);

/*
 * -----------------------------------------------------------------
 * Sensitivity solution extraction routines
 * -----------------------------------------------------------------
 * The following functions can be called to obtain the sensitivity
 * variables after a successful integration step.
 * 
 * IDAGetSens and IDAGetSens1 return all the sensitivity vectors
 *   or only one of them, respectively, at the same time as that at 
 *   which IDASolve returned the solution.
 *   The array of output vectors or output vector ySout must be
 *   allocated by the user.
 *
 * IDAGetSensDky1 computes the kth derivative of the is-th
 *   sensitivity (is=1, 2, ..., Ns) of the y function at time t,
 *   where tn-hu <= t <= tn, tn denotes the current internal time
 *   reached, and hu is the last internal step size successfully
 *   used by the solver. The user may request k=0, 1, ..., qu,
 *   where qu is the current order.
 *   The is-th sensitivity derivative vector is returned in dky.
 *   This vector must be allocated by the caller. It is only legal
 *   to call this function after a successful return from IDASolve
 *   with sensitivity computations enabled.
 *   Arguments have the same meaning as in IDADGetky.
 *
 * IDAGetSensDky computes the k-th derivative of all
 *   sensitivities of the y function at time t. It repeatedly calls
 *   IDAGetSensDky. The argument dkyS must be a pointer to
 *   N_Vector and must be allocated by the user to hold at least Ns
 *   vectors.
 *
 * Return values are similar to those of IDAGetDky. Additionally,
 * these functions can return IDA_NO_SENS if sensitivities were
 * not computed and IDA_BAD_IS if is < 0 or is >= Ns.
 * -----------------------------------------------------------------
 */


SUNDIALS_EXPORT int IDAGetSens(void *ida_mem, realtype *tret, N_Vector *yySout);
SUNDIALS_EXPORT int IDAGetSens1(void *ida_mem, realtype *tret, int is, N_Vector yySret);

SUNDIALS_EXPORT int IDAGetSensDky(void *ida_mem, realtype t, int k, N_Vector *dkyS);
SUNDIALS_EXPORT int IDAGetSensDky1(void *ida_mem, realtype t, int k, int is, N_Vector dkyS);


/*
 * -----------------------------------------------------------------
 * Consistent sensitivity IC calculation optional outputs
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetSensConsistentIC(void *ida_mem, N_Vector *yyS0, N_Vector *ypS0);

/*
 * -----------------------------------------------------------------
 * Forward sensitivity optional output extraction routines
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the integration of sensitivities.
 * -----------------------------------------------------------------
 * IDAGetSensNumResEvals returns the number of calls to the
 *                       sensitivity residual function.
 * IDAGetNumResEvalsSens returns the number of calls to the
 *                       user res routine due to finite difference
 *                       evaluations of the sensitivity equations.
 * IDAGetSensNumErrTestFails returns the number of local error
 *                           test failures for sensitivity variables.
 * IDAGetSensNumLinSolvSetups returns the number of calls made
 *                            to the linear solver's setup routine
 *                            due to sensitivity computations.
 * IDAGetSensErrWeights returns the sensitivity error weight
 *                      vectors. The user need not allocate space
 *                      for ewtS.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetSensNumResEvals(void *ida_mem, long int *nresSevals);
SUNDIALS_EXPORT int IDAGetNumResEvalsSens(void *ida_mem, long int *nresevalsS);
SUNDIALS_EXPORT int IDAGetSensNumErrTestFails(void *ida_mem, long int *nSetfails);
SUNDIALS_EXPORT int IDAGetSensNumLinSolvSetups(void *ida_mem, long int *nlinsetupsS);
SUNDIALS_EXPORT int IDAGetSensErrWeights(void *ida_mem, N_Vector_S eSweight);

/*
 * -----------------------------------------------------------------
 * As a convenience, the following function provides the
 * optional outputs in a group.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetSensStats(void *ida_mem, long int *nresSevals, 
                                    long int *nresevalsS, 
                                    long int *nSetfails, 
                                    long int *nlinsetupsS);

/*
 * ----------------------------------------------------------------
 * Sensitivity nonlinear solver optional output extraction functions          
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetSensNumNonlinSolvIters(void *ida_mem, long int *nSniters);
SUNDIALS_EXPORT int IDAGetSensNumNonlinSolvConvFails(void *ida_mem, 
                                                     long int *nSncfails);
SUNDIALS_EXPORT int IDAGetSensNonlinSolvStats(void *ida_mem, 
                                              long int *nSniters, 
                                              long int *nSncfails);


/*
 * -----------------------------------------------------------------
 * Quadrature sensitivity optional output extraction routines
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs and
 * statistics related to the integration of quadrature sensitivitiess.
 * -----------------------------------------------------------------
 * IDAGetQuadSensNumRhsEvals returns the number of calls to the
 *       user function fQS defining the right hand side of the 
 *       quadrature sensitivity equations.
 * IDAGetQuadSensNumErrTestFails returns the number of local error
 *       test failures for quadrature sensitivity variables.
 * IDAGetQuadSensErrWeights returns the vector of error weights
 *       for the quadrature sensitivity variables. The user must
 *       allocate space for ewtQS.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetQuadSensNumRhsEvals(void *ida_mem, long int *nrhsQSevals);
SUNDIALS_EXPORT int IDAGetQuadSensNumErrTestFails(void *ida_mem, long int *nQSetfails);
SUNDIALS_EXPORT int IDAGetQuadSensErrWeights(void *ida_mem, N_Vector *eQSweight);

/*
 * -----------------------------------------------------------------
 * As a convenience, the following function provides the above
 * optional outputs in a group.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetQuadSensStats(void *ida_mem,
                                          long int *nrhsQSevals,
                                          long int *nQSetfails);


/*
 * -----------------------------------------------------------------
 * Quadrature Sensitivity solution extraction routines
 * -----------------------------------------------------------------
 * The following functions can be called to obtain the sensitivity
 * variables after a successful integration step.
 * 
 * IDAGetQuadSens and IDAGetQuadSens1 return all the sensitivity 
 *   vectors or only one of them, respectively, at the same time 
 *   as that at which IDASolve returned the solution.
 *   The array of output vectors or output vector yQSout must be
 *   allocated by the user.
 *
 * IDAGetQuadSensDky1 computes the kth derivative of the is-th
 *   sensitivity (is=1, 2, ..., Ns) of the quadrature function at 
 *   time t, where  tn - hu <= t <= tn,  tn denotes the current 
 *   internal  time reached  and hu is the  last internal  
 *   successfully step size. The user may request k=0,..., qu, 
 *   where qu is the current order.
 *
 *   The is-th sensitivity derivative vector is returned in dky.
 *   This vector must be allocated by the caller. It is only legal
 *   to call this function after a successful return from IDASolve
 *   with sensitivity computations enabled.
 *   Arguments have the same meaning as in IDADGetky.
 *
 * IDAGetQuadSensDky computes the k-th derivative of all
 *   sensitivities of the y function at time t. It repeatedly calls
 *   IDAGetQuadSensDky. The argument dkyS must be a pointer to
 *   N_Vector and must be allocated by the user to hold at least Ns
 *   vectors.
 *
 * Return values are similar to those of IDAGetDky. Additionally,
 * these functions can return IDA_NO_SENS if sensitivities were
 * not computed and IDA_BAD_IS if is < 0 or is >= Ns.
 * -----------------------------------------------------------------
 */


SUNDIALS_EXPORT int IDAGetQuadSens(void *ida_mem, realtype *tret, N_Vector *yyQSout);
SUNDIALS_EXPORT int IDAGetQuadSens1(void *ida_mem, realtype *tret, int is, N_Vector yyQSret);

SUNDIALS_EXPORT int IDAGetQuadSensDky(void *ida_mem, realtype t, int k, N_Vector *dkyQS);
SUNDIALS_EXPORT int IDAGetQuadSensDky1(void *ida_mem, realtype t, int k, int is, N_Vector dkyQS);


/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with an IDAS return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *IDAGetReturnFlagName(long int flag);

/*
 * ----------------------------------------------------------------
 * Function : IDAFree                                             
 * ----------------------------------------------------------------
 * IDAFree frees the problem memory IDA_mem allocated by          
 * IDAInit.  Its only argument is the pointer idamem            
 * returned by IDAInit.                                         
 * ----------------------------------------------------------------
 */

SUNDIALS_EXPORT void IDAFree(void **ida_mem);

/*
 * -----------------------------------------------------------------
 * Function : IDAQuadFree
 * -----------------------------------------------------------------
 * IDAQuadFree frees the problem memory in ida_mem allocated
 * for quadrature integration. Its only argument is the pointer
 * ida_mem returned by IDACreate.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void IDAQuadFree(void *ida_mem);

/*
 * -----------------------------------------------------------------
 * Function : IDASensFree
 * -----------------------------------------------------------------
 * IDASensFree frees the problem memory in ida_mem allocated
 * for sensitivity analysis. Its only argument is the pointer
 * ida_mem returned by IDACreate.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void IDASensFree(void *ida_mem);

/*
 * -----------------------------------------------------------------
 * Function : IDAQuadSensFree
 * -----------------------------------------------------------------
 * IDAQuadSensFree frees the problem memory in ida_mem allocated
 * for quadrature sensitivity analysis. Its only argument is the 
 * pointer ida_mem returned by IDACreate.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void IDAQuadSensFree(void* ida_mem);

/* 
 * =================================================================
 *
 * INITIALIZATION AND DEALLOCATION FUNCTIONS FOR BACKWARD PROBLEMS 
 *
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * IDAAdjInit
 * -----------------------------------------------------------------
 * IDAAdjInit specifies some parameters for ASA, initializes ASA
 * and allocates space for the adjoint memory structure.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAAdjInit(void *ida_mem, long int steps, int interp);
/*
 * -----------------------------------------------------------------
 * IDAAdjReInit
 * -----------------------------------------------------------------
 * IDAAdjReInit reinitializes the IDAS memory structure for ASA,
 * assuming that the number of steps between check points and the
 * type of interpolation remained unchanged. The list of check points
 * (and associated memory) is deleted. The list of backward problems
 * is kept (however, new backward problems can be added to this list
 * by calling IDACreateB). The IDAS memory for the forward and 
 * backward problems can be reinitialized separately by calling 
 * IDAReInit and IDAReInitB, respectively.
 * NOTE: if a entirely new list of backward problems is desired,
 *   then simply free the adjoint memory (by calling IDAAdjFree)
 *   and reinitialize ASA with IDAAdjReInit 
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAAdjReInit(void *ida_mem);

/*
 * -----------------------------------------------------------------
 * IDAAdjFree
 * -----------------------------------------------------------------
 * IDAAdjFree frees the memory allocated by IDAAdjInit.
 * It is typically called by IDAFree.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void IDAAdjFree(void *ida_mem);

/*
 * =================================================================
 *
 * OPTIONAL INPUT FUNCTIONS FOR BACKWARD PROBLEMS
 *
 * =================================================================
 */


/*
 * =================================================================
 *
 * Interfaces to IDAS functions for setting-up backward problems.
 *
 * =================================================================
 */

SUNDIALS_EXPORT int IDACreateB(void *ida_mem, int *which);

SUNDIALS_EXPORT int IDAInitB(void *ida_mem, int which, IDAResFnB resB,
                             realtype tB0, N_Vector yyB0, N_Vector ypB0);

SUNDIALS_EXPORT int IDAInitBS(void *ida_mem, int which, IDAResFnBS resS,
                              realtype tB0, N_Vector yyB0, N_Vector ypB0);

SUNDIALS_EXPORT int IDAReInitB(void *ida_mem, int which,
			       realtype tB0, N_Vector yyB0, N_Vector ypB0);

SUNDIALS_EXPORT int IDASStolerancesB(void *ida_mem, int which, 
                                     realtype relTolB, realtype absTolB);
SUNDIALS_EXPORT int IDASVtolerancesB(void *ida_mem, int which, 
                                     realtype relTolB, N_Vector absTolB);

SUNDIALS_EXPORT int IDAQuadInitB(void *ida_mem, int which, 
                                 IDAQuadRhsFnB rhsQB, N_Vector yQB0);

SUNDIALS_EXPORT int IDAQuadInitBS(void *ida_mem, int which, 
                                  IDAQuadRhsFnBS rhsQS, N_Vector yQB0);

SUNDIALS_EXPORT int IDAQuadReInitB(void *ida_mem, int which, N_Vector yQB0);

SUNDIALS_EXPORT int IDAQuadSStolerancesB(void *ida_mem, int which,
                                         realtype reltolQB, realtype abstolQB);
SUNDIALS_EXPORT int IDAQuadSVtolerancesB(void *ida_mem, int which,
                                         realtype reltolQB, N_Vector abstolQB);

/*
 * ----------------------------------------------------------------
 * The following functions computes consistent initial conditions
 * for the backward problems.
 * ----------------------------------------------------------------
 * Function : IDACalcICB                                         
 * ----------------------------------------------------------------
 * IDACalcICB calculates corrected initial conditions for a DAE  
 * backward system (index-one in semi-implicit form).
 * ----------------------------------------------------------------
 * Function : IDACalcICBS
 * ----------------------------------------------------------------
 * IDACalcICBS calculates corrected initial conditions for a DAE  
 * backward problems that also depends on  the sensitivities.
 *
 * They use Newton iteration combined with a Linesearch algorithm. 
 *
 * Calling IDACalcICB(S) is optional. It is only necessary when the   
 * initial conditions do not solve the given system.  I.e., if    
 * yB0 and ypB0 are known to satisfy the backward problem, then       
 * a call to IDACalcICB is NOT necessary (for index-one problems). 
 *
 * Any call to IDACalcICB(S) should precede the call(s) to 
 * IDASolveB for the given problem. 
 * 
 * The functions compute the algebraic components of y and 
 * differential components of y', given the differential        
 * components of y.  This option requires that the N_Vector id was 
 * set through a call to IDASetIdB specifying the differential and 
 * algebraic components.      
 *
 * The arguments to IDACalcICB(S) are as follows:
 *                                                                
 * ida_mem  is the pointer to IDA memory returned by IDACreate.                        
 *
 * which    is the index of the backward problem returned by
 *          IDACreateB
 *
 * tout1    is the first value of t at which a soluton will be      
 *          requested (from IDASolveB).  (This is needed here to     
 *          determine the direction of integration and rough 
 *          scale  in the independent variable t.)  
 *
 * yy0      state variables y and y' corresponding to the initial
 * yp0      time at which the backward problem is (re)started.
 *
 * yyS0     sensitivities variables corresponding to the initial
 * ypS0     time at which the backward problem is (re)started. 
 *
 * Return value is a int flag. For more information see IDACalcIC.
*/

SUNDIALS_EXPORT int IDACalcICB (void *ida_mem, int which, realtype tout1,
                                N_Vector yy0, N_Vector yp0);

SUNDIALS_EXPORT int IDACalcICBS(void *ida_mem, int which, realtype tout1,
                                N_Vector yy0, N_Vector yp0,
                                N_Vector *yyS0, N_Vector *ypS0); 


/*
 * =================================================================
 *
 * MAIN SOLVER FUNCTIONS FOR FORWARD PROBLEMS
 *
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * IDASolveF
 * -----------------------------------------------------------------
 * IDASolveF integrates towards tout and returns solution into yret
 * and ypret.
 *
 * In the same time, it stores check point data every 'steps'.
 *
 * IDASolveF can be called repeatedly by the user.
 *
 * ncheckPtr represents the number of check points stored so far.
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASolveF(void *ida_mem, realtype tout, 
                              realtype *tret, 
                              N_Vector yret, N_Vector ypret, 
                              int itask, int *ncheckPtr);


/*
 * -----------------------------------------------------------------
 * IDASolveB
 * -----------------------------------------------------------------
 * IDASolveB performs the integration of all backward problems 
 * specified through calls to IDACreateB through a sequence of 
 * forward-backward runs in between consecutive check points. It can 
 * be called either in IDA_NORMAL or IDA_ONE_STEP mode. After a 
 * successful return from IDASolveB, the solution and quadrature 
 * variables at the current return time for any given backward 
 * problem can be obtained by calling IDAGetB and IDAGetQuadB.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASolveB(void *ida_mem, realtype tBout, int itaskB);


/*
 * =================================================================
 *
 * OPTIONAL INPUT FUNCTIONS FOR BACKWARD PROBLEMS
 *
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * IDAAdjSetNoSensi
 * -----------------------------------------------------------------
 * Disables the forward sensitivity analysis in IDASolveF.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAAdjSetNoSensi(void *ida_mem);

/*
 * -----------------------------------------------------------------
 * Optional input functions for backward problems
 * -----------------------------------------------------------------
 * These functions are just wrappers around the corresponding
 * functions from the forward module, with some particularizations 
 * for the backward integration.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDASetUserDataB(void *ida_mem, int which, void *user_dataB);

SUNDIALS_EXPORT int IDASetMaxOrdB(void *ida_mem, int which, int maxordB);

SUNDIALS_EXPORT int IDASetMaxNumStepsB(void *ida_mem, int which, long int mxstepsB);

SUNDIALS_EXPORT int IDASetInitStepB(void *ida_mem, int which, realtype hinB);

SUNDIALS_EXPORT int IDASetMaxStepB(void *ida_mem, int which, realtype hmaxB);

SUNDIALS_EXPORT int IDASetSuppressAlgB(void *ida_mem, int which, 
                                       booleantype suppressalgB);
SUNDIALS_EXPORT int IDASetIdB(void *ida_mem, int which, N_Vector idB);

SUNDIALS_EXPORT int IDASetConstraintsB(void *ida_mem, int which, 
                                       N_Vector constraintsB);

SUNDIALS_EXPORT int IDASetQuadErrConB(void *ida_mem, int which, int errconQB);

/*
 * =================================================================
 *
 * EXTRACTION AND DENSE OUTPUT FUNCTIONS FOR BACKWARD PROBLEMS
 *
 * =================================================================
 */
  
/*
 * -----------------------------------------------------------------
 * IDAGetB and IDAGetQuadB
 * -----------------------------------------------------------------
 * Extraction functions for the solution and quadratures for a given 
 * backward problem. They return their corresponding output vector
 * at the current time reached by the integration of the backward
 * problem. To obtain the solution or quadratures associated with
 * a given backward problem at some other time within the last 
 * integration step (dense output), first obtain a pointer to the
 * proper IDAS memory by calling IDAGetAdjIDABmem and then use it
 * to call IDAGetDky and IDAGetQuadDky.  
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetB(void* ida_mem, int which, realtype *tret,
                            N_Vector yy, N_Vector yp);

SUNDIALS_EXPORT int IDAGetQuadB(void *ida_mem, int which, 
                                realtype *tret, N_Vector qB);

/*
 * =================================================================
 *
 * OPTIONAL OUTPUT FUNCTIONS FOR BACKWARD PROBLEMS
 *
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * IDAGetAdjIDABmem
 * -----------------------------------------------------------------
 * IDAGetAdjIDABmem returns a (void *) pointer to the IDAS
 * memory allocated for the backward problem. This pointer can
 * then be used to call any of the IDAGet* IDAS routines to
 * extract optional output for the backward integration phase.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void *IDAGetAdjIDABmem(void *ida_mem, int which);

/* 
 * -----------------------------------------------------------------
 * IDAGetConsistentICB
 * -----------------------------------------------------------------
 * IDAGetConsistentIC returns the consistent initial conditions
 * computed by IDACalcICB or IDCalcICBS
 */
SUNDIALS_EXPORT int IDAGetConsistentICB(void *ida_mem, int which, 
                                        N_Vector yyB0, N_Vector ypB0);

/*
 * -----------------------------------------------------------------
 * IDAGetAdjY
 * -----------------------------------------------------------------
 * Returns the interpolated forward solution at time t. This
 * function is a wrapper around the interpType-dependent internal
 * function.
 * The calling function must allocate space for yy and yp.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetAdjY(void *ida_mem, realtype t, 
                               N_Vector yy, N_Vector yp);


/*
 * -----------------------------------------------------------------
 * IDAGetAdjCheckPointsInfo
 * -----------------------------------------------------------------
 * Loads an array of nckpnts structures of type IDAadjCheckPointRec
 * defined below.
 *
 * The user must allocate space for ckpnt (ncheck+1).
 * -----------------------------------------------------------------
 */

  typedef struct {
  void *my_addr;
  void *next_addr;
  realtype t0;
  realtype t1;
  long int nstep;
  int order;
  realtype step;
  } IDAadjCheckPointRec;

SUNDIALS_EXPORT int IDAGetAdjCheckPointsInfo(void *ida_mem, 
                                             IDAadjCheckPointRec *ckpnt);

/*
 * -----------------------------------------------------------------
 * IDAGetAdjDataPointHermite
 * -----------------------------------------------------------------
 * Returns the 2 vectors stored for cubic Hermite interpolation at
 * the data point 'which'. The user must allocate space for yy and
 * yd. 
 *
 * Returns IDA_MEM_NULL if ida_mem is NULL, IDA_ILL_INPUT if the 
 * interpolation type previously specified is not IDA_HERMITE or
 * IDA_SUCCESS otherwise.
 *
 * -----------------------------------------------------------------
 * IDAGetAdjDataPointPolynomial
 * -----------------------------------------------------------------
 * Returns the vector stored for polynomial interpolation at the
 * data point 'which'. The user must allocate space for y.
 *
 * Returns IDA_MEM_NULL if ida_mem is NULL, IDA_ILL_INPUT if the 
 * interpolation type previously specified is not IDA_POLYNOMIAL or
 * IDA_SUCCESS otherwise.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetAdjDataPointHermite(void *ida_mem, int which,
                                              realtype *t, N_Vector yy, N_Vector yd);

SUNDIALS_EXPORT int IDAGetAdjDataPointPolynomial(void *ida_mem, int which,
                                                 realtype *t, int *order, 
                                                 N_Vector y);
  
  
/*
 * -----------------------------------------------------------------
 * IDAGetAdjCurrentCheckPoint
 * -----------------------------------------------------------------
 * Returns the address of the 'active' check point.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int IDAGetAdjCurrentCheckPoint(void *ida_mem, void **addr);


#ifdef __cplusplus
}
#endif

#endif
