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
 * Common header file for the direct linear solvers in CVODES.
 *
 * Part I contains type definitions and function prototypes for 
 * using a CVDLS linear solver on forward problems (IVP 
 * integration and/or FSA)
 *
 * Part II contains type definitions and function prototypes for 
 * using a CVDLS linear solver on adjoint (backward) problems
 * -----------------------------------------------------------------
 */

#ifndef _CVSDLS_H
#define _CVSDLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 * C V S D I R E C T     C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * CVSDIRECT return values 
 * -----------------------------------------------------------------
 */

#define CVDLS_SUCCESS           0
#define CVDLS_MEM_NULL         -1
#define CVDLS_LMEM_NULL        -2
#define CVDLS_ILL_INPUT        -3
#define CVDLS_MEM_FAIL         -4

/* Additional last_flag values */

#define CVDLS_JACFUNC_UNRECVR  -5
#define CVDLS_JACFUNC_RECVR    -6

/* Return values for the adjoint module */

#define CVDLS_NO_ADJ           -101
#define CVDLS_LMEMB_NULL       -102

/*
 * =================================================================
 * PART I:  F O R W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * FUNCTION TYPES
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Type: CVDlsDenseJacFn
 * -----------------------------------------------------------------
 *
 * A dense Jacobian approximation function Jac must be of type 
 * CVDlsDenseJacFn. Its parameters are:
 *
 * N   is the problem size.
 *
 * Jac is the dense matrix (of type DlsMat) that will be loaded
 *     by a CVDlsDenseJacFn with an approximation to the Jacobian 
 *     matrix J = (df_i/dy_j) at the point (t,y). 
 *
 * t   is the current value of the independent variable.
 *
 * y   is the current value of the dependent variable vector,
 *     namely the predicted value of y(t).
 *
 * fy  is the vector f(t,y).
 *
 * user_data is a pointer to user data - the same as the user_data
 *     parameter passed to CVodeSetFdata.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for
 * vectors of length N which can be used by a CVDlsDenseJacFn
 * as temporary storage or work space.
 *
 * A CVDlsDenseJacFn should return 0 if successful, a positive 
 * value if a recoverable error occurred, and a negative value if 
 * an unrecoverable error occurred.
 *
 * -----------------------------------------------------------------
 *
 * NOTE: The following are two efficient ways to load a dense Jac:         
 * (1) (with macros - no explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = DENSE_COL(Jac,j);                                 
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * (2) (without macros - explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = (Jac->data)[j];                                   
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
 * efficient in general.  It is only appropriate for use in small 
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively 
 *     (see cvode.h). The unit roundoff is available as 
 *     UNIT_ROUNDOFF defined in sundials_types.h.
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*CVDlsDenseJacFn)(long int N, realtype t,
			       N_Vector y, N_Vector fy, 
			       DlsMat Jac, void *user_data,
			       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
/*
 * -----------------------------------------------------------------
 * Type: CVDlsBandJacFn
 * -----------------------------------------------------------------
 *
 * A band Jacobian approximation function Jac must have the
 * prototype given below. Its parameters are:
 *
 * N is the length of all vector arguments.
 *
 * mupper is the upper half-bandwidth of the approximate banded
 * Jacobian. This parameter is the same as the mupper parameter
 * passed by the user to the linear solver initialization function.
 *
 * mlower is the lower half-bandwidth of the approximate banded
 * Jacobian. This parameter is the same as the mlower parameter
 * passed by the user to the linear solver initialization function.
 *
 * t is the current value of the independent variable.
 *
 * y is the current value of the dependent variable vector,
 *      namely the predicted value of y(t).
 *
 * fy is the vector f(t,y).
 *
 * Jac is the band matrix (of type DlsMat) that will be loaded
 * by a CVDlsBandJacFn with an approximation to the Jacobian matrix
 * Jac = (df_i/dy_j) at the point (t,y).
 * Three efficient ways to load J are:
 *
 * (1) (with macros - no explicit data structure references)
 *    for (j=0; j < n; j++) {
 *       col_j = BAND_COL(Jac,j);
 *       for (i=j-mupper; i <= j+mlower; i++) {
 *         generate J_ij = the (i,j)th Jacobian element
 *         BAND_COL_ELEM(col_j,i,j) = J_ij;
 *       }
 *     }
 *
 * (2) (with BAND_COL macro, but without BAND_COL_ELEM macro)
 *    for (j=0; j < n; j++) {
 *       col_j = BAND_COL(Jac,j);
 *       for (k=-mupper; k <= mlower; k++) {
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k
 *         col_j[k] = J_ij;
 *       }
 *     }
 *
 * (3) (without macros - explicit data structure references)
 *     offset = Jac->smu;
 *     for (j=0; j < n; j++) {
 *       col_j = ((Jac->data)[j])+offset;
 *       for (k=-mupper; k <= mlower; k++) {
 *         generate J_ij = the (i,j)th Jacobian element, i=j+k
 *         col_j[k] = J_ij;
 *       }
 *     }
 * Caution: Jac->smu is generally NOT the same as mupper.
 *
 * The BAND_ELEM(A,i,j) macro is appropriate for use in small
 * problems in which efficiency of access is NOT a major concern.
 *
 * user_data is a pointer to user data - the same as the user_data
 *          parameter passed to CVodeSetFdata.
 *
 * NOTE: If the user's Jacobian routine needs other quantities,
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through
 *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively
 *     (see cvode.h). The unit roundoff is available as
 *     UNIT_ROUNDOFF defined in sundials_types.h
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for
 * vectors of length N which can be used by a CVDlsBandJacFn
 * as temporary storage or work space.
 *
 * A CVDlsBandJacFn should return 0 if successful, a positive value
 * if a recoverable error occurred, and a negative value if an 
 * unrecoverable error occurred.
 * -----------------------------------------------------------------
 */

typedef int (*CVDlsBandJacFn)(long int N, long int mupper, long int mlower,
			      realtype t, N_Vector y, N_Vector fy, 
			      DlsMat Jac, void *user_data,
			      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * -----------------------------------------------------------------
 * EXPORTED FUNCTIONS 
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVDLS linear solver
 * -----------------------------------------------------------------
 *
 * CVDlsSetDenseJacFn specifies the dense Jacobian approximation
 * routine to be used for a direct dense linear solver.
 *
 * CVDlsSetBandJacFn specifies the band Jacobian approximation
 * routine to be used for a direct band linear solver.
 *
 * By default, a difference quotient approximation, supplied with
 * the solver is used.
 *
 * The return value is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the CVODE memory was NULL
 *    CVDLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn jac);
SUNDIALS_EXPORT int CVDlsSetBandJacFn(void *cvode_mem, CVDlsBandJacFn jac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVSDIRECT linear solver
 * -----------------------------------------------------------------
 *
 * CVDlsGetWorkSpace   returns the real and integer workspace used
 *                     by the direct linear solver.
 * CVDlsGetNumJacEvals returns the number of calls made to the
 *                     Jacobian evaluation routine jac.
 * CVDlsGetNumRhsEvals returns the number of calls to the user
 *                     f routine due to finite difference Jacobian
 *                     evaluation.
 * CVDlsGetLastFlag    returns the last error flag set by any of
 *                     the CVSDIRECT interface functions.
 *
 * The return value of CVDlsGet* is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the CVODES memory was NULL
 *    CVDLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDlsGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals);
SUNDIALS_EXPORT int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS);
SUNDIALS_EXPORT int CVDlsGetLastFlag(void *cvode_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a CVSDIRECT return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *CVDlsGetReturnFlagName(long int flag);

/*
 * =================================================================
 * PART II:  B A C K W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * FUNCTION TYPES
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Type: CVDlsDenseJacFnB
 * -----------------------------------------------------------------
 * A dense Jacobian approximation function jacB for the adjoint
 * (backward) problem must have the prototype given below. 
 * -----------------------------------------------------------------
 */

typedef int (*CVDlsDenseJacFnB)(long int nB, realtype t,
				N_Vector y, 
				N_Vector yB, N_Vector fyB,
				DlsMat JB, void *user_dataB, 
				N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * Type: CVDlsDenseJacFnBS
 * -----------------------------------------------------------------
 * A dense Jacobian approximation function jacBS for the adjoint
 * (backward) problem, sensitivity-dependent case,  must have the
 *  prototype given below. 
 * -----------------------------------------------------------------
 */

typedef int (*CVDlsDenseJacFnBS)(long int nB, realtype t,
				N_Vector y, N_Vector *yS,
				N_Vector yB, N_Vector fyB,
				DlsMat JB, void *user_dataB, 
				N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * Type : CVDlsBandJacFnB
 * -----------------------------------------------------------------
 * A band Jacobian approximation function jacB for the adjoint 
 * (backward) problem must have the prototype given below. 
 * -----------------------------------------------------------------
 */

typedef int (*CVDlsBandJacFnB)(long int nB, long int mupperB, long int mlowerB,
			       realtype t, N_Vector y,
			       N_Vector yB, N_Vector fyB,
			       DlsMat JB, void *user_dataB, 
			       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * Type : CVDlsBandJacFnBS
 * -----------------------------------------------------------------
 * A band Jacobian approximation function jacBS for the adjoint 
 * (backward) problem, sensitivity-dependent case, must have the
 * prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVDlsBandJacFnBS)(long int nB, long int mupperB, long int mlowerB,
			       realtype t, N_Vector y, N_Vector *yS,
			       N_Vector yB, N_Vector fyB,
			       DlsMat JB, void *user_dataB,
			       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * EXPORTED FUNCTIONS 
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Functions: CVDlsSet*JacFnB and CVDlsSet*JacFnBS
 * -----------------------------------------------------------------
 * CVDlsSetDenseJacFnB and CVDlsSetBandJacFnB specify the dense and
 * band Jacobian functions, respectively, to be used by a
 * CVSDIRECT linear solver for the backward integration phase, when
 * the backward problem does not depend on forward sensitivities.
 * CVDlsSetDenseJacFnBS and CVDlsSetBandJacFnBS specify the Jacobian
 * functions when the backward problem does depend on sensitivities.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDlsSetDenseJacFnB(void *cvode_mem, int which,
                                        CVDlsDenseJacFnB jacB);

SUNDIALS_EXPORT int CVDlsSetDenseJacFnBS(void *cvode_mem, int which,
                                         CVDlsDenseJacFnBS jacBS);

SUNDIALS_EXPORT int CVDlsSetBandJacFnB(void *cvode_mem, int which,
                                       CVDlsBandJacFnB jacB);

SUNDIALS_EXPORT int CVDlsSetBandJacFnBS(void *cvode_mem, int which,
                                        CVDlsBandJacFnBS jacBS);

#ifdef __cplusplus
}
#endif

#endif
