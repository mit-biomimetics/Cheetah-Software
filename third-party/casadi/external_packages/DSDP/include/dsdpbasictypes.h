#ifndef __DSDP_BASIC_TYPES
#define __DSDP_BASIC_TYPES
/*!
  \file dsdpbasictypes.h
  \brief Solver, solution types, termination codes, 
*/

/*! 
\typedef struct DSDP_C* DSDP;
\brief An implementation of
the dual-scaling algorithm for semidefinite programming.
*/
typedef struct DSDP_C* DSDP;

/*! 
\typedef enum DSDPTruth
\brief Boolean variables
*/
typedef enum { DSDP_FALSE = 0, /*!< 0*/ DSDP_TRUE = 1/*!< 1 */} DSDPTruth;

/*! 
\typedef enum DSDPDualFactorMatrix
\brief  DSDP requires two instances of the data structures S.
*/
typedef enum {
  DUAL_FACTOR           =  1, /*!< First instance for dual variable S */
  PRIMAL_FACTOR         =  2  /*!< Second instance used to compute X */
} DSDPDualFactorMatrix;

typedef enum { DSDPAlways=1, DSDPNever=2, DSDPInfeasible=0} DSDPPenalty;

/*! 
\typedef enum DSDPSolutionType
\brief Formulations (P) and (D) can be feasible and bounded, feasible
and unbounded, or infeasible.
\sa DSDPGetSolutionType()
*/
typedef enum {/* converged */
  DSDP_PDUNKNOWN          = 0, /*!<  Not sure whether (D) or (P) is feasible, check y bounds */
  DSDP_PDFEASIBLE         = 1, /*!<  Both (D) and (P) are feasible and bounded */
  DSDP_UNBOUNDED          = 3, /*!<  (D) is unbounded and (P) is infeasible  */
  DSDP_INFEASIBLE         = 4  /*!<  (D) in infeasible and (P) is unbounded  */
} DSDPSolutionType;

/*! 
\typedef enum DSDPTerminationReason
\brief There are many reasons to terminate the solver.
\sa DSDPStopReason()
*/
typedef enum {
  DSDP_CONVERGED          =  1, /*!< Good news: Solution found. */
  DSDP_INFEASIBLE_START   = -6, /*!< The initial points y and r imply that S is not positive*/
  DSDP_SMALL_STEPS        = -2, /*!< Short step lengths created by numerical difficulties prevent progress */
  DSDP_INDEFINITE_SCHUR_MATRIX  = -8,  /*!< Theoretically this matrix is positive definite */
  DSDP_MAX_IT             = -3, /*!< Reached maximum number of iterations */
  DSDP_NUMERICAL_ERROR   = -9,  /*!< Another numerical error occurred. Check solution */
  DSDP_UPPERBOUND         =  5, /*!< Objective (DD) big enough to stop */
  DSDP_USER_TERMINATION   =  7,  /*!< DSDP didn't stop it, did you? */
  CONTINUE_ITERATING      =  0  /*!< Don't Stop */ } DSDPTerminationReason;

extern int DSDPSetConvergenceFlag(DSDP,DSDPTerminationReason);

#endif
