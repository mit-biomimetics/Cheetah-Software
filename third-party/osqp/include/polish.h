/* Solution polish based on assuming the active set */
#ifndef POLISH_H
# define POLISH_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


# include "types.h"

/**
 * Solution polish: Solve equality constrained QP with assumed active
 *constraints
 * @param  work Workspace
 * @return      Exitflag:  0: Factorization successfull
 *                         1: Factorization unsuccessfull
 */
c_int polish(OSQPWorkspace *work);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef POLISH_H
