/* ========================================================================= */
/* === amd_internal.h ====================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,                                      */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* This file is for internal use in AMD itself, and does not normally need to
 * be included in user code (it is included in UMFPACK, however).   All others
 * should use amd.h instead.
 */

/* ========================================================================= */
/* === NDEBUG ============================================================== */
/* ========================================================================= */

/*
 * Turning on debugging takes some work (see below).   If you do not edit this
 * file, then debugging is always turned off, regardless of whether or not
 * -DNDEBUG is specified in your compiler options.
 *
 * If AMD is being compiled as a mexFunction, then MATLAB is defined,
 * and mxAssert is used instead of assert.  If debugging is not enabled, no
 * MATLAB include files or functions are used.  Thus, the AMD library libamd.a
 * can be safely used in either a stand-alone C program or in another
 * mexFunction, without any change.
 */

/*
    AMD will be exceedingly slow when running in debug mode.  The next three
    lines ensure that debugging is turned off.
*/
#ifndef NDEBUG
#define NDEBUG
#endif

#include "glob_opts.h"

/*
    To enable debugging, uncomment the following line:
#undef NDEBUG
*/

/* ------------------------------------------------------------------------- */
/* ANSI include files */
/* ------------------------------------------------------------------------- */

/* from stdlib.h:  size_t, malloc, free, realloc, and calloc */
#include <stdlib.h>

#if !defined(NPRINT) || !defined(NDEBUG)
/* from stdio.h:  printf.  Not included if NPRINT is defined at compile time.
 * fopen and fscanf are used when debugging. */
#include <stdio.h>
#endif

/* from limits.h:  INT_MAX and LONG_MAX */
#include <limits.h>

/* from math.h: sqrt */
#include <math.h>

/* ------------------------------------------------------------------------- */
/* MATLAB include files (only if being used in or via MATLAB) */
/* ------------------------------------------------------------------------- */

#ifdef MATLAB
#include "matrix.h"
#include "mex.h"
#endif

/* ------------------------------------------------------------------------- */
/* basic definitions */
/* ------------------------------------------------------------------------- */

#ifdef FLIP
#undef FLIP
#endif

#ifdef MAX
#undef MAX
#endif

#ifdef MIN
#undef MIN
#endif

#ifdef EMPTY
#undef EMPTY
#endif

#ifdef GLOBAL
#undef GLOBAL
#endif

#ifdef PRIVATE
#undef PRIVATE
#endif

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) ((i < EMPTY) ? FLIP (i) : (i))

/* for integer MAX/MIN, or for c_floats when we don't care how NaN's behave: */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* logical expression of p implies q: */
#define IMPLIES(p,q) (!(p) || (q))

/* Note that the IBM RS 6000 xlc predefines TRUE and FALSE in <types.h>. */
/* The Compaq Alpha also predefines TRUE and FALSE. */
#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

#define TRUE (1)
#define FALSE (0)
#define PRIVATE static
#define GLOBAL
#define EMPTY (-1)

/* Note that Linux's gcc 2.96 defines NULL as ((void *) 0), but other */
/* compilers (even gcc 2.95.2 on Solaris) define NULL as 0 or (0).  We */
/* need to use the ANSI standard value of 0. */
#ifdef NULL
#undef NULL
#endif

#define NULL 0

/* largest value of size_t */
#ifndef SIZE_T_MAX
#ifdef SIZE_MAX
/* C99 only */
#define SIZE_T_MAX SIZE_MAX
#else
#define SIZE_T_MAX ((size_t) (-1))
#endif
#endif

/* ------------------------------------------------------------------------- */
/* integer type for AMD: int or SuiteSparse_long */
/* ------------------------------------------------------------------------- */

#include "amd.h"

#if defined (DLONG) || defined (ZLONG)

#define Int SuiteSparse_long
#define ID  SuiteSparse_long_id
#define Int_MAX SuiteSparse_long_max

#define AMD_order amd_l_order
#define AMD_defaults amd_l_defaults
#define AMD_control amd_l_control
#define AMD_info amd_l_info
#define AMD_1 amd_l1
#define AMD_2 amd_l2
#define AMD_valid amd_l_valid
#define AMD_aat amd_l_aat
#define AMD_postorder amd_l_postorder
#define AMD_post_tree amd_l_post_tree
#define AMD_dump amd_l_dump
#define AMD_debug amd_l_debug
#define AMD_debug_init amd_l_debug_init
#define AMD_preprocess amd_l_preprocess

#else

#define Int int
#define ID "%d"
#define Int_MAX INT_MAX

#define AMD_order amd_order
#define AMD_defaults amd_defaults
#define AMD_control amd_control
#define AMD_info amd_info
#define AMD_1 amd_1
#define AMD_2 amd_2
#define AMD_valid amd_valid
#define AMD_aat amd_aat
#define AMD_postorder amd_postorder
#define AMD_post_tree amd_post_tree
#define AMD_dump amd_dump
#define AMD_debug amd_debug
#define AMD_debug_init amd_debug_init
#define AMD_preprocess amd_preprocess

#endif

/* ------------------------------------------------------------------------- */
/* AMD routine definitions (not user-callable) */
/* ------------------------------------------------------------------------- */

GLOBAL size_t AMD_aat
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Len [ ],
    Int Tp [ ],
    c_float Info [ ]
) ;

GLOBAL void AMD_1
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int P [ ],
    Int Pinv [ ],
    Int Len [ ],
    Int slen,
    Int S [ ],
    c_float Control [ ],
    c_float Info [ ]
) ;

GLOBAL void AMD_postorder
(
    Int nn,
    Int Parent [ ],
    Int Npiv [ ],
    Int Fsize [ ],
    Int Order [ ],
    Int Child [ ],
    Int Sibling [ ],
    Int Stack [ ]
) ;

GLOBAL Int AMD_post_tree
(
    Int root,
    Int k,
    Int Child [ ],
    const Int Sibling [ ],
    Int Order [ ],
    Int Stack [ ]
#ifndef NDEBUG
    , Int nn
#endif
) ;

GLOBAL void AMD_preprocess
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Rp [ ],
    Int Ri [ ],
    Int W [ ],
    Int Flag [ ]
) ;

/* ------------------------------------------------------------------------- */
/* debugging definitions */
/* ------------------------------------------------------------------------- */

#ifndef NDEBUG

/* from assert.h:  assert macro */
#include <assert.h>

#ifndef EXTERN
#define EXTERN extern
#endif

EXTERN Int AMD_debug ;

GLOBAL void AMD_debug_init ( char *s ) ;

GLOBAL void AMD_dump
(
    Int n,
    Int Pe [ ],
    Int Iw [ ],
    Int Len [ ],
    Int iwlen,
    Int pfree,
    Int Nv [ ],
    Int Next [ ],
    Int Last [ ],
    Int Head [ ],
    Int Elen [ ],
    Int Degree [ ],
    Int W [ ],
    Int nel
) ;

#ifdef ASSERT
#undef ASSERT
#endif

/* Use mxAssert if AMD is compiled into a mexFunction */
#ifdef MATLAB
#define ASSERT(expression) (mxAssert ((expression), ""))
#else
#define ASSERT(expression) (assert (expression))
#endif

#define AMD_DEBUG0(params) { SUITESPARSE_PRINTF (params) ; }
#define AMD_DEBUG1(params) { if (AMD_debug >= 1) SUITESPARSE_PRINTF (params) ; }
#define AMD_DEBUG2(params) { if (AMD_debug >= 2) SUITESPARSE_PRINTF (params) ; }
#define AMD_DEBUG3(params) { if (AMD_debug >= 3) SUITESPARSE_PRINTF (params) ; }
#define AMD_DEBUG4(params) { if (AMD_debug >= 4) SUITESPARSE_PRINTF (params) ; }

#else

/* no debugging */
#define ASSERT(expression)
#define AMD_DEBUG0(params)
#define AMD_DEBUG1(params)
#define AMD_DEBUG2(params)
#define AMD_DEBUG3(params)
#define AMD_DEBUG4(params)

#endif
