#ifndef _CS_H
#define _CS_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#define CS_VER 2                    /* CSparse Version */
#define CS_SUBVER 2
#define CS_SUBSUB 4
#define CS_DATE "Nov 30, 2009"     /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006-2009"

/* --- primary CSparse routines and data structures ------------------------- */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    int nzmax ;     /* maximum number of entries */
    int m ;         /* number of rows */
    int n ;         /* number of columns */
    int *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    int nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

#ifdef __cplusplus
extern "C" {
#endif

  cs *cs_add (const cs *A, const cs *B, double alpha, double beta) ;
  int cs_cholsol (int order, const cs *A, double *b) ;
  cs *cs_compress (const cs *T) ;
  int cs_dupl (cs *A) ;
  int cs_entry (cs *T, int i, int j, double x) ;
  int cs_gaxpy (const cs *A, const double *x, double *y) ;
  cs *cs_load (FILE *f) ;
  int cs_lusol (int order, const cs *A, double *b, double tol) ;
  cs *cs_multiply (const cs *A, const cs *B) ;
  double cs_norm (const cs *A) ;
  int cs_print (const cs *A, int brief) ;
  int cs_qrsol (int order, const cs *A, double *b) ;
  cs *cs_transpose (const cs *A, int values) ;
  /* utilities */
  void *cs_calloc (int n, size_t size) ;
  void *cs_free (void *p) ;
  void *cs_realloc (void *p, int n, size_t size, int *ok) ;
  cs *cs_spalloc (int m, int n, int nzmax, int values, int triplet) ;
  cs *cs_spfree (cs *A) ;
  int cs_sprealloc (cs *A, int nzmax) ;
  void *cs_malloc (int n, size_t size) ;

  /* --- secondary CSparse routines and data structures ----------------------- */
  typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
  {
      int *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
      int *q ;        /* fill-reducing column permutation for LU and QR */
      int *parent ;   /* elimination tree for Cholesky and QR */
      int *cp ;       /* column pointers for Cholesky, row counts for QR */
      int *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
      int m2 ;        /* # of rows for QR, after adding fictitious rows */
      double lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
      double unz ;    /* # entries in U for LU; in R for QR */
  } css ;

  typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
  {
      cs *L ;         /* L for LU and Cholesky, V for QR */
      cs *U ;         /* U for LU, R for QR, not used for Cholesky */
      int *pinv ;     /* partial pivoting for LU */
      double *B ;     /* beta [0..n-1] for QR */
  } csn ;

  typedef struct cs_dmperm_results    /* cs_dmperm or cs_scc output */
  {
      int *p ;        /* size m, row permutation */
      int *q ;        /* size n, column permutation */
      int *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
      int *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
      int nb ;        /* # of blocks in fine dmperm decomposition */
      int rr [5] ;    /* coarse row decomposition */
      int cc [5] ;    /* coarse column decomposition */
  } csd ;

  int *cs_amd (int order, const cs *A) ;
  csn *cs_chol (const cs *A, const css *S) ;
  csd *cs_dmperm (const cs *A, int seed) ;
  int cs_droptol (cs *A, double tol) ;
  int cs_dropzeros (cs *A) ;
  int cs_happly (const cs *V, int i, double beta, double *x) ;
  int cs_ipvec (const int *p, const double *b, double *x, int n) ;
  int cs_lsolve (const cs *L, double *x) ;
  int cs_ltsolve (const cs *L, double *x) ;
  csn *cs_lu (const cs *A, const css *S, double tol) ;
  cs *cs_permute (const cs *A, const int *pinv, const int *q, int values) ;
  int *cs_pinv (const int *p, int n) ;
  int cs_pvec (const int *p, const double *b, double *x, int n) ;
  csn *cs_qr (const cs *A, const css *S) ;
  css *cs_schol (int order, const cs *A) ;
  css *cs_sqr (int order, const cs *A, int qr) ;
  cs *cs_symperm (const cs *A, const int *pinv, int values) ;
  int cs_updown (cs *L, int sigma, const cs *C, const int *parent) ;
  int cs_usolve (const cs *U, double *x) ;
  int cs_utsolve (const cs *U, double *x) ;
  /* utilities */
  css *cs_sfree (css *S) ;
  csn *cs_nfree (csn *N) ;
  csd *cs_dfree (csd *D) ;

  /* --- tertiary CSparse routines -------------------------------------------- */
  int *cs_counts (const cs *A, const int *parent, const int *post, int ata) ;
  double cs_cumsum (int *p, int *c, int n) ;
  int cs_dfs (int j, cs *G, int top, int *xi, int *pstack, const int *pinv) ;
  int cs_ereach (const cs *A, int k, const int *parent, int *s, int *w) ;
  int *cs_etree (const cs *A, int ata) ;
  int cs_fkeep (cs *A, int (*fkeep) (int, int, double, void *), void *other) ;
  double cs_house (double *x, double *beta, int n) ;
  int cs_leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
       int *ancestor, int *jleaf) ;
  int *cs_maxtrans (const cs *A, int seed) ;
  int *cs_post (const int *parent, int n) ;
  int *cs_randperm (int n, int seed) ;
  int cs_reach (cs *G, const cs *B, int k, int *xi, const int *pinv) ;
  int cs_scatter (const cs *A, int j, double beta, int *w, double *x, int mark,
      cs *C, int nz) ;
  csd *cs_scc (cs *A) ;
  int cs_spsolve (cs *G, const cs *B, int k, int *xi, double *x,
      const int *pinv, int lo) ;
  int cs_tdfs (int j, int k, int *head, const int *next, int *post,
      int *stack) ;
  /* utilities */
  csd *cs_dalloc (int m, int n) ;
  csd *cs_ddone (csd *D, cs *C, void *w, int ok) ;
  cs *cs_done (cs *C, void *w, void *x, int ok) ;
  int *cs_idone (int *p, cs *C, void *w, int ok) ;
  csn *cs_ndone (csn *N, cs *C, void *w, void *x, int ok) ;

#ifdef __cplusplus
}
#endif

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))
#endif
