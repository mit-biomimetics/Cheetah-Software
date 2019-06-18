#if !defined(DSDPLAPACK)
#define DSDPLAPACK
/*!
\file dsdplapack.h
\brief DSDP uses BLAS and LAPACK for many of its operations.
*/

typedef long int ffinteger;
/*
typedef int ffinteger;
*/
/*
#define  __DSDP_NONAMEMANGLING
#undef  __DSDP_NONAMEMANGLING
*/


#ifdef __cplusplus
#define  __DSDP_NONAMEMANGLING
#endif


#ifdef __DSDP_NONAMEMANGLING
#define EXTERN_C_BEGIN extern "C" {
#define EXTERN_C_END }
#else
#define EXTERN_C_BEGIN
#define EXTERN_C_END
#endif

/*
#define NOUNDERBLAS
#define CAPSBLAS
*/

#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
#define dstev   DSTEV
#define dpotrf  DPOTRF
#define dtrsm   DTRSM
#define dsyev   DSYEV
#define dpotrs  DPOTRS
#define daxpy   DAXPY
#define dgemv   DGEMV
#define dscal   DSCAL
#define dger    DGER
#define dsymv   DSYMV
#define dasum   DASUM
#define ddot    DDOT
#define dnrm2   DNRM2
#define dspmv   DSPMV
#define dspr    DSPR
#define dpptrf  DPPTRF
#define dpptrs  DPPTRS
#define dtpsv   DTPSV
#define dspevd  DSPEVD
#define dtrsv   DTRSV
#define dsyr    DSYR
#define dtrmv   DTRMV
#define dpotri  DPOTRI
#define dpptri  DPPTRI
#define dsyevx  DSYEVX
#define dsyevd  DSYEVD
#define dspevx  DSPEVX
#define dsyevr  DSYEVR
#define dstevr  DSTEVR
#endif
#endif

#ifndef NOUNDERBLAS
#ifdef CAPSBLAS
#define dstev   DSTEV_
#define dpotrf  DPOTRF_
#define dtrsm   DTRSM_
#define dsyev   DSYEV_
#define dpotrs  DPOTRS_
#define daxpy   DAXPY_
#define dgemv   DGEMV_
#define dscal   DSCAL_
#define dger    DGER_
#define dsymv   DSYMV_
#define ddot    DDOT_
#define dnrm2   DNRM2_
#define dasum   DASUM_
#define dspmv   DSPMV_
#define dspr    DSPR_
#define dpptrf  DPPTRF_
#define dpptrs  DPPTRS_
#define dtpsv   DTPSV_
#define dspevd  DSPEVD_
#define dtrsv   DTRSV_
#define dsyr    DSYR_
#define dtrmv   DTRMV_
#define dpotri  DPOTRI_
#define dpptri  DPPTRI_
#define dsyevx  DSYEVX_
#define dsyevd  DSYEVD_
#define dspevx  DSPEVX_
#define dsyevr  DSYEVR_
#define dstevr  DSTEVR_
#endif
#endif

#ifdef NOUNDERBLAS
#ifndef CAPSBLAS
#define dstev   dstev
#define dpotrf  dpotrf
#define dtrsm   dtrsm
#define dsyev   dsyev
#define dpotrs  dpotrs
#define daxpy   daxpy
#define dgemv   dgemv
#define dscal   dscal
#define dger    dger
#define dsymv   dsymv
#define dasum   dasum
#define ddot    ddot
#define dnrm2   dnrm2
#define dspmv   dspmv
#define dspr    dspr
#define dpptrf  dpptrf
#define dpptrs  dpptrs
#define dtpsv   dtpsv
#define dspevd  dspevd
#define dtrsv   dtrsv
#define dsyr    dsyr
#define dtrmv   dtrmv
#define dpotri  dpotri
#define dpptri  dpptri
#define dsyevx  dsyevx
#define dsyevd  dsyevd
#define dspevx  dspevx
#define dsyevr  dsyevr
#define dstevr  dstevr
#endif
#endif

#ifndef NOUNDERBLAS
#ifndef CAPSBLAS
#define dstev   dstev_
#define dpotrf  dpotrf_
#define dtrsm   dtrsm_
#define dsyev   dsyev_
#define dpotrs  dpotrs_
#define daxpy   daxpy_
#define dgemv   dgemv_
#define dscal   dscal_
#define dger    dger_
#define dsymv   dsymv_
#define dasum   dasum_
#define ddot    ddot_
#define dnrm2   dnrm2_
#define dspmv   dspmv_
#define dspr    dspr_
#define dpptrf  dpptrf_
#define dpptrs  dpptrs_
#define dtpsv   dtpsv_
#define dspevd  dspevd_
#define dtrsv   dtrsv_
#define dsyr    dsyr_
#define dtrmv   dtrmv_
#define dpotri  dpotri_
#define dpptri  dpptri_
#define dsyevx  dsyevx_
#define dsyevd  dsyevd_
#define dspevx  dspevx_
#define dsyevr  dsyevr_
#define dstevr  dstevr_
#endif
#endif

EXTERN_C_BEGIN

void dpotrs(char*,ffinteger*,ffinteger*,double*,ffinteger*,double*,ffinteger*,ffinteger*);                             /* Cholesky Solve */
void dpotrf(char*,ffinteger*,double*,ffinteger*,ffinteger*);                                               /* Cholesky Factor */
void dtrsm(char*,char*,char*,char*,ffinteger*,ffinteger*,double*,double*,ffinteger*,double*,ffinteger*);         /* Cholesky trianglular solve */ 
void dsyev(char*,char*,ffinteger*,double*,ffinteger*,double*,double*,ffinteger*,ffinteger*);                     /* Compute eigenvalues/vectors */
void dstev(char*,ffinteger*,double*,double*,double*,ffinteger*,double*,ffinteger*);                     /* Compute eigenvalues/vectors */
void dgemv(char*,ffinteger*,ffinteger*,double*,double*,ffinteger*,double*,ffinteger*,double*,double*,ffinteger*);
void dspmv(char*,ffinteger*,double*,double*,double*,ffinteger*,double*,double*,ffinteger*);

void dspr(char*,ffinteger*,double*,double*,ffinteger*,double*);
void dpptrs(char*,ffinteger*,ffinteger*,double*,double*,ffinteger*,ffinteger*);
void dpptrf(char*,ffinteger*,double*,ffinteger*);

void dtrsv(char*,char*,char*,ffinteger*,double*,ffinteger*,double*,ffinteger*);
void dsyr(char*,ffinteger*,double*,double*,ffinteger*,double*,ffinteger*);
void dtrmv(char*,char*,char*,ffinteger*,double*,ffinteger*,double*,ffinteger*);

void dtpsv(char*,char*,char*,ffinteger*,double*,double*,ffinteger*);
void dger(ffinteger*,ffinteger*,double*,double*,ffinteger*,double*,ffinteger*,double*,ffinteger*);
void dsymv(char*,ffinteger*,double*,double*,ffinteger*,double*,ffinteger*, double*,double*,ffinteger*);
void dspevd(char*,char*,ffinteger*,double*,double*,double*,ffinteger*,double*,ffinteger*,ffinteger*,ffinteger*,ffinteger*);

double dasum(ffinteger*,double*,ffinteger*);
void dscal(ffinteger*,double*,double*,ffinteger*);
void daxpy(ffinteger*,double*,double*,ffinteger*,double*,ffinteger*);
double ddot(ffinteger*,double*,ffinteger*,double*,ffinteger*);
double dnrm2(ffinteger*,double*,ffinteger*);

void dpotri(char*,ffinteger*,double*,ffinteger*,ffinteger*);
void dpptri(char*,ffinteger*,double*,ffinteger*);

void dsyevx(char*, char*, char*, ffinteger*,double*, ffinteger*, double*, double*, ffinteger*, ffinteger*, double*, ffinteger*, double *,double*, ffinteger *, double *, ffinteger*, ffinteger*, ffinteger*, ffinteger*);

void dspevx(char*, char*, char*, ffinteger*,double*, double*, double*, ffinteger*, ffinteger*, double*, ffinteger*, double *,double*, ffinteger *, double *, ffinteger*, ffinteger*, ffinteger*);

void dsdevx(char*, char*, char*, ffinteger*,double*, double*, double*, ffinteger*, ffinteger*, double*, ffinteger*, double *,double*, ffinteger *, double *, ffinteger*, ffinteger*, ffinteger*, ffinteger*);


void dsyevr(char*, char*, char*, ffinteger*, double*, ffinteger*, double*, double*, ffinteger*, ffinteger*,double*, ffinteger*, double*, double*, ffinteger*, ffinteger*, double*, ffinteger*, ffinteger*, ffinteger*, ffinteger*);

void dstevr(char*, char*, ffinteger*, double*, double*, double*, double*, ffinteger*, ffinteger*,double*, ffinteger*, double*, double*, ffinteger*, ffinteger*, double*, ffinteger*, ffinteger*, ffinteger*, ffinteger*);

EXTERN_C_END

#endif
