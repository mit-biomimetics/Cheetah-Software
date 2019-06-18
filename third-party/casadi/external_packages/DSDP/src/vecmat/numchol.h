#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef enum {
  CfcOk=0,
  CfcSpace,   /* fail to allocate required space */
  CfcIndef    /* indefinity is detected          */
} cfc_sta;

typedef struct {
  int    mrow;     /* number of rows allocated                    */
  int    nrow;     /* number of rows used                         */
    
  int    snnz;     /* number of indices for nonzeros in S         */
  int    *shead;   /* position of first nonzero in row i of S     */
  int    *ssize;   /* number of non-zeros in row i of S below     */
                   /* the diagonal                                */
  int    *ssub;    /* column index buffer for non-zeros in S      */
  double *diag;    /* diagonal matrix D in the factorization      */
  double *sqrtdiag;/* sqrt o diagonal matrix D in the factorization */
    
  int    unnz;     /* number of nonzeros in the upper factor      */
  int    ujnz;     /* number of column indices in the compressed  */
                     /* indices buffer ujsub                        */ 
  int    *ujbeg;   /* beginning position of indices in row i of U */
  int    *uhead;   /* position of first nonzero in row i of U     */
  int    *ujsze;   /* number of indices in row i of U             */ 
  int    *usub;    /* compressed column index buffer of U         */
  double *uval;    /* nonzero values in factor U                  */
    
  int    *perm;    /* permutation order                           */
  int    *invp;    /* inverse order of perm                       */
    
  int    nsnds;    /* number of supernodes                        */
  int    *subg;    /* index of the first column in supernode i    */
  int    ndens;    /* numer of dense rows                         */
  int    nsndn;    /* number supernodes in dense rows             */
  int    *dhead;   /* pointer first column in each dense row      */
  int    *dsub;    /* indices in dense rows                       */
  int    *dbeg;    /* beginning of column index                   */ 
  int    sdens;    /* separate dense row                          */

  int    alldense;

  double tolpiv;
  int cachesize;
  int cacheunit;

  /* New */
  int n;
  int *iw;
  double *rw;
  int factor;
} chfac;


typedef struct {
    int idep;
    int last;
    int most;
    int cure;
    int loca;
    int lowp;
    int ntot;
    
    int *head;
    int *port;
    int *fwrd;
    int *bwrd;
  } xlist;
  
typedef struct {
    int nnod;
    int nn0;
    int raft;
    int head;
    int last;
    int ntot;
    
    int *adjn;
    int *rbeg;
    int *rexs;
    int *rlen;
    int *rend;
    int *pres;
    int *succ;
  } order;

typedef enum {
  OptFound=0,
  SysError=100,
  OutOfSpc,CholErr
} xcode;

#if !defined (min)
#define min(a,b) ((a <= b)? (a) : (b))
#endif
#if !defined (max)
#define max(a,b) ((a >= b)? (a) : (b))
#endif
#if !defined (sign)
#define sign(a) ((a<0)? (-1) : (1))
#endif
#if !defined (TRUE)
#define TRUE 1
#endif
#if !defined (FALSE)
#define FALSE 0
#endif

#include "sdpfun.h"
