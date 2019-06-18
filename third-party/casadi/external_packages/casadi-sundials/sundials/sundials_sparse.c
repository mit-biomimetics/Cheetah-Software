/*
 * -----------------------------------------------------------------
 * $Revision: 4433 $
 * $Date: 2015-03-23 18:24:01 -0700 (Mon, 23 Mar 2015) $
 * -----------------------------------------------------------------
 * Programmers: Carol Woodward @ LLNL
 *              Daniel R. Reynolds @ SMU
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
 * This is the implementation file for operations on the SUNDIALS
 * sparse matrix structure.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_sparse.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* 
Creates a new (empty) sparse matrix of a desired size and nonzero density.
Returns NULL if a memory allocation error occurred.
*/
SlsMat NewSparseMat(int M, int N, int NNZ)
{
  SlsMat A;

  if ( (M <= 0) || (N <= 0) ) return(NULL);

  A = NULL;
  A = (SlsMat) malloc(sizeof(struct _SlsMat));
  if (A==NULL) return (NULL);
  
  A->data = (realtype *) malloc(NNZ * sizeof(realtype));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }
  A->rowvals = (int *) malloc(NNZ * sizeof(int));
  if (A->rowvals == NULL) {
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }
  A->colptrs = (int *) malloc((N+1) * sizeof(int));
  if (A->colptrs == NULL) {
    free(A->rowvals);
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }

  A->M = M;
  A->N = N;
  A->NNZ = NNZ;
  /* A->colptrs[N] = NNZ; */
  A->colptrs[N] = 0;

  return(A);
}

/* 
Utility routine to create a new sparse matrix out of an 
existing dense or band matrix.  Returns NULL if a memory 
allocation error occurred.
*/
SlsMat SlsConvertDls(DlsMat A)
{
  int i, j, nnz;
  realtype dtmp;
  SlsMat As = NULL;
  
  /* proceed according to A's type (dense/band) */
  if (A->type == SUNDIALS_DENSE) {

    /* determine total number of nonzeros */
    nnz = 0;
    for (j=0; j<A->N; j++)
      for (i=0; i<A->M; i++)
	nnz += (DENSE_ELEM(A,i,j) != 0.0);

    /* allocate sparse matrix */
    As = NewSparseMat(A->M, A->N, nnz);
    if (As == NULL)  return NULL;

    /* copy nonzeros from A into As */
    nnz = 0;
    for (j=0; j<A->N; j++) {
      As->colptrs[j] = nnz;
      for (i=0; i<A->M; i++) {
	dtmp = DENSE_ELEM(A,i,j);
	if ( dtmp != 0.0 ) { 
	  As->rowvals[nnz] = i;
	  As->data[nnz++] = dtmp;
	}
      }
    }
    As->colptrs[A->N] = nnz;

  } else { /* SUNDIALS_BAND */

    /* determine total number of nonzeros */
    nnz = 0;
    for (j=0; j<A->N; j++)
      for (i=j-(A->mu); i<j+(A->ml); i++)
	nnz += (BAND_ELEM(A,i,j) != 0.0);

    /* allocate sparse matrix */
    As = NewSparseMat(A->M, A->N, nnz);
    if (As == NULL)  return NULL;

    /* copy nonzeros from A into As */
    nnz = 0;
    for (j=0; j<A->N; j++) {
      As->colptrs[j] = nnz;
      for (i=j-(A->mu); i<j+(A->ml); i++) {
	dtmp = BAND_ELEM(A,i,j);
	if ( dtmp != 0.0 ) { 
	  As->rowvals[nnz] = i;
	  As->data[nnz++] = dtmp;
	}
      }
    }
    As->colptrs[A->N] = nnz;

  }

  return(As);
}


/* 
Frees memory and deletes the structure for an existing sparse matrix.
*/
void DestroySparseMat(SlsMat A)
{
  if (A->data) {
    free(A->data);  
    A->data = NULL;
  }
  if (A->rowvals) {
    free(A->rowvals);
    A->rowvals = NULL;
  }
  if (A->colptrs) {
    free(A->colptrs);
    A->colptrs = NULL;
  }
  free(A); A = NULL;
}


/* 
Sets all entries in a sparse matrix to zero.
*/
void SlsSetToZero(SlsMat A)
{
  int i;

  for (i=0; i<A->NNZ; i++) {
    A->data[i] = ZERO;
    A->rowvals[i] = 0;
  }

  for (i=0; i<A->N; i++) {
    A->colptrs[i] = 0;
  }
  /* A->colptrs[A->N] = A->NNZ; */
  A->colptrs[A->N] = 0;

}


/* 
Copies the sparse matrix A into sparse matrix B.  It is assumed 
that A and B have the same dimensions, but we account for the 
situation in which B has fewer nonzeros than A.
*/
void CopySparseMat(SlsMat A, SlsMat B)
{
  int i;
  int A_nz = A->colptrs[A->N];

  /* ensure that B is allocated with at least as 
     much memory as we have nonzeros in A */
  if (B->NNZ < A_nz) {
    B->rowvals = realloc(B->rowvals, A_nz*sizeof(int));
    B->data = realloc(B->data, A_nz*sizeof(realtype));
    B->NNZ = A_nz;
  }

  /* zero out B so that copy works correctly */
  SlsSetToZero(B);

  /* copy the data and row indices over */
  for (i=0; i<A_nz; i++){
    B->data[i] = A->data[i];
    B->rowvals[i] = A->rowvals[i];
  }

  /* copy the column pointers over */
  for (i=0; i<A->N; i++) {
    B->colptrs[i] = A->colptrs[i];
  }
  B->colptrs[A->N] = A_nz;

}


/* 
Scales a sparse matrix A by the coefficient b.
*/
void ScaleSparseMat(realtype b, SlsMat A)
{
  int i;

  for (i=0; i<A->colptrs[A->N]; i++){
    A->data[i] = b * (A->data[i]);
  }
}


/* 
Adds the identity to a sparse matrix.  Resizes A if necessary 
(if it had 0-valued diagonal entries).
*/
/* void AddIdentitySparseMat(SlsMat A) */
/* { */
/*   int j, i, p, M, N, nz; */
/*   int *w, *Cp, *Ap, *Ai, *Ci; */
/*   realtype *x, *Ax, *Cx; */
/*   SlsMat C; */

/*   M = A->M; */
/*   N = A->N; */

/*   w = (int *)malloc(M * sizeof(int)); */
/*   x = (realtype *)malloc(M * sizeof(realtype)); */
/*   C = NewSparseMat(A->M, A->N, (A->NNZ)+M); */

/*   Cp = C->colptrs; */
/*   Ci = C->rowvals; */
/*   Cx = C->data; */
/*   Ap = A->colptrs; */
/*   Ai = A->rowvals; */
/*   Ax = A->data; */

/*   /\* Initialize values *\/ */
/*   nz = 0; */
/*   for (j=0; j<M; j++) { */
/*     w[j] = 0; */
/*     x[j] = 0.0; */
/*   } */

/*   for (j=0; j<N; j++) { */
/*     Cp[j] = nz; */
/*     for (p=Ap[j]; p<Ap[j+1]; p++) { */
/*       i = Ai[p]; */
/*       w[i] = j+1; */
/*       Ci[nz] = i; */
/*       nz++; */
/*       x[i] = Ax[p]; */
/*     } */
/*     if (w[j] < j+1) { */
/*       Ci[nz] = j; */
/*       nz++; */
/*       x[j] = 1.0; */
/*     } else { */
/*       x[j] += 1.0; */
/*     } */
/*     for (p=Cp[j]; p<nz; p++) { */
/*       Cx[p] = x[Ci[p]]; */
/*     } */

/*   } */
/*   Cp[N] = nz; */
  
/*   A->N = C->N; */
/*   A->M = C->M; */
/*   A->NNZ = C->NNZ; */

/*   if (A->data) { */
/*     free(A->data);   */
/*     A->data = C->data; */
/*     C->data = NULL; */
/*   } */
/*   if (A->rowvals) { */
/*     free(A->rowvals); */
/*     A->rowvals = C->rowvals; */
/*     C->rowvals = NULL; */
/*   } */
/*   if (A->colptrs) { */
/*     free(A->colptrs); */
/*     A->colptrs = C->colptrs; */
/*     C->colptrs = NULL; */
/*   } */

/*   DestroySparseMat(C);  */
/*   free(w); */
/*   free(x); */

/*   /\*  Reallocate the new matrix to get rid of any extra space *\/ */
/*   ReallocSparseMat(A); */

/* } */


/* 
Adds 1 to every diagonal entry of A.  Works for general [rectangular] 
matrices and handles potentially increased size if A does not currently 
contain a value on the diagonal.
*/
void AddIdentitySparseMat(SlsMat A)
{
  int j, i, p, nz, newmat, found;
  int *w, *Ap, *Ai, *Cp, *Ci;
  realtype *x, *Ax, *Cx;
  SlsMat C;

  /* determine if A already contains values on the diagonal (hence 
     memory allocation necessary)*/
  newmat=0;
  for (j=0; j < SUNMIN(A->N,A->M); j++) {
    /* scan column of A, searching for diagonal value */
    found = 0;
    for (i=A->colptrs[j]; i<A->colptrs[j+1]; i++)
      if (A->rowvals[i] == j) {
	found = 1;
	break;
      }
    /* if no diagonal found, signal new matrix */
    if (!found) {
      newmat=1;
      break;
    }
  }

  /* perform operation */

  /*   case 1: A already contains a diagonal */
  if (!newmat) {

    /* iterate through columns, adding 1.0 to diagonal */
    for (j=0; j < SUNMIN(A->N,A->M); j++)
      for (i=A->colptrs[j]; i<A->colptrs[j+1]; i++)
	if (A->rowvals[i] == j) 
	  A->data[i] += ONE;

  /*   case 2: A does not already contain a diagonal */
  } else {

    /* create work arrays for row indices and nonzero column values */
    w = (int *) malloc(A->M * sizeof(int));
    x = (realtype *) malloc(A->M * sizeof(realtype));

    /* create new matrix for sum (overestimate nnz as sum of each) */
    C = NewSparseMat(A->M, A->N, (A->colptrs[A->N]) + SUNMIN(A->M, A->N));

    /* access data from CSR structures (return if failure) */
    Cp = Ci = Ap = Ai = NULL;
    Cx = Ax = NULL;
    if (C->colptrs)  Cp = C->colptrs;
    else  return;
    if (C->rowvals)  Ci = C->rowvals;
    else  return;
    if (C->data)     Cx = C->data;
    else  return;
    if (A->colptrs)  Ap = A->colptrs;
    else  return;
    if (A->rowvals)  Ai = A->rowvals;
    else  return;
    if (A->data)     Ax = A->data;
    else  return;

    /* initialize total nonzero count */
    nz = 0;

    /* iterate through columns */
    for (j=0; j<A->N; j++) {

      /* set current column pointer to current # nonzeros */
      Cp[j] = nz;

      /* clear out temporary arrays for this column */
      for (i=0; i<A->M; i++) {
	w[i] = 0;
	x[i] = 0.0;
      }

      /* iterate down column of A, collecting nonzeros */
      for (p=Ap[j]; p<Ap[j+1]; p++) {
	w[Ai[p]] += 1;       /* indicate that row is filled */
	x[Ai[p]] = Ax[p];    /* collect value */
      }

      /* add identity to this column */
      if (j < A->M) {
	w[j] += 1;     /* indicate that row is filled */
	x[j] += ONE;   /* update value */
      }

      /* fill entries of C with this column's data */
      for (i=0; i<A->M; i++) {
	if ( w[i] > 0 ) { 
	  Ci[nz] = i;  
	  Cx[nz++] = x[i];
	}
      }
    }

    /* indicate end of data */
    Cp[A->N] = nz;

    /* update A's structure with C's values; nullify C's pointers */
    A->NNZ = C->NNZ;

    if (A->data)
      free(A->data);  
    A->data = C->data;
    C->data = NULL;

    if (A->rowvals)
      free(A->rowvals);
    A->rowvals = C->rowvals;
    C->rowvals = NULL;

    if (A->colptrs)
      free(A->colptrs);
    A->colptrs = C->colptrs;
    C->colptrs = NULL;

    /* clean up */
    DestroySparseMat(C); 
    free(w);
    free(x);

    /* reallocate the new matrix to remove extra space */
    ReallocSparseMat(A);
  }

}


/* 
Add two sparse matrices: A = A+B.  Handles potentially increased size if 
matrices have different sparsity patterns.  Returns 0 if successful, and
1 if unsuccessful (in which case A is left unchanged).
*/
int SlsAddMat(SlsMat A, SlsMat B)
{
  int j, i, p, nz, newmat;
  int *w, *Ap, *Ai, *Bp, *Bi, *Cp, *Ci;
  realtype *x, *Ax, *Bx, *Cx;
  SlsMat C;

  /* ensure that matrix dimensions agree */
  if ((A->M != B->M) || (A->N != B->N))
    return(-1);

  /* create work arrays for row indices and nonzero column values */
  w = (int *) malloc(A->M * sizeof(int));
  x = (realtype *) malloc(A->M * sizeof(realtype));

  /* determine if A already contains the sparsity pattern of B */
  newmat=0;
  for (j=0; j<A->N; j++) {

    /* clear work array */
    for (i=0; i<A->M; i++)  w[i] = 0;

    /* scan column of A, incrementing w by one */
    for (i=A->colptrs[j]; i<A->colptrs[j+1]; i++)
      w[A->rowvals[i]] += 1;

    /* scan column of B, decrementing w by one */
    for (i=B->colptrs[j]; i<B->colptrs[j+1]; i++)
      w[B->rowvals[i]] -= 1;

    /* if any entry of w is negative, A doesn't contain B's sparsity */
    for (i=0; i<A->M; i++)
      if (w[i] < 0) {
	newmat = 1;
	break;
      }
    if (newmat) break;

  }

  /* perform operation */

  /*   case 1: A already contains sparsity pattern of B */
  if (!newmat) {

    /* iterate through columns, adding matrices */
    for (j=0; j<A->N; j++) {

      /* clear work array */
      for (i=0; i<A->M; i++)
	x[i] = ZERO;

      /* scan column of B, updating work array */
      for (i=B->colptrs[j]; i<B->colptrs[j+1]; i++)
	x[B->rowvals[i]] = B->data[i];

      /* scan column of A, updating entries appropriately array */
      for (i=A->colptrs[j]; i<A->colptrs[j+1]; i++)
	A->data[i] += x[A->rowvals[i]];

    }

  /*   case 2: A does not already contain B's sparsity */
  } else {

    /* create new matrix for sum (overestimate nnz as sum of each) */
    C = NewSparseMat(A->M, A->N, (A->colptrs[A->N])+(B->colptrs[B->N]));

    /* access data from CSR structures (return if failure) */
    Cp = Ci = Ap = Ai = Bp = Bi = NULL;
    Cx = Ax = Bx = NULL;
    if (C->colptrs)  Cp = C->colptrs;
    else  return(-1);
    if (C->rowvals)  Ci = C->rowvals;
    else  return(-1);
    if (C->data)     Cx = C->data;
    else  return(-1);
    if (A->colptrs)  Ap = A->colptrs;
    else  return(-1);
    if (A->rowvals)  Ai = A->rowvals;
    else  return(-1);
    if (A->data)     Ax = A->data;
    else  return(-1);
    if (B->colptrs)  Bp = B->colptrs;
    else  return(-1);
    if (B->rowvals)  Bi = B->rowvals;
    else  return(-1);
    if (B->data)     Bx = B->data;
    else  return(-1);

    /* initialize total nonzero count */
    nz = 0;

    /* iterate through columns */
    for (j=0; j<C->N; j++) {

      /* set current column pointer to current # nonzeros */
      Cp[j] = nz;

      /* clear out temporary arrays for this column */
      for (i=0; i<C->M; i++) {
	w[i] = 0;
	x[i] = 0.0;
      }

      /* iterate down column of A, collecting nonzeros */
      for (p=Ap[j]; p<Ap[j+1]; p++) {
	w[Ai[p]] += 1;       /* indicate that row is filled */
	x[Ai[p]] = Ax[p];    /* collect value */
      }

      /* iterate down column of B, collecting nonzeros */
      for (p=Bp[j]; p<Bp[j+1]; p++) {
	w[Bi[p]] += 1;       /* indicate that row is filled */
	x[Bi[p]] += Bx[p];   /* collect value */
      }

      /* fill entries of C with this column's data */
      for (i=0; i<C->M; i++) {
	if ( w[i] > 0 ) { 
	  Ci[nz] = i;  
	  Cx[nz++] = x[i];
	}
      }
    }

    /* indicate end of data */
    Cp[A->N] = nz;

    /* update A's structure with C's values; nullify C's pointers */
    A->NNZ = C->NNZ;

    free(A->data);  
    A->data = C->data;
    C->data = NULL;

    free(A->rowvals);
    A->rowvals = C->rowvals;
    C->rowvals = NULL;

    free(A->colptrs);
    A->colptrs = C->colptrs;
    C->colptrs = NULL;

    /* clean up */
    DestroySparseMat(C); 

    /* reallocate the new matrix to remove extra space */
    ReallocSparseMat(A);

  }

  /* clean up */
  free(w);
  free(x);

  /* return success */
  return(0);
}


/* 
Resizes the memory allocated for a given sparse matrix, shortening 
it down to the number of actual nonzero entries.
*/
void ReallocSparseMat(SlsMat A)
{
  int nzmax; 

  nzmax = A->colptrs[A->N];
  A->rowvals = realloc(A->rowvals, nzmax*sizeof(int));
  A->data = realloc(A->data, nzmax*sizeof(realtype));
  A->NNZ = nzmax;
}


/* 
Computes y=A*x, where A is a sparse matrix of dimension MxN, x is a 
realtype array of length N, and y is a realtype array of length M. 
Returns 0 if successful, and 1 if unsuccessful (failed memory access).
*/
int SlsMatvec(SlsMat A, realtype *x, realtype *y)
{
  int j, i;
  int *Ap, *Ai;
  realtype *Ax;

  /* access data from CSR structure (return if failure) */
  if (A->colptrs)  Ap = A->colptrs;
  else  return(-1);
  if (A->rowvals)  Ai = A->rowvals;
  else  return(-1);
  if (A->data)     Ax = A->data;
  else  return(-1);

  /* ensure that vectors are non-NULL */
  if ((x == NULL) || (y == NULL))
    return(-1);

  /* initialize result */
  for (i=0; i<A->M; i++)
    y[i] = 0.0;

  /* iterate through matrix columns */
  for (j=0; j<A->N; j++) {

    /* iterate down column of A, performing product */
    for (i=Ap[j]; i<Ap[j+1]; i++)
      y[Ai[i]] += Ax[i]*x[j];

  }

  /* return success */
  return(0);
}


/* 
Prints the nonzero entries of a sparse matrix to screen.
*/
void PrintSparseMat(SlsMat A)
{
  int i,j, M, N, NNZ;
  int *colptrs;

  colptrs = A->colptrs;

  M = A->M;
  N = A->N;
  NNZ = A->NNZ;

  printf("\n");
  
  printf("%d by %d NNZ: %d \n", M, N, NNZ);
  for (j=0; j < A->N; j++) {
    printf("  col %d : locations %d to %d\n", j, colptrs[j], colptrs[j+1]-1);
    for (i = colptrs[j]; i < colptrs[j+1]; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%d  %12Lg  ", A->rowvals[i], A->data[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%d  %12g  ", A->rowvals[i], A->data[i]);
#else
      printf("%d  %12g  ", A->rowvals[i], A->data[i]);
#endif
    }
    printf("\n");
  }
  printf("\n");
    
}

