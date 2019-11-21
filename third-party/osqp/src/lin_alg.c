#include "lin_alg.h"


/* VECTOR FUNCTIONS ----------------------------------------------------------*/


void vec_add_scaled(c_float       *c,
                    const c_float *a,
                    const c_float *b,
                    c_int          n,
                    c_float        sc) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] =  a[i] + sc * b[i];
  }
}

c_float vec_scaled_norm_inf(const c_float *S, const c_float *v, c_int l) {
  c_int   i;
  c_float abs_Sv_i;
  c_float max = 0.0;

  for (i = 0; i < l; i++) {
    abs_Sv_i = c_absval(S[i] * v[i]);

    if (abs_Sv_i > max) max = abs_Sv_i;
  }
  return max;
}

c_float vec_norm_inf(const c_float *v, c_int l) {
  c_int   i;
  c_float abs_v_i;
  c_float max = 0.0;

  for (i = 0; i < l; i++) {
    abs_v_i = c_absval(v[i]);

    if (abs_v_i > max) max = abs_v_i;
  }
  return max;
}

c_float vec_norm_inf_diff(const c_float *a, const c_float *b, c_int l) {
  c_float nmDiff = 0.0, tmp;
  c_int   i;

  for (i = 0; i < l; i++) {
    tmp = c_absval(a[i] - b[i]);

    if (tmp > nmDiff) nmDiff = tmp;
  }
  return nmDiff;
}

c_float vec_mean(const c_float *a, c_int n) {
  c_float mean = 0.0;
  c_int   i;

  for (i = 0; i < n; i++) {
    mean += a[i];
  }
  mean /= (c_float)n;

  return mean;
}

void int_vec_set_scalar(c_int *a, c_int sc, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] = sc;
  }
}

void vec_set_scalar(c_float *a, c_float sc, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] = sc;
  }
}

void vec_add_scalar(c_float *a, c_float sc, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] += sc;
  }
}

void vec_mult_scalar(c_float *a, c_float sc, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] *= sc;
  }
}

#ifndef EMBEDDED
c_float* vec_copy(c_float *a, c_int n) {
  c_float *b;
  c_int    i;

  b = c_malloc(n * sizeof(c_float));

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }

  return b;
}

#endif // end EMBEDDED


void prea_int_vec_copy(const c_int *a, c_int *b, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void prea_vec_copy(const c_float *a, c_float *b, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void vec_ew_recipr(const c_float *a, c_float *b, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    b[i] = (c_float)1.0 / a[i];
  }
}

c_float vec_prod(const c_float *a, const c_float *b, c_int n) {
  c_float prod = 0.0;
  c_int   i; // Index

  for (i = 0; i < n; i++) {
    prod += a[i] * b[i];
  }

  return prod;
}

void vec_ew_prod(const c_float *a, const c_float *b, c_float *c, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] = b[i] * a[i];
  }
}

#if EMBEDDED != 1
void vec_ew_sqrt(c_float *a, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] = c_sqrt(a[i]);
  }
}

void vec_ew_max(c_float *a, c_int n, c_float max_val) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] = c_max(a[i], max_val);
  }
}

void vec_ew_min(c_float *a, c_int n, c_float min_val) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] = c_min(a[i], min_val);
  }
}

void vec_ew_max_vec(const c_float *a, const c_float *b, c_float *c, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] = c_max(a[i], b[i]);
  }
}

void vec_ew_min_vec(const c_float *a, const c_float *b, c_float *c, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] = c_min(a[i], b[i]);
  }
}

#endif // EMBEDDED != 1


/* MATRIX FUNCTIONS ----------------------------------------------------------*/

/* multiply scalar to matrix */
void mat_mult_scalar(csc *A, c_float sc) {
  c_int i, nnzA;

  nnzA = A->p[A->n];

  for (i = 0; i < nnzA; i++) {
    A->x[i] *= sc;
  }
}

void mat_premult_diag(csc *A, const c_float *d) {
  c_int j, i;

  for (j = 0; j < A->n; j++) {                // Cycle over columns
    for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row in the column
      A->x[i] *= d[A->i[i]];                  // Scale by corresponding element
                                              // of d for row i
    }
  }
}

void mat_postmult_diag(csc *A, const c_float *d) {
  c_int j, i;

  for (j = 0; j < A->n; j++) {                // Cycle over columns j
    for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row i in column j
      A->x[i] *= d[j];                        // Scale by corresponding element
                                              // of d for column j
    }
  }
}

void mat_vec(const csc *A, const c_float *x, c_float *y, c_int plus_eq) {
  c_int i, j;

  if (!plus_eq) {
    // y = 0
    for (i = 0; i < A->m; i++) {
      y[i] = 0;
    }
  }

  // if A is empty
  if (A->p[A->n] == 0) {
    return;
  }

  if (plus_eq == -1) {
    // y -=  A*x
    for (j = 0; j < A->n; j++) {
      for (i = A->p[j]; i < A->p[j + 1]; i++) {
        y[A->i[i]] -= A->x[i] * x[j];
      }
    }
  } else {
    // y +=  A*x
    for (j = 0; j < A->n; j++) {
      for (i = A->p[j]; i < A->p[j + 1]; i++) {
        y[A->i[i]] += A->x[i] * x[j];
      }
    }
  }
}

void mat_tpose_vec(const csc *A, const c_float *x, c_float *y,
                   c_int plus_eq, c_int skip_diag) {
  c_int i, j, k;

  if (!plus_eq) {
    // y = 0
    for (i = 0; i < A->n; i++) {
      y[i] = 0;
    }
  }

  // if A is empty
  if (A->p[A->n] == 0) {
    return;
  }

  if (plus_eq == -1) {
    // y -=  A*x
    if (skip_diag) {
      for (j = 0; j < A->n; j++) {
        for (k = A->p[j]; k < A->p[j + 1]; k++) {
          i     = A->i[k];
          y[j] -= i == j ? 0 : A->x[k] * x[i];
        }
      }
    } else {
      for (j = 0; j < A->n; j++) {
        for (k = A->p[j]; k < A->p[j + 1]; k++) {
          y[j] -= A->x[k] * x[A->i[k]];
        }
      }
    }
  } else {
    // y +=  A*x
    if (skip_diag) {
      for (j = 0; j < A->n; j++) {
        for (k = A->p[j]; k < A->p[j + 1]; k++) {
          i     = A->i[k];
          y[j] += i == j ? 0 : A->x[k] * x[i];
        }
      }
    } else {
      for (j = 0; j < A->n; j++) {
        for (k = A->p[j]; k < A->p[j + 1]; k++) {
          y[j] += A->x[k] * x[A->i[k]];
        }
      }
    }
  }
}

#if EMBEDDED != 1
void mat_inf_norm_cols(const csc *M, c_float *E) {
  c_int j, ptr;

  // Initialize zero max elements
  for (j = 0; j < M->n; j++) {
    E[j] = 0.;
  }

  // Compute maximum across columns
  for (j = 0; j < M->n; j++) {
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      E[j] = c_max(c_absval(M->x[ptr]), E[j]);
    }
  }
}

void mat_inf_norm_rows(const csc *M, c_float *E) {
  c_int i, j, ptr;

  // Initialize zero max elements
  for (j = 0; j < M->m; j++) {
    E[j] = 0.;
  }

  // Compute maximum across rows
  for (j = 0; j < M->n; j++) {
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      i    = M->i[ptr];
      E[i] = c_max(c_absval(M->x[ptr]), E[i]);
    }
  }
}

void mat_inf_norm_cols_sym_triu(const csc *M, c_float *E) {
  c_int   i, j, ptr;
  c_float abs_x;

  // Initialize zero max elements
  for (j = 0; j < M->n; j++) {
    E[j] = 0.;
  }

  // Compute maximum across columns
  // Note that element (i, j) contributes to
  // -> Column j (as expected in any matrices)
  // -> Column i (which is equal to row i for symmetric matrices)
  for (j = 0; j < M->n; j++) {
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      i     = M->i[ptr];
      abs_x = c_absval(M->x[ptr]);
      E[j]  = c_max(abs_x, E[j]);

      if (i != j) {
        E[i] = c_max(abs_x, E[i]);
      }
    }
  }
}

#endif /* if EMBEDDED != 1 */


c_float quad_form(const csc *P, const c_float *x) {
  c_float quad_form = 0.;
  c_int   i, j, ptr;                                // Pointers to iterate over
                                                    // matrix: (i,j) a element
                                                    // pointer

  for (j = 0; j < P->n; j++) {                      // Iterate over columns
    for (ptr = P->p[j]; ptr < P->p[j + 1]; ptr++) { // Iterate over rows
      i = P->i[ptr];                                // Row index

      if (i == j) {                                 // Diagonal element
        quad_form += (c_float).5 * P->x[ptr] * x[i] * x[i];
      }
      else if (i < j) {                             // Off-diagonal element
        quad_form += P->x[ptr] * x[i] * x[j];
      }
      else {                                        // Element in lower diagonal
                                                    // part
#ifdef PRINTING
        c_eprint("quad_form matrix is not upper triangular");
#endif /* ifdef PRINTING */
        return OSQP_NULL;
      }
    }
  }
  return quad_form;
}
