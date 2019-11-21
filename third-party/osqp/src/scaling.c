#include "scaling.h"

#if EMBEDDED != 1


// Set values lower than threshold SCALING_REG to 1
void limit_scaling(c_float *D, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    D[i] = D[i] < MIN_SCALING ? 1.0 : D[i];
    D[i] = D[i] > MAX_SCALING ? MAX_SCALING : D[i];
  }
}

/**
 * Compute infinite norm of the colums of the KKT matrix without forming it
 *
 * The norm is stored in the vector v = (D, E)
 *
 * @param P        Cost matrix
 * @param A        Contraints matrix
 * @param D        Norm of columns related to variables
 * @param D_temp_A Temporary vector for norm of columns of A
 * @param E        Norm of columns related to constraints
 * @param n        Dimension of KKT matrix
 */
void compute_inf_norm_cols_KKT(const csc *P, const csc *A,
                               c_float *D, c_float *D_temp_A,
                               c_float *E, c_int n) {
  // First half
  //  [ P ]
  //  [ A ]
  mat_inf_norm_cols_sym_triu(P, D);
  mat_inf_norm_cols(A, D_temp_A);
  vec_ew_max_vec(D, D_temp_A, D, n);

  // Second half
  //  [ A']
  //  [ 0 ]
  mat_inf_norm_rows(A, E);
}

c_int scale_data(OSQPWorkspace *work) {
  // Scale KKT matrix
  //
  //    [ P   A']
  //    [ A   0 ]
  //
  // with diagonal matrix
  //
  //  S = [ D    ]
  //      [    E ]
  //


  c_int   i;          // Iterations index
  c_int   n, m;       // Number of constraints and variables
  c_float c_temp;     // Cost function scaling
  c_float inf_norm_q; // Infinity norm of q

  n = work->data->n;
  m = work->data->m;

  // Initialize scaling to 1
  work->scaling->c = 1.0;
  vec_set_scalar(work->scaling->D,    1., work->data->n);
  vec_set_scalar(work->scaling->Dinv, 1., work->data->n);
  vec_set_scalar(work->scaling->E,    1., work->data->m);
  vec_set_scalar(work->scaling->Einv, 1., work->data->m);


  for (i = 0; i < work->settings->scaling; i++) {
    //
    // First Ruiz step
    //

    // Compute norm of KKT columns
    compute_inf_norm_cols_KKT(work->data->P, work->data->A,
                              work->D_temp, work->D_temp_A,
                              work->E_temp, n);

    // Set to 1 values with 0 norms (avoid crazy scaling)
    limit_scaling(work->D_temp, n);
    limit_scaling(work->E_temp, m);

    // Take square root of norms
    vec_ew_sqrt(work->D_temp, n);
    vec_ew_sqrt(work->E_temp, m);

    // Divide scalings D and E by themselves
    vec_ew_recipr(work->D_temp, work->D_temp, n);
    vec_ew_recipr(work->E_temp, work->E_temp, m);

    // Equilibrate matrices P and A and vector q
    // P <- DPD
    mat_premult_diag(work->data->P, work->D_temp);
    mat_postmult_diag(work->data->P, work->D_temp);

    // A <- EAD
    mat_premult_diag(work->data->A, work->E_temp);
    mat_postmult_diag(work->data->A, work->D_temp);

    // q <- Dq
    vec_ew_prod(work->D_temp,     work->data->q, work->data->q,    n);

    // Update equilibration matrices D and E
    vec_ew_prod(work->scaling->D, work->D_temp,  work->scaling->D, n);
    vec_ew_prod(work->scaling->E, work->E_temp,  work->scaling->E, m);

    //
    // Cost normalization step
    //

    // Compute avg norm of cols of P
    mat_inf_norm_cols_sym_triu(work->data->P, work->D_temp);
    c_temp = vec_mean(work->D_temp, n);

    // Compute inf norm of q
    inf_norm_q = vec_norm_inf(work->data->q, n);

    // If norm_q == 0, set it to 1 (ignore it in the scaling)
    // NB: Using the same function as with vectors here
    limit_scaling(&inf_norm_q, 1);

    // Compute max between avg norm of cols of P and inf norm of q
    c_temp = c_max(c_temp, inf_norm_q);

    // Limit scaling (use same function as with vectors)
    limit_scaling(&c_temp, 1);

    // Invert scaling c = 1 / cost_measure
    c_temp = 1. / c_temp;

    // Scale P
    mat_mult_scalar(work->data->P, c_temp);

    // Scale q
    vec_mult_scalar(work->data->q, c_temp, n);

    // Update cost scaling
    work->scaling->c *= c_temp;
  }


  // Store cinv, Dinv, Einv
  work->scaling->cinv = 1. / work->scaling->c;
  vec_ew_recipr(work->scaling->D, work->scaling->Dinv, work->data->n);
  vec_ew_recipr(work->scaling->E, work->scaling->Einv, work->data->m);


  // Scale problem vectors l, u
  vec_ew_prod(work->scaling->E, work->data->l, work->data->l, work->data->m);
  vec_ew_prod(work->scaling->E, work->data->u, work->data->u, work->data->m);

  return 0;
}

#endif // EMBEDDED

c_int unscale_data(OSQPWorkspace *work) {
  // Unscale cost
  mat_mult_scalar(work->data->P, work->scaling->cinv);
  mat_premult_diag(work->data->P, work->scaling->Dinv);
  mat_postmult_diag(work->data->P, work->scaling->Dinv);
  vec_mult_scalar(work->data->q, work->scaling->cinv, work->data->n);
  vec_ew_prod(work->scaling->Dinv, work->data->q, work->data->q, work->data->n);

  // Unscale constraints
  mat_premult_diag(work->data->A, work->scaling->Einv);
  mat_postmult_diag(work->data->A, work->scaling->Dinv);
  vec_ew_prod(work->scaling->Einv, work->data->l, work->data->l, work->data->m);
  vec_ew_prod(work->scaling->Einv, work->data->u, work->data->u, work->data->m);

  return 0;
}

// // Scale solution
// c_int scale_solution(OSQPWorkspace * work){
//
//     // primal
//     vec_ew_prod(work->scaling->Dinv, work->solution->x, work->data->n);
//
//     // dual
//     vec_ew_prod(work->scaling->Einv, work->solution->y, work->data->m);
//     vec_mult_scalar(work->solution->y, work->scaling->c, work->data->m);
//
//     return 0;
// }


c_int unscale_solution(OSQPWorkspace *work) {
  // primal
  vec_ew_prod(work->scaling->D,
              work->solution->x,
              work->solution->x,
              work->data->n);

  // dual
  vec_ew_prod(work->scaling->E,
              work->solution->y,
              work->solution->y,
              work->data->m);
  vec_mult_scalar(work->solution->y, work->scaling->cinv, work->data->m);

  return 0;
}
