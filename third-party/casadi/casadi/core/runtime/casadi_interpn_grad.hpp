// NOLINT(legal/copyright)
// SYMBOL "interpn_grad"
template<typename T1>
void casadi_interpn_grad(T1* grad, casadi_int ndim, const T1* grid, const casadi_int* offset, const T1* values, const T1* x, const casadi_int* lookup_mode, casadi_int m, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  T1 *alpha, *coeff, *v;
  casadi_int *index, *corner;
  casadi_int i;
  // Quick return
  if (!grad) return;
  // Work vectors
  alpha = w; w += ndim;
  coeff = w; w += ndim;
  v = w; w+= m;
  index = iw; iw += ndim;
  corner = iw; iw += ndim;

  // Left index and fraction of interval
  casadi_interpn_weights(ndim, grid, offset, x, alpha, index, lookup_mode);
  // Loop over all corners, add contribution to output
  casadi_fill_casadi_int(corner, ndim, 0);
  casadi_fill(grad, ndim*m, 0.);
  do {
    casadi_int i, j;
    // Get coefficients
    casadi_fill(v, m, 0.);
    casadi_interpn_interpolate(v, ndim, offset, values,
      alpha, index, corner, coeff, m);
    // Propagate to alpha
    for (i=ndim-1; i>=0; --i) {
      if (corner[i]) {
        for (j=0; j<m; ++j) {
          grad[i*m+j] += v[j]*coeff[i];
          v[j] *= alpha[i];
        }
      } else {
        for (j=0; j<m; ++j) {
          grad[i*m+j] -= v[j]*coeff[i];
          v[j] *= 1-alpha[i];
        }
      }
    }
  } while (casadi_flip(corner, ndim));
  // Propagate to x
  for (i=0; i<ndim; ++i) {
    casadi_int k, j;
    const T1* g;
    T1 delta;
    g = grid + offset[i];
    j = index[i];
    delta =  g[j+1]-g[j];
    for (k=0;k<m;++k) grad[k] /= delta;
    grad += m;
  }
}
