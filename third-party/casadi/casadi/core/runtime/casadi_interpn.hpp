// NOLINT(legal/copyright)
// SYMBOL "interpn"
template<typename T1>
void casadi_interpn(T1* res, casadi_int ndim, const T1* grid, const casadi_int* offset, const T1* values, const T1* x, const casadi_int* lookup_mode, casadi_int m, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  // Work vectors
  T1* alpha;
  casadi_int *index, *corner;
  alpha = w; w += ndim;
  index = iw; iw += ndim;
  corner = iw; iw += ndim;
  // Left index and fraction of interval
  casadi_interpn_weights(ndim, grid, offset, x, alpha, index, lookup_mode);
  // Loop over all corners, add contribution to output
  casadi_fill_casadi_int(corner, ndim, 0);
  casadi_fill(res, m, 0.0);
  do {
    T1* coeff = 0;
    casadi_interpn_interpolate(res, ndim, offset, values,
      alpha, index, corner, coeff, m);
  } while (casadi_flip(corner, ndim));
}
