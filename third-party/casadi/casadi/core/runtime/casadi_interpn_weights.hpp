// NOLINT(legal/copyright)
// SYMBOL "interpn_weights"
template<typename T1>
void casadi_interpn_weights(casadi_int ndim, const T1* grid, const casadi_int* offset, const T1* x, T1* alpha, casadi_int* index, const casadi_int* lookup_mode) { // NOLINT(whitespace/line_length)
  // Left index and fraction of interval
  casadi_int i;
  for (i=0; i<ndim; ++i) {
    casadi_int ng, j;
    T1 xi;
    const T1* g;
    // Grid point
    xi = x ? x[i] : 0;
    // Grid
    g = grid + offset[i];
    ng = offset[i+1]-offset[i];
    // Find left index
    j = index[i] = casadi_low(xi, g, ng, lookup_mode[i]);
    // Get interpolation/extrapolation alpha
    alpha[i] = (xi-g[j])/(g[j+1]-g[j]);
  }
}
