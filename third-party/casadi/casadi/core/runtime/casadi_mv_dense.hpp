// NOLINT(legal/copyright)
// SYMBOL "mv_dense"
template<typename T1>
void casadi_mv_dense(const T1* x, casadi_int nrow_x, casadi_int ncol_x,
    const T1* y, T1* z, casadi_int tr) {
  casadi_int i, j;
  if (!x || !y || !z) return;
  if (tr) {
    for (i=0; i<ncol_x; ++i) {
      for (j=0; j<nrow_x; ++j) {
        z[i] += *x++ * y[j];
      }
    }
  } else {
    for (i=0; i<ncol_x; ++i) {
      for (j=0; j<nrow_x; ++j) {
        z[j] += *x++ * y[i];
      }
    }
  }
}
