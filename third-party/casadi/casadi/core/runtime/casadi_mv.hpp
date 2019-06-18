// NOLINT(legal/copyright)
// SYMBOL "mv"
template<typename T1>
void casadi_mv(const T1* x, const casadi_int* sp_x, const T1* y, T1* z, casadi_int tr) {
  casadi_int ncol_x, i, el;
  const casadi_int *colind_x, *row_x;
  if (!x || !y || !z) return;
  // Get sparsities
  ncol_x = sp_x[1];
  colind_x = sp_x+2; row_x = sp_x + 2 + ncol_x+1;
  if (tr) {
    // loop over the columns of x
    for (i=0; i<ncol_x; ++i) {
      // loop over the non-zeros of x
      for (el=colind_x[i]; el<colind_x[i+1]; ++el) {
        z[i] += x[el] * y[row_x[el]];
      }
    }
  } else {
    // loop over the columns of x
    for (i=0; i<ncol_x; ++i) {
      // loop over the non-zeros of x
      for (el=colind_x[i]; el<colind_x[i+1]; ++el) {
        z[row_x[el]] += x[el] * y[i];
      }
    }
  }
}
