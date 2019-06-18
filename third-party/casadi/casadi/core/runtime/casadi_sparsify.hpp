// NOLINT(legal/copyright)
// SYMBOL "sparsify"
template<typename T1, typename T2>
void casadi_sparsify(const T1* x, T2* y, const casadi_int* sp_y, casadi_int tr) {
  casadi_int nrow_y, ncol_y, i, el;
  const casadi_int *colind_y, *row_y;
  nrow_y = sp_y[0];
  ncol_y = sp_y[1];
  colind_y = sp_y+2; row_y = sp_y+ncol_y+3;
  if (tr) {
    for (i=0; i<ncol_y; ++i) {
      for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {
        *y++ = CASADI_CAST(T2, x[i + row_y[el]*ncol_y]);
      }
    }
  } else {
    for (i=0; i<ncol_y; ++i) {
      for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {
        *y++ = CASADI_CAST(T2, x[row_y[el] + i*nrow_y]);
      }
    }
  }
}
