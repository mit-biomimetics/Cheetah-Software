// NOLINT(legal/copyright)
// SYMBOL "bilin"
template<typename T1>
T1 casadi_bilin(const T1* A, const casadi_int* sp_A, const T1* x, const T1* y) {
  casadi_int ncol_A, cc, rr, el;
  const casadi_int *colind_A, *row_A;
  T1 ret;
  // Get sparsities
  ncol_A = sp_A[1];
  colind_A = sp_A+2; row_A = sp_A + 2 + ncol_A+1;
  // Return value
  ret=0;
  // Loop over the columns of A
  for (cc=0; cc<ncol_A; ++cc) {
    // Loop over the nonzeros of A
    for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
      // Get the row
      rr=row_A[el];
      // Add contribution
      ret += x[rr]*A[el]*y[cc];
    }
  }
  return ret;
}
