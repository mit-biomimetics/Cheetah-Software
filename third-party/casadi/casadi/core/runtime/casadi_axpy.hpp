// NOLINT(legal/copyright)
// SYMBOL "axpy"
template<typename T1>
void casadi_axpy(casadi_int n, T1 alpha, const T1* x, T1* y) {
  casadi_int i;
  if (!x || !y) return;
  for (i=0; i<n; ++i) *y++ += alpha**x++;
}
