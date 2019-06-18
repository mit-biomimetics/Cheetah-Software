// NOLINT(legal/copyright)
// SYMBOL "scal"
template<typename T1>
void casadi_scal(casadi_int n, T1 alpha, T1* x) {
  casadi_int i;
  if (!x) return;
  for (i=0; i<n; ++i) *x++ *= alpha;
}
