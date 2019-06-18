// NOLINT(legal/copyright)
// SYMBOL "fill"
template<typename T1>
void casadi_fill(T1* x, casadi_int n, T1 alpha) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}
