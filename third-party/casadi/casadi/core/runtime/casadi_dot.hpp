// NOLINT(legal/copyright)
// SYMBOL "dot"
template<typename T1>
T1 casadi_dot(casadi_int n, const T1* x, const T1* y) {
  casadi_int i;
  T1 r = 0;
  for (i=0; i<n; ++i) r += *x++ * *y++;
  return r;
}
