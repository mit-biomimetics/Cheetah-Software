// NOLINT(legal/copyright)
// SYMBOL "norm_1"
template<typename T1>
T1 casadi_norm_1(casadi_int n, const T1* x) {
  casadi_int i;
  T1 ret = 0;
  if (x) {
    for (i=0; i<n; ++i) ret += fabs(*x++);
  }
  return ret;
}
