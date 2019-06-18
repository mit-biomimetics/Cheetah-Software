// NOLINT(legal/copyright)
// SYMBOL "polyval"
template<typename T1>
T1 casadi_polyval(const T1* p, casadi_int n, T1 x) {
  casadi_int i;
  T1 r=p[0];
  for (i=1; i<=n; ++i) {
    r = r*x + p[i];
  }
  return r;
}
