// NOLINT(legal/copyright)

// C-REPLACE "std::numeric_limits<double>::infinity()" "INF"

// SYMBOL "mmin"
template<typename T1>
T1 casadi_mmin(const T1* x, casadi_int n, casadi_int is_dense) {
  casadi_int i;
  T1 r = is_dense ? std::numeric_limits<double>::infinity() : 0;
  if (!x) return r;
  for (i=0; i<n; ++i) r = fmin(r, x[i]);
  return r;
}

// SYMBOL "mmax"
template<typename T1>
T1 casadi_mmax(const T1* x, casadi_int n, casadi_int is_dense) {
  casadi_int i;
  T1 r = is_dense ? -std::numeric_limits<double>::infinity() : 0;
  if (!x) return r;
  for (i=0; i<n; ++i) r = fmax(r, x[i]);
  return r;
}
