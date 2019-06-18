// NOLINT(legal/copyright)
// SYMBOL "max_viol"
template<typename T1>
T1 casadi_max_viol(casadi_int n, const T1* x, const T1* lb, const T1* ub) {
  T1 r;
  casadi_int i;
  const T1 zero = 0;
  r = 0;
  for (i=0; i<n; ++i) {
    T1 x_i, lb_i, ub_i;
    x_i = x ? *x++ : zero;
    lb_i = lb ? *lb++ : zero;
    ub_i = ub ? *ub++ : zero;
    r = fmax(r, fmax(x_i-ub_i, zero));
    r = fmax(r, fmax(lb_i-x_i, zero));
  }
  return r;
}
