// NOLINT(legal/copyright)
// SYMBOL "norm_2"
template<typename T1>
T1 casadi_norm_2(casadi_int n, const T1* x) {
  return sqrt(casadi_dot(n, x, x));
}
