// NOLINT(legal/copyright)
// SYMBOL "iamax"
template<typename T1>
casadi_int casadi_iamax(casadi_int n, const T1* x, casadi_int inc_x) {
  casadi_int largest_index, i;
  T1 t, largest_value;
  largest_value = -1.0;
  largest_index = -1;
  for (i=0; i<n; ++i) {
    t = fabs(*x);
    x += inc_x;
    if (t>largest_value) {
      largest_value = t;
      largest_index = i;
    }
  }
  return largest_index;
}
